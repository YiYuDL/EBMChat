#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 15:42:41 2025

@author: user
"""

import requests
from typing import Optional, Dict, List
from Bio import Entrez, Medline
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
import os
import time
from pathlib import Path
import signal

def pmid2doi(pmid):
    pmid = int(pmid)
    email_for_unpaywall="yuyi689@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=pmid,retmax = 1)
    record = Entrez.read(handle)
    handle.close()
    idlist = record['IdList']
    # print(idlist)
    # print(record['Count'])
    
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",retmode="text")
    records = Medline.parse(handle)     
    records = list(records)
    handle.close()
    
    record = records[0]
    doi = record.get("LID", "?").split()[0]
    return doi

HEADERS = {
    "User-Agent": "pdf-link-finder/1.0 (+mailto:you@example.com)"  # Crossref 建议包含可联系邮箱
}

def from_unpaywall(doi: str, email: str) -> Optional[str]:
    url = f"https://api.unpaywall.org/v2/{doi}"
    r = requests.get(url, params={"email": email}, timeout=15)
    if r.status_code == 200:
        data = r.json()
        loc = data.get("best_oa_location") or {}
        pdf = loc.get("url_for_pdf")
        if pdf:
            return pdf
        # 退而求其次：其他 OA 位置中找 pdf
        for loc in data.get("oa_locations") or []:
            if loc.get("url_for_pdf"):
                return loc["url_for_pdf"]
    return None

def from_openalex(doi: str) -> Optional[str]:
    url = f"https://api.openalex.org/works/https://doi.org/{doi}"
    r = requests.get(url, timeout=15)
    if r.status_code == 200:
        j = r.json()
        # 先查 best_oa_location
        loc = (j.get("best_oa_location") or {})
        if loc.get("pdf_url"):
            return loc["pdf_url"]
        # 再查 primary_location / locations
        p = j.get("primary_location") or {}
        if p.get("pdf_url"):
            return p["pdf_url"]
        for loc in j.get("locations") or []:
            if loc.get("pdf_url"):
                return loc["pdf_url"]
    return None

def from_crossref_links(doi: str) -> Optional[str]:
    url = f"https://api.crossref.org/works/{doi}"
    r = requests.get(url, headers=HEADERS, timeout=15)
    if r.status_code == 200:
        msg = r.json().get("message", {})
        for link in msg.get("link", []) or []:
            # 过滤出 content-type 为 PDF 的直链
            if str(link.get("content-type", "")).lower() == "application/pdf" and link.get("URL"):
                return link["URL"]
    return None

def via_doi_content_negotiation(doi: str) -> Optional[str]:
    # 少数出版社支持对 DOI 直接做 Accept: application/pdf 的内容协商，可能 303 跳转到 PDF
    try:
        r = requests.get(f"https://doi.org/{doi}", headers={"Accept": "application/pdf", **HEADERS}, timeout=15, allow_redirects=True)
        ct = r.headers.get("Content-Type", "").lower()
        if "application/pdf" in ct and r.url.lower().endswith(".pdf"):
            return r.url
    except requests.RequestException:
        pass
    return None

def from_europe_pmc(doi: str) -> Optional[str]:
    # 命中 PMC/仓库时可有 fullTextUrlList
    api = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    r = requests.get(api, params={"query": f"DOI:{doi}", "format": "json", "resultType": "core"}, timeout=15)
    if r.status_code == 200:
        res = (r.json().get("resultList") or {}).get("result") or []
        for item in res:
            for u in (item.get("fullTextUrlList", {}) or {}).get("fullTextUrl", []) or []:
                if str(u.get("documentStyle", "")).lower() == "pdf" and u.get("url"):
                    return u["url"]
    return None

def from_semanticscholar(doi: str, api_key: Optional[str] = None) -> Optional[str]:
    # S2 Graph API: openAccessPdf.url
    url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{doi}"
    params = {"fields": "openAccessPdf"}
    headers = {"x-api-key": api_key} if api_key else {}
    r = requests.get(url, params=params, headers=headers, timeout=15)
    if r.status_code == 200:
        data = r.json()
        pdf = (data.get("openAccessPdf") or {}).get("url")
        if pdf:
            return pdf
    return None

def get_pdf_link_by_pmid(pmid) -> Dict[str, Optional[str]]:
    """
    返回一个包含来源与链接的结果；按最稳妥到最“尝试”的顺序探测。
    仅返回开放获取或明确可直达的 PDF 链接；不绕过付费墙。
    """
    doi = pmid2doi(pmid)
    email_for_unpaywall="yuyi689@gmail.com"
    checks: List = [
        ("unpaywall", lambda: from_unpaywall(doi, email_for_unpaywall)),
        ("openalex", lambda: from_openalex(doi)),
        ("crossref_links", lambda: from_crossref_links(doi)),
        ("europe_pmc", lambda: from_europe_pmc(doi)),
        ("semantic_scholar", lambda: from_semanticscholar(doi, None)),
        ("doi_content_negotiation", lambda: via_doi_content_negotiation(doi)),
    ]
    for name, fn in checks:
        url = fn()
        if url:
            return url
    return None

def download_by_link(link, pmid):
    
    parent_dir = Path.cwd().parent
    source_dir = parent_dir / "source"

    if not source_dir.exists():
        source_dir.mkdir()
        print(f"source folder: {source_dir} is created")
    else:
        print(f"source folder: {source_dir} is existed")
    
    download_dir = source_dir
    pmid = str(pmid)
    Path(download_dir).mkdir(parents=True, exist_ok=True)

    target_filename = f"{pmid}.pdf"
    target_path = os.path.join(download_dir, target_filename)
    total_timeout = 90  # 建议放宽到 60~120 秒，视网速而定
    poll_interval = 0.3
    stable_seconds = 1.5

    # Firefox 配置（确保直接下载）
    options = Options()
    options.add_argument("--headless")
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.dir", download_dir)
    options.set_preference("browser.download.useDownloadDir", True)
    options.set_preference("browser.download.manager.showWhenStarting", False)
    options.set_preference("browser.download.alwaysOpenPanel", False)
    options.set_preference("browser.helperApps.neverAsk.saveToDisk",
                           "application/pdf,application/octet-stream,application/x-pdf,application/vnd.pdf")
    options.set_preference("pdfjs.disabled", True)  # 禁用内置 PDF 预览

    driver = None
    try:
        driver = webdriver.Firefox(options=options)
        driver.set_page_load_timeout(20)
    except Exception as e:
        print(f"浏览器启动失败: {e}")
        return False

    def hard_kill():
        try:
            # 先尝试正常退出
            try:
                driver.quit()
            except Exception:
                pass
            # 再兜底杀进程
            if hasattr(driver, "service") and hasattr(driver.service, "process") and driver.service.process:
                os.kill(driver.service.process.pid, signal.SIGKILL)
        except Exception:
            pass

    start = time.time()
    initial = set(os.listdir(download_dir))
    candidate_pdf_path = None
    last_size = None
    last_change_ts = None

    try:
        print("正在发起PDF下载请求...")
        driver.get(link)

        while time.time() - start < total_timeout:
            current = set(os.listdir(download_dir))
            new_files = sorted(current - initial)

            # 在所有新出现的 .pdf 中选择“最新修改”的作为候选
            pdf_candidates = []
            for f in new_files:
                if f.startswith("."):
                    continue
                if f.lower().endswith(".pdf"):
                    p = os.path.join(download_dir, f)
                    if os.path.isfile(p):
                        pdf_candidates.append((os.path.getmtime(p), p))

            if pdf_candidates:
                # 最新的那个
                candidate_pdf_path = max(pdf_candidates)[1]
                part_path = candidate_pdf_path + ".part"

                size = os.path.getsize(candidate_pdf_path)
                now = time.time()

                # 更新“文件大小变化”的计时
                if last_size != size:
                    last_size = size
                    last_change_ts = now

                # 完成判定：无 .part，且大小稳定 stable_seconds，且大于 0
                no_part = not os.path.exists(part_path)
                stable_enough = (last_change_ts is not None) and (now - last_change_ts >= stable_seconds)
                if size > 0 and no_part and stable_enough:
                    print(f"下载完成：{os.path.basename(candidate_pdf_path)} | 大小 {size} bytes")
                    break

            time.sleep(poll_interval)

    except Exception as e:
        print(f"下载监控异常: {e}")

    finally:
        # 只有在确认为“完成”后再处理改名；否则直接退出
        try:
            if candidate_pdf_path and os.path.exists(candidate_pdf_path):
                part_still_exists = os.path.exists(candidate_pdf_path + ".part")
                size = os.path.getsize(candidate_pdf_path)
            else:
                part_still_exists = True  # 没找到候选，视为未完成
                size = 0

            if (time.time() - start) < total_timeout and not part_still_exists and size > 0:
                # 安全重命名：如果目标名与源名不同再执行
                if os.path.abspath(candidate_pdf_path) != os.path.abspath(target_path):
                    try:
                        os.replace(candidate_pdf_path, target_path)
                        print(f"重命名为: {target_filename}")
                    except Exception as e:
                        print(f"重命名失败: {e}")
                        hard_kill()
                        return False
                else:
                    print(f"文件已是目标名: {target_filename}")
                hard_kill()
                return True
            else:
                print("未检测到完整PDF（可能仍在下载或超时）")
                hard_kill()
                return False
        except Exception as e:
            print(f"收尾阶段异常: {e}")
            hard_kill()
            return False



