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

    '''input your own email address'''
    Entrez.email = 'yourownemial@gmail.com'

    handle = Entrez.esearch(db="pubmed", term=pmid,retmax = 1)
    record = Entrez.read(handle)
    handle.close()
    idlist = record['IdList']


    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",retmode="text")
    records = Medline.parse(handle)
    records = list(records)
    handle.close()

    record = records[0]
    doi = record.get("LID", "?").split()[0]
    return doi

HEADERS = {
    "User-Agent": "pdf-link-finder/1.0 (+mailto:you@example.com)"
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

        for loc in data.get("oa_locations") or []:
            if loc.get("url_for_pdf"):
                return loc["url_for_pdf"]
    return None

def from_openalex(doi: str) -> Optional[str]:
    url = f"https://api.openalex.org/works/https://doi.org/{doi}"
    r = requests.get(url, timeout=15)
    if r.status_code == 200:
        j = r.json()

        loc = (j.get("best_oa_location") or {})
        if loc.get("pdf_url"):
            return loc["pdf_url"]

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

            if str(link.get("content-type", "")).lower() == "application/pdf" and link.get("URL"):
                return link["URL"]
    return None

def via_doi_content_negotiation(doi: str) -> Optional[str]:

    try:
        r = requests.get(f"https://doi.org/{doi}", headers={"Accept": "application/pdf", **HEADERS}, timeout=15, allow_redirects=True)
        ct = r.headers.get("Content-Type", "").lower()
        if "application/pdf" in ct and r.url.lower().endswith(".pdf"):
            return r.url
    except requests.RequestException:
        pass
    return None

def from_europe_pmc(doi: str) -> Optional[str]:

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
    Return a PDF link for a PMID by checking sources from most reliable to fallback options.

    Only returns open-access or directly reachable PDF links and does not bypass paywalls.
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
    total_timeout = 90
    poll_interval = 0.3
    stable_seconds = 1.5


    options = Options()
    options.add_argument("--headless")
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.dir", download_dir)
    options.set_preference("browser.download.useDownloadDir", True)
    options.set_preference("browser.download.manager.showWhenStarting", False)
    options.set_preference("browser.download.alwaysOpenPanel", False)
    options.set_preference("browser.helperApps.neverAsk.saveToDisk",
                           "application/pdf,application/octet-stream,application/x-pdf,application/vnd.pdf")
    options.set_preference("pdfjs.disabled", True)

    driver = None
    try:
        driver = webdriver.Firefox(options=options)
        driver.set_page_load_timeout(20)
    except Exception as e:
        print(f"浏览器启动失败: {e}")
        return False

    def hard_kill():
        try:

            try:
                driver.quit()
            except Exception:
                pass

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


            pdf_candidates = []
            for f in new_files:
                if f.startswith("."):
                    continue
                if f.lower().endswith(".pdf"):
                    p = os.path.join(download_dir, f)
                    if os.path.isfile(p):
                        pdf_candidates.append((os.path.getmtime(p), p))

            if pdf_candidates:

                candidate_pdf_path = max(pdf_candidates)[1]
                part_path = candidate_pdf_path + ".part"

                size = os.path.getsize(candidate_pdf_path)
                now = time.time()


                if last_size != size:
                    last_size = size
                    last_change_ts = now


                no_part = not os.path.exists(part_path)
                stable_enough = (last_change_ts is not None) and (now - last_change_ts >= stable_seconds)
                if size > 0 and no_part and stable_enough:
                    print(f"下载完成：{os.path.basename(candidate_pdf_path)} | 大小 {size} bytes")
                    break

            time.sleep(poll_interval)

    except Exception as e:
        print(f"下载监控异常: {e}")

    finally:

        try:
            if candidate_pdf_path and os.path.exists(candidate_pdf_path):
                part_still_exists = os.path.exists(candidate_pdf_path + ".part")
                size = os.path.getsize(candidate_pdf_path)
            else:
                part_still_exists = True
                size = 0

            if (time.time() - start) < total_timeout and not part_still_exists and size > 0:

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
