#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 17:03:54 2025

@author: user
"""

#shi
import sys
import os
from pathlib import Path

# 获取项目根目录（当前文件的父目录的父目录）
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# 添加到系统路径
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
#shi

from pathlib import Path
from typing import Union
from utils.download import get_pdf_link_by_pmid,download_by_link
from utils.pmid2citation import pmid2citation
from utils.answer import RAG
# import os

# 假设以下两个函数已在你的项目中实现/可导入：
# def RAG(pdf_path: str, question: str) -> str: ...
# def download(pmid: str) -> Union[bool, str, Path, None]: ...
# - download 返回值容错说明：
#   True/False 表示下载成功/失败；
#   str/Path 表示下载后文件或目录路径；
#   None/其他 表示未知，需要通过文件是否出现来判断。

# SOURCE_DIR = Path('/home/user/source')
# if not os.path.exists(SOURCE_DIR):
#     os.makedirs(SOURCE_DIR, exist_ok=True)

def answer_gen(pmid: str, question: str) -> str:
    """
    根据 PMID 与问题生成答案：
    1) 若 source/PMID.pdf 存在，直接调用 RAG(PMID.pdf, question)
    2) 若不存在，则尝试 download(PMID)
    3) 下载成功则再次调用 RAG(PMID.pdf, question)
    4) 下载仍不成功则返回提示信息
    """
    # pmid = citation.split(',')[-1].strip().rstrip('.')
    citation = pmid2citation(pmid)
    
    pmid_str = str(pmid).strip()
    pdf_name = f"{pmid_str}.pdf"
    SOURCE_DIR = Path("/home/user/source").mkdir(parents=True, exist_ok=True) or Path("/home/user/source")
    pdf_path = SOURCE_DIR / pdf_name
    
    # 情况 1：本地已有
    if pdf_path.exists():
        output = RAG(pmid, question)
        # return output
        print('ben di you')
        print('\n')
        output = output + '\n\n' + 'Citation: ' + citation
        print(output)
        return output
        
    else:
    # 情况 2：尝试下载
        try:
            link = get_pdf_link_by_pmid(pmid)
            dl_result = download_by_link(link,pmid)
            # dl_result = download(pmid_str)
        except Exception:
            dl_result = None
    
        # 判断下载结果是否成功，并定位最终 PDF 路径
        success = False
        final_pdf_path = pdf_path
    
        if isinstance(dl_result, bool):
            success = dl_result
        elif isinstance(dl_result, (str, Path)):
            p = Path(dl_result)
            # 若返回的是目录，则拼上目标文件名；若是文件，直接使用
            if p.is_dir():
                candidate = p / pdf_name
            else:
                candidate = p
            if candidate.suffix.lower() == ".pdf" and candidate.exists():
                final_pdf_path = candidate
                success = True
    
        # 即便 download 返回不明，也再检查目标位置是否已出现文件
        if not success and pdf_path.exists():
            final_pdf_path = pdf_path
            success = True

    # 情况 3：下载成功 -> 执行 RAG
        if success:
            # print('xia zai cheng gong')
            try:
                
                output = RAG(pmid, question)
                output = output + '\n\n' + 'Citation: ' + citation
                print('\n')
                print(output)
                
                return output
                
            except:
                print(f'please check whether the downloaded file is named as {pdf_name}')
                
                return None
            
        else:
            print(f'please download the literature with PMID: {pmid} manually, then use Answer_Gen tool to answer the question')
            
            return None
    # 情况 4：下载失败 -> 返回提示
    # print('please download the literature with PMID manually, then use Answer_Gen tool to answer the question')
    # return "please download the literature with PMID manually, then use Answer_Gen tool to answer the question"

# Question = 'what are recommended treatments for metastatic colorectal cancer'
# Citation = 'Metastatic colorectal cancer: ESMO Clinical Practice Guideline for diagnosis, treatment and follow-up. Guideline, Ann Oncol, 2023, https://doi.org/10.1016/j.annonc.2022.10.003, 36307056.'
# Citation1 = 'Prostate Cancer, Version 4.2023, NCCN Clinical Practice Guidelines in Oncology. Guideline, J Natl Compr Canc Netw, 2023, https://doi.org/10.6004/jnccn.2023.0050, 37856213.'
# Question1 = 'What are recommended treatments of metastatic prostate cancer'

# list1 = [39846222,39738041,37856213,36307056]
# Answer_Gen(list1[0], Question1)