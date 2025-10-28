#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 17:26:39 2025

@author: user
"""

import os
import re
from typing import Dict, Any, List, Optional
from typing_extensions import TypedDict
from langchain.tools import StructuredTool

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

from tools.evidence import evid_eval
from tools.query import query_opt
from tools.query import generate_query
from tools.search import pm_search
from tools.answer import answer_gen
from tools.chat import result_show


def _q2q_wrapper(question: str) -> str:
    return str(generate_query(question))

def _qopt_wrapper(query: str) -> str:
    return str(query_opt(query))

def _pm_wrapper(query: str) -> str:
    return str(pm_search(query))

def _evid_wrapper(search_results: str, question: str) -> str:
    return str(evid_eval(search_results, question))

def _ans_wrapper(pmid: str, question: str) -> str:
    return str(answer_gen(pmid, question))

def _show_wrapper(keyword: str, file_type: str) -> str:
    # 这里保持你工具的签名：keyword(str), file_type(str)
    # 你在场景里会把 memory['searched_result'] 替换到 file_type 里
    return str(result_show(keyword, file_type))

Question2Query = StructuredTool.from_function(
    func=_q2q_wrapper,
    name="Question2Query",
    description="Step 1/5. Input: question(str). Output: initial_query(str)."
)
Query_Opt = StructuredTool.from_function(
    func=_qopt_wrapper,
    name="Query_Opt",
    description="Step 2/5. Input: query(str). Output: optimized_query(str)."
)
PM_Search = StructuredTool.from_function(
    func=_pm_wrapper,
    name="PM_Search",
    description="Step 3/5. Input: query(str). Output: search_results(str)."
)
Evid_Eval = StructuredTool.from_function(
    func=_evid_wrapper,
    name="Evid_Eval",
    description="Step 4/5. Input: search_results(str), question(str). Output: pmid(str)."
)
Answer_Gen = StructuredTool.from_function(
    func=_ans_wrapper,
    name="Answer_Gen",
    description="Step 5/5. Input: pmid(str), question(str). Output: final_answer(str)."
)
Result_Show = StructuredTool.from_function(
    func=_show_wrapper,
    name="Result_Show",
    description="Show specific results. Input: keyword(str), file_type(str). Output: shown_results(str)."
)