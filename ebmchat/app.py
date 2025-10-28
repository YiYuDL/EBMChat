#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 28 21:51:07 2025

@author: user
"""

#shi
import sys
import os
from pathlib import Path

# 获取项目根目录（当前文件的父目录的父目录）
PROJECT_ROOT = Path(__file__).resolve().parent

# 添加到系统路径
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
#shi

from agents.agent import apps,AgentState
from langchain_core.messages import HumanMessage, AIMessage
# from agent import AgentState


def EBMChat(question: str, thread_id: str) -> str:
    """
    Input: question (str), thread_id (str)
    Logic:
    - If this is the first call for this thread_id, execute the initial fixed Pipeline;
    - Otherwise enter ReAct mode where the model decides whether/how to use tools;
    - Returns the text content of the final AIMessage.
    """
    # Wrap the current user question as a HumanMessage for input
    input_state = {
        "messages": [HumanMessage(content=question)]
    }
    # Use persistent memory via thread_id
    result_state: AgentState = apps.invoke(
        input_state,
        config={"configurable": {"thread_id": thread_id}}
    )

    # Extract the last AIMessage as output
    msgs = result_state.get("messages", [])
    ai_msgs = [m for m in msgs if isinstance(m, AIMessage)]
    if ai_msgs:
        return ai_msgs[-1].content if isinstance(ai_msgs[-1].content, str) else str(ai_msgs[-1].content)
    # If no AIMessage found, try returning the last message's text
    if msgs:
        last = msgs[-1]
        return last.content if isinstance(last.content, str) else str(last.content)
    return ""

#EBMChat("What are treatments for advanced hepatocellular carcinoma", "thread_001")
#EBMChat("please show me the initial query and optimized query from the last question", "thread_001")
#EBMChat("please show me the Guideline information from the generated table file", "thread_001")