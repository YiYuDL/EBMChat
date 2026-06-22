"""
Created on Sun Sep 28 21:51:07 2025

@author: user
"""


import sys
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from agents.agent import apps,AgentState
from langchain_core.messages import HumanMessage, AIMessage


def EBMChat(question: str, thread_id: str) -> str:
    """
    Input: question (str), thread_id (str)
    Logic:
    - If this is the first call for this thread_id, execute the initial fixed Pipeline;
    - Otherwise enter ReAct mode where the model decides whether/how to use tools;
    - Returns the text content of the final AIMessage.
    """

    input_state = {
        "messages": [HumanMessage(content=question)]
    }

    result_state: AgentState = apps.invoke(
        input_state,
        config={"configurable": {"thread_id": thread_id}}
    )


    msgs = result_state.get("messages", [])
    ai_msgs = [m for m in msgs if isinstance(m, AIMessage)]
    if ai_msgs:
        return ai_msgs[-1].content if isinstance(ai_msgs[-1].content, str) else str(ai_msgs[-1].content)

    if msgs:
        last = msgs[-1]
        return last.content if isinstance(last.content, str) else str(last.content)
    return ""
