"""
Created on Sun Aug 28 22:48:57 2025

@author: user
"""

import os
import re
from typing import Dict, Any, List, Optional
from typing_extensions import TypedDict


from langchain_openai import ChatOpenAI
from langchain.tools import StructuredTool
from langchain_core.messages import HumanMessage, AIMessage, SystemMessage, ToolMessage, AnyMessage
from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver
from agents.make_tools import Question2Query, Query_Opt, PM_Search, Evid_Eval, Answer_Gen, Result_Show


TOOLS = {
    "Question2Query": Question2Query,
    "Query_Opt": Query_Opt,
    "PM_Search": PM_Search,
    "Evid_Eval": Evid_Eval,
    "Answer_Gen": Answer_Gen,
    "Result_Show": Result_Show,
}


TOOL_OUTPUT_TO_MEMORY_KEY = {
    "Question2Query": "initial_query",
    "Query_Opt": "optimized_query",
    "PM_Search": "searched_result",
    "Evid_Eval": "pmid",
    "Answer_Gen": "final_answer",
    "Result_Show": "shown_results"
}


class AgentState(TypedDict, total=False):

    messages: List[AnyMessage]

    has_pipeline_run: bool

    memory: Dict[str, Any]

    last_tool_inputs: Dict[str, Any]
    last_tool_outputs: Dict[str, Any]


llm = ChatOpenAI(model="gpt-4.1", temperature=0.2)


REACT_SYSTEM_PROMPT = """You are an assistant with tools. For a new thread, the system will already run the 5-step pipeline automatically. For subsequent turns in the same thread, decide which tool(s) to call based on the user's request and the conversation memory.

Available tools:
1) Question2Query(question:str) -> initial_query:str
2) Query_Opt(query:str) -> optimized_query:str
3) PM_Search(query:str) -> search_results:str
4) Evid_Eval(search_results:str, question:str) -> pmid:str
5) Answer_Gen(pmid:str, question:str) -> final_answer:str
6) Result_Show(keyword:str, file_type:str) -> shown_results:str

Memory variables are available. If you need to pass a previous variable to a tool, reference memory items using a dollar-sign syntax:
- $question_last        -> last question used in the pipeline
- $initial_query        -> result from Question2Query
- $optimized_query      -> result from Query_Opt
- $searched_result      -> result from PM_Search
- $pmid                 -> result from Evid_Eval
- $final_answer         -> result from Answer_Gen

For example, if the user asks: "please show me the searched result for Guidelines in the csv file from last step",
call:
Result_Show: {"keyword": "Guidelines", "file_type": "$searched_result"}

Return a concise answer. Only show raw content when explicitly asked; otherwise summarize results. If no tool call is needed, just answer directly.
"""

def render_memory_for_prompt(mem: Dict[str, Any]) -> str:


    def short(s: Any, n: int = 1000) -> str:
        try:
            text = str(s)
        except Exception:
            text = repr(s)
        if len(text) > n:
            return text[:n] + "...[truncated]"
        return text

    parts = []
    for k in ["question_last","initial_query","optimized_query","searched_result","pmid","final_answer","shown_results"]:
        if k in mem:
            parts.append(f"{k}: {short(mem[k], 800)}")
    if not parts:
        return "No memory yet."
    return "\n".join(parts)


def route_start(state: AgentState) -> str:
    if not state.get("has_pipeline_run", False):
        return "pipeline"
    return "react_llm"


def pipeline_node(state: AgentState) -> AgentState:
    messages = state.get("messages", [])

    user_msgs = [m for m in messages if isinstance(m, HumanMessage)]
    if not user_msgs:

        return state
    question = user_msgs[-1].content if isinstance(user_msgs[-1].content, str) else str(user_msgs[-1].content)

    memory = state.get("memory", {})
    trace = []


    q2q_in = {"question": question}
    q2q_out = TOOLS["Question2Query"].invoke(q2q_in)
    trace.append({"tool": "Question2Query", "input": q2q_in, "output": q2q_out})


    qopt_in = {"query": q2q_out}
    qopt_out = TOOLS["Query_Opt"].invoke(qopt_in)
    trace.append({"tool": "Query_Opt", "input": qopt_in, "output": qopt_out})


    pm_in = {"query": qopt_out}
    pm_out = TOOLS["PM_Search"].invoke(pm_in)
    trace.append({"tool": "PM_Search", "input": pm_in, "output": pm_out})


    evid_in = {"search_results": pm_out, "question": question}
    evid_out = TOOLS["Evid_Eval"].invoke(evid_in)
    trace.append({"tool": "Evid_Eval", "input": evid_in, "output": evid_out})


    ans_in = {"pmid": evid_out, "question": question}
    ans_out = TOOLS["Answer_Gen"].invoke(ans_in)
    trace.append({"tool": "Answer_Gen", "input": ans_in, "output": ans_out})


    memory.update({
        "question_last": question,
        "initial_query": q2q_out,
        "optimized_query": qopt_out,
        "searched_result": pm_out,
        "pmid": evid_out,
        "final_answer": ans_out,
        "pipeline_trace": trace,
    })


    new_messages = messages + [AIMessage(content=ans_out)]

    return {
        **state,
        "messages": new_messages,
        "has_pipeline_run": True,
        "memory": memory,

        "last_tool_inputs": {},
        "last_tool_outputs": {},
    }


def llm_node(state: AgentState) -> AgentState:
    messages = state.get("messages", [])
    memory = state.get("memory", {})


    mem_text = render_memory_for_prompt(memory)
    system_text = REACT_SYSTEM_PROMPT + "\n\nMemory snapshot:\n" + mem_text

    llm_with_tools = llm.bind_tools(list(TOOLS.values()))

    llm_messages = [SystemMessage(content=system_text)] + messages
    ai_msg = llm_with_tools.invoke(llm_messages)

    new_messages = messages + [ai_msg]
    return {
        **state,
        "messages": new_messages
    }


def tool_exec_node(state: AgentState) -> AgentState:
    messages = state.get("messages", [])
    memory = state.get("memory", {})
    last_tool_inputs = state.get("last_tool_inputs", {}).copy()
    last_tool_outputs = state.get("last_tool_outputs", {}).copy()

    if not messages:
        return state


    ai_msgs = [m for m in messages if isinstance(m, AIMessage) and m.tool_calls]
    if not ai_msgs:
        return state
    last_ai = ai_msgs[-1]

    tool_result_messages = []
    for tc in last_ai.tool_calls:
        name = tc["name"]
        args = dict(tc.get("args", {}))


        for k, v in list(args.items()):
            if isinstance(v, str):
                m = re.fullmatch(r"\$(\w+)", v.strip())
                if m:
                    var_name = m.group(1)
                    if var_name in memory:
                        args[k] = memory[var_name]


        try:
            tool = TOOLS.get(name)
            if tool is None:
                result = f"[ToolError] Unknown tool: {name}"
            else:
                result = tool.invoke(args)
        except Exception as e:
            result = f"[ToolError] {type(e).__name__}: {e}"


        last_tool_inputs[name] = args
        last_tool_outputs[name] = result


        tm = ToolMessage(
            content=str(result),
            name=name,
            tool_call_id=tc["id"]
        )
        tool_result_messages.append(tm)


        mem_key = TOOL_OUTPUT_TO_MEMORY_KEY.get(name)
        if mem_key:
            memory[mem_key] = result


    new_messages = messages + tool_result_messages

    return {
        **state,
        "messages": new_messages,
        "memory": memory,
        "last_tool_inputs": last_tool_inputs,
        "last_tool_outputs": last_tool_outputs,
    }


def llm_next_step(state: AgentState) -> str:
    messages = state.get("messages", [])
    if not messages:
        return "end"
    last = messages[-1]
    if isinstance(last, AIMessage) and getattr(last, "tool_calls", None):
        return "tools"
    return "end"


builder = StateGraph(AgentState)

builder.add_node("pipeline", pipeline_node)
builder.add_node("react_llm", llm_node)
builder.add_node("tools", tool_exec_node)

builder.add_conditional_edges(
    START,
    route_start,
    {
        "pipeline": "pipeline",
        "react_llm": "react_llm",
    }
)


builder.add_conditional_edges(
    "react_llm",
    llm_next_step,
    {
        "tools": "tools",
        "end": END
    }
)

builder.add_edge("tools", "react_llm")


checkpointer = MemorySaver()
apps = builder.compile(checkpointer=checkpointer)
