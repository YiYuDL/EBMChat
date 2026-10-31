"""OpenAI model configuration shared by the EBMChat agent and RAG pipeline."""

import os
from functools import lru_cache

from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from openai import OpenAI


OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
LLM_MODEL = os.getenv("EBMCHAT_LLM_MODEL", "gpt-4.1")
LLM_MAX_TOKENS = int(os.getenv("EBMCHAT_LLM_MAX_TOKENS", "32768"))
LLM_TEMPERATURE = float(os.getenv("EBMCHAT_LLM_TEMPERATURE", "0.2"))

EMBEDDING_MODEL = os.getenv("EBMCHAT_EMBEDDING_MODEL", "text-embedding-3-large")
EMBEDDING_BASE_URL = os.getenv(
    "EBMCHAT_EMBEDDING_BASE_URL", "https://api.openai.com/v1"
)
EMBEDDING_API_KEY = os.getenv("EBMCHAT_EMBEDDING_API_KEY") or OPENAI_API_KEY


def create_chat_model(temperature: float | None = None) -> ChatOpenAI:
    """Create the OpenAI chat model used throughout EBMChat."""

    return ChatOpenAI(
        model=LLM_MODEL,
        api_key=OPENAI_API_KEY,
        temperature=LLM_TEMPERATURE if temperature is None else temperature,
        max_tokens=LLM_MAX_TOKENS,
        max_retries=2,
    )


def create_openai_client() -> OpenAI:
    """Create the OpenAI client used by query generation."""

    return OpenAI(api_key=OPENAI_API_KEY)


@lru_cache(maxsize=1)
def get_embeddings() -> OpenAIEmbeddings:
    """Return one shared OpenAI embedding client for indexing and querying."""

    return OpenAIEmbeddings(
        model=EMBEDDING_MODEL,
        base_url=EMBEDDING_BASE_URL,
        api_key=EMBEDDING_API_KEY,
    )
