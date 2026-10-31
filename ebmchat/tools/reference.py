"""Targeted search over the bibliography of a previously selected paper."""

from __future__ import annotations

import hashlib
import math
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from langchain_core.messages import HumanMessage, SystemMessage

from utils.model_runtime import (
    EMBEDDING_API_KEY,
    create_chat_model,
    get_embeddings,
)


OPENALEX_BASE_URL = os.getenv(
    "EBMCHAT_OPENALEX_BASE_URL", "https://api.openalex.org"
).rstrip("/")
CONTACT_EMAIL = os.getenv("EBMCHAT_CONTACT_EMAIL", "").strip()
OPENALEX_API_KEY = os.getenv("EBMCHAT_OPENALEX_API_KEY", "").strip()

_TOKEN_RE = re.compile(r"[a-z0-9]+")
_STOPWORDS = {
    "a", "an", "and", "are", "as", "at", "based", "be", "by", "for",
    "from", "in", "is", "it", "literature", "of", "on", "or", "paper",
    "reference", "references", "the", "this", "to", "was", "were", "what",
    "which", "with",
}


def _http_session() -> requests.Session:
    retry = Retry(
        total=3,
        connect=3,
        read=3,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    contact = CONTACT_EMAIL or "contact-not-configured"
    session.headers.update(
        {"User-Agent": f"EBMChat-Ref-Search/1.0 ({contact})"}
    )
    return session


def _openalex_get(
    session: requests.Session,
    path: str,
    params: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    query = dict(params or {})
    if CONTACT_EMAIL:
        query["mailto"] = CONTACT_EMAIL
    if OPENALEX_API_KEY:
        query["api_key"] = OPENALEX_API_KEY
    response = session.get(
        f"{OPENALEX_BASE_URL}/{path.lstrip('/')}",
        params=query,
        timeout=(10, 45),
    )
    if response.status_code == 404:
        raise ValueError(f"OpenAlex could not resolve {path!r}.")
    response.raise_for_status()
    return response.json()


def _chunks(values: List[str], size: int) -> Iterable[List[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def _openalex_id(value: str) -> str:
    return str(value).rstrip("/").rsplit("/", 1)[-1]


def _extract_pmid(ids: Dict[str, Any]) -> Optional[str]:
    value = (ids or {}).get("pmid")
    if not value:
        return None
    match = re.search(r"(\d+)$", str(value).rstrip("/"))
    return match.group(1) if match else None


def _extract_doi(work: Dict[str, Any]) -> Optional[str]:
    value = work.get("doi") or (work.get("ids") or {}).get("doi")
    if not value:
        return None
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", str(value), flags=re.I)


def _abstract_from_inverted_index(index: Any) -> str:
    if not isinstance(index, dict) or not index:
        return ""
    positioned: List[tuple[int, str]] = []
    for token, positions in index.items():
        if not isinstance(positions, list):
            continue
        positioned.extend((int(position), str(token)) for position in positions)
    positioned.sort(key=lambda item: item[0])
    return " ".join(token for _, token in positioned)


def _evidence_hierarchy(title: str, work_type: str) -> str:
    lowered = title.lower()
    if "guideline" in lowered or "practice guideline" in lowered:
        return "Guideline"
    if "systematic review" in lowered or "meta-analysis" in lowered:
        return "SR&MA"
    trial_terms = (
        "randomized", "randomised", "phase 1", "phase 2", "phase 3",
        "phase i", "phase ii", "phase iii", "clinical trial",
    )
    if any(term in lowered for term in trial_terms):
        return "RCT"
    return work_type or "Other"


def _work_to_row(work: Dict[str, Any]) -> Dict[str, Any]:
    ids = work.get("ids") or {}
    title = work.get("display_name") or work.get("title") or ""
    work_type = work.get("type") or ""
    primary_location = work.get("primary_location") or {}
    source = primary_location.get("source") or {}
    return {
        "PMID": _extract_pmid(ids),
        "DOI": _extract_doi(work),
        "OpenAlex_ID": _openalex_id(work.get("id", "")),
        "Title": title,
        "Abstract": _abstract_from_inverted_index(
            work.get("abstract_inverted_index")
        ),
        "Year": work.get("publication_year"),
        "Type": work_type,
        "Evidence_Hierarchy": _evidence_hierarchy(title, work_type),
        "Journal": source.get("display_name") or "",
    }


def fetch_openalex_references(parent_pmid: str) -> List[Dict[str, Any]]:
    """Return OpenAlex metadata for works cited by a PubMed article."""

    parent_pmid = str(parent_pmid).strip()
    if not re.fullmatch(r"\d+", parent_pmid):
        raise ValueError("parent_pmid must contain digits only.")

    session = _http_session()
    parent = _openalex_get(session, f"works/pmid:{parent_pmid}")
    referenced_ids = [
        _openalex_id(value)
        for value in (parent.get("referenced_works") or [])
        if value
    ]
    if not referenced_ids:
        raise ValueError(
            f"No references were available in OpenAlex for PMID {parent_pmid}."
        )

    works: List[Dict[str, Any]] = []
    for batch in _chunks(referenced_ids, 50):
        payload = _openalex_get(
            session,
            "works",
            params={
                "filter": f"openalex_id:{'|'.join(batch)}",
                "per-page": 100,
            },
        )
        works.extend(payload.get("results") or [])

    rows = [_work_to_row(work) for work in works]
    reference_positions = {
        openalex_id: position
        for position, openalex_id in enumerate(referenced_ids)
    }
    rows.sort(
        key=lambda row: reference_positions.get(
            row["OpenAlex_ID"], len(reference_positions)
        )
    )
    return rows


def _tokenize(text: str) -> List[str]:
    return [
        token
        for token in _TOKEN_RE.findall((text or "").lower())
        if len(token) > 1 and token not in _STOPWORDS
    ]


def _bm25_scores(question: str, documents: List[str]) -> List[float]:
    tokenized = [_tokenize(document) for document in documents]
    query_terms = list(dict.fromkeys(_tokenize(question)))
    if not query_terms or not tokenized:
        return [0.0] * len(documents)

    document_count = len(tokenized)
    average_length = sum(len(tokens) for tokens in tokenized) / document_count
    document_frequency = {
        term: sum(term in set(tokens) for tokens in tokenized)
        for term in query_terms
    }
    k1 = 1.5
    b = 0.75
    scores: List[float] = []
    for tokens in tokenized:
        counts = Counter(tokens)
        length_normalizer = 1 - b + b * len(tokens) / max(average_length, 1)
        score = 0.0
        for term in query_terms:
            frequency = counts.get(term, 0)
            if not frequency:
                continue
            frequency_docs = document_frequency[term]
            inverse_frequency = math.log(
                1 + (document_count - frequency_docs + 0.5) / (frequency_docs + 0.5)
            )
            score += inverse_frequency * (
                frequency * (k1 + 1)
            ) / (frequency + k1 * length_normalizer)
        scores.append(score)

    maximum = max(scores, default=0.0)
    if maximum <= 0:
        return [0.0] * len(scores)
    return [score / maximum for score in scores]


def _cosine_scores(question: str, documents: List[str]) -> List[float]:
    embedding_model = get_embeddings()
    document_vectors = np.asarray(
        embedding_model.embed_documents(documents), dtype=float
    )
    question_vector = np.asarray(
        embedding_model.embed_query(question), dtype=float
    )
    document_norms = np.linalg.norm(document_vectors, axis=1)
    question_norm = np.linalg.norm(question_vector)
    denominator = document_norms * question_norm
    similarities = np.divide(
        document_vectors @ question_vector,
        denominator,
        out=np.zeros_like(document_norms),
        where=denominator != 0,
    )
    return [float(value) for value in similarities]


def _rerank_with_llm(
    question: str,
    candidates: List[Dict[str, Any]],
    pool_size: int = 15,
) -> Optional[str]:
    """Choose the best PMID from a small retrieved candidate pool."""

    pool = candidates[: max(1, min(pool_size, len(candidates)))]
    candidate_text = "\n".join(
        f"{position}. PMID {row['PMID']}: {row['Title']}"
        for position, row in enumerate(pool, start=1)
    )
    allowed_pmids = {str(row["PMID"]) for row in pool}
    try:
        model = create_chat_model(temperature=0).bind(max_tokens=2048)
        response = model.invoke(
            [
                SystemMessage(
                    content=(
                        "Select the single cited paper that most directly "
                        "answers the question. Return only its PMID digits, "
                        "with no explanation."
                    )
                ),
                HumanMessage(
                    content=(
                        f"Question: {question}\n\n"
                        f"Candidate cited papers:\n{candidate_text}"
                    )
                ),
            ]
        )
        selected_text = str(response.content).strip()
        match = re.fullmatch(r"\s*(\d+)\s*", selected_text)
        if match and match.group(1) in allowed_pmids:
            return match.group(1)
    except Exception as exc:
        print(f"Reference reranking was unavailable: {exc}")
    return None


def ref_search(parent_pmid: str, question: str, top_k: int = 5) -> Dict[str, Any]:
    """Find references within a selected paper that best answer a follow-up question."""

    question = str(question).strip()
    if not question:
        raise ValueError("question must not be empty.")
    top_k = max(1, min(int(top_k), 20))

    rows = fetch_openalex_references(parent_pmid)
    candidates = [row for row in rows if row.get("PMID") and row.get("Title")]
    if not candidates:
        raise ValueError(
            f"No PubMed-linked references were available for PMID {parent_pmid}."
        )

    documents = [
        f"{row['Title']} {row.get('Abstract') or ''}".strip()
        for row in candidates
    ]
    lexical_scores = _bm25_scores(question, documents)

    if EMBEDDING_API_KEY:
        semantic_scores = _cosine_scores(question, documents)
        combined_scores = [
            0.9 * semantic + 0.1 * lexical
            for semantic, lexical in zip(semantic_scores, lexical_scores)
        ]
        ranking_method = "openai_embedding_cosine_plus_bm25"
    else:
        semantic_scores = [None] * len(candidates)
        combined_scores = lexical_scores
        ranking_method = "bm25_fallback_no_embedding_key"

    for row, semantic, lexical, combined in zip(
        candidates, semantic_scores, lexical_scores, combined_scores
    ):
        row["Semantic_Relevance"] = semantic
        row["Lexical_Relevance"] = float(lexical)
        row["Relevance"] = float(combined)

    candidates.sort(key=lambda row: row["Relevance"], reverse=True)
    selected_pmid = _rerank_with_llm(question, candidates)
    if selected_pmid:
        selected_index = next(
            index
            for index, row in enumerate(candidates)
            if str(row["PMID"]) == selected_pmid
        )
        selected = candidates.pop(selected_index)
        candidates.insert(0, selected)
        ranking_method += "_llm_rerank"
    else:
        selected = candidates[0]
        ranking_method += "_llm_unavailable"

    for rank, row in enumerate(candidates, start=1):
        row["Final_Rank"] = rank
        row["LLM_Selected"] = str(row["PMID"]) == str(selected["PMID"])

    table_dir = Path.cwd() / "table"
    table_dir.mkdir(parents=True, exist_ok=True)
    query_hash = hashlib.sha256(question.encode("utf-8")).hexdigest()[:10]
    table_path = table_dir / f"references_{parent_pmid}_{query_hash}.csv"
    pd.DataFrame(candidates).to_csv(table_path, index=False)

    visible_fields = (
        "PMID", "DOI", "Title", "Year", "Evidence_Hierarchy", "Journal",
        "Semantic_Relevance", "Lexical_Relevance", "Relevance",
        "Final_Rank", "LLM_Selected",
    )
    top_candidates = [
        {field: row.get(field) for field in visible_fields}
        for row in candidates[:top_k]
    ]
    return {
        "parent_pmid": str(parent_pmid),
        "question": question,
        "reference_count": len(rows),
        "pubmed_reference_count": len(candidates),
        "selected_pmid": selected["PMID"],
        "selected_doi": selected.get("DOI"),
        "selected_title": selected["Title"],
        "relevance": selected["Relevance"],
        "ranking_method": ranking_method,
        "references_file": str(table_path),
        "top_candidates": top_candidates,
        "source": "openalex_referenced_works",
    }
