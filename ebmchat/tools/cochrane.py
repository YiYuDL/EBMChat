"""Search Cochrane Review records and save normalized results as CSV."""

from __future__ import annotations

import html
import io
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd
import requests
from Bio import Medline
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from tools.query import cochrane_query


CROSSREF_BASE_URL = os.getenv(
    "EBMCHAT_CROSSREF_BASE_URL", "https://api.crossref.org"
).rstrip("/")
CONTACT_EMAIL = os.getenv("EBMCHAT_CONTACT_EMAIL", "").strip()
COCHRANE_ISSN = "1465-1858"
COCHRANE_JOURNAL = "Cochrane Database of Systematic Reviews"

_COCHRANE_ID_RE = re.compile(r"\b(CD\d{6,})\b", re.IGNORECASE)
_VERSION_RE = re.compile(r"\.pub(\d+)\b", re.IGNORECASE)
_HTML_TAG_RE = re.compile(r"<[^>]+>")


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
        {"User-Agent": f"EBMChat-Coc-Search/1.0 ({contact})"}
    )
    return session


def _chunks(values: List[str], size: int) -> Iterable[List[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def _first(value: Any, default: str = "") -> str:
    if isinstance(value, list):
        return str(value[0]) if value else default
    return str(value) if value is not None else default


def _plain_search_terms(query: str) -> str:
    """Turn Boolean syntax into terms suitable for Crossref relevance search."""

    terms = re.sub(r"\b(?:AND|OR|NOT)\b", " ", query, flags=re.IGNORECASE)
    terms = re.sub(r"[()\[\]{}\"]", " ", terms)
    return re.sub(r"\s+", " ", terms).strip()


def _cochrane_id(doi: str) -> str:
    match = _COCHRANE_ID_RE.search(doi or "")
    return match.group(1).upper() if match else ""


def _version(doi: str) -> int:
    match = _VERSION_RE.search(doi or "")
    return int(match.group(1)) if match else 0


def _published_year(item: Dict[str, Any]) -> Optional[int]:
    for key in ("published", "published-online", "published-print", "created"):
        parts = (item.get(key) or {}).get("date-parts") or []
        if parts and parts[0]:
            try:
                return int(parts[0][0])
            except (TypeError, ValueError):
                continue
    return None


def _clean_abstract(value: Any) -> str:
    text = _first(value)
    if not text:
        return ""
    text = _HTML_TAG_RE.sub(" ", text)
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def _author_text(authors: Any) -> str:
    names = []
    for author in authors or []:
        if not isinstance(author, dict):
            continue
        name = " ".join(
            part for part in (author.get("given"), author.get("family")) if part
        )
        if name:
            names.append(name)
    return "; ".join(names)


def _deduplicate_versions(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Keep the newest registered version of each Cochrane Review."""

    selected: Dict[str, Dict[str, Any]] = {}
    for item in items:
        doi = str(item.get("DOI") or "")
        review_id = _cochrane_id(doi) or doi.lower()
        candidate = dict(item)
        candidate["_search_score"] = float(item.get("score") or 0.0)
        existing = selected.get(review_id)
        if existing is None:
            selected[review_id] = candidate
            continue

        maximum_score = max(
            float(existing.get("_search_score") or 0.0),
            float(candidate.get("_search_score") or 0.0),
        )
        if _version(doi) > _version(str(existing.get("DOI") or "")):
            candidate["_search_score"] = maximum_score
            selected[review_id] = candidate
        else:
            existing["_search_score"] = maximum_score

    return sorted(
        selected.values(),
        key=lambda item: float(item.get("_search_score") or 0.0),
        reverse=True,
    )


def _search_crossref(
    session: requests.Session, query: str, max_results: int
) -> List[Dict[str, Any]]:
    params: Dict[str, Any] = {
        "query.bibliographic": _plain_search_terms(query),
        # Fetch extra records because protocol and review versions share a CD id.
        "rows": min(max(max_results * 2, 20), 1000),
    }
    if CONTACT_EMAIL:
        params["mailto"] = CONTACT_EMAIL
    response = session.get(
        f"{CROSSREF_BASE_URL}/journals/{COCHRANE_ISSN}/works",
        params=params,
        timeout=(10, 45),
    )
    response.raise_for_status()
    items = (response.json().get("message") or {}).get("items") or []
    return _deduplicate_versions(items)[:max_results]


def _resolve_pmids(
    session: requests.Session, dois: List[str]
) -> Dict[str, str]:
    resolved: Dict[str, str] = {}
    for batch in _chunks([doi for doi in dois if doi], 100):
        params: Dict[str, Any] = {
            "ids": ",".join(batch),
            "format": "json",
            "tool": "EBMChat",
        }
        if CONTACT_EMAIL:
            params["email"] = CONTACT_EMAIL
        response = session.get(
            "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
            params=params,
            timeout=(10, 45),
        )
        response.raise_for_status()
        for record in response.json().get("records") or []:
            doi = str(record.get("doi") or record.get("requested-id") or "")
            pmid = record.get("pmid")
            if doi and pmid:
                resolved[doi.lower()] = str(pmid)
    return resolved


def _fetch_pubmed_metadata(
    session: requests.Session, pmids: List[str]
) -> Dict[str, Dict[str, Any]]:
    metadata: Dict[str, Dict[str, Any]] = {}
    for batch in _chunks([pmid for pmid in pmids if pmid], 100):
        params: Dict[str, Any] = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "medline",
            "retmode": "text",
            "tool": "EBMChat",
        }
        if CONTACT_EMAIL:
            params["email"] = CONTACT_EMAIL
        response = session.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params=params,
            timeout=(10, 60),
        )
        response.raise_for_status()
        for record in Medline.parse(io.StringIO(response.text)):
            pmid = str(record.get("PMID") or "")
            if pmid:
                metadata[pmid] = record
    return metadata


def _evidence_hierarchy(title: str, doi: str) -> str:
    lowered = title.lower()
    if "protocol" in lowered or _version(doi) == 0:
        return "Protocol"
    if "withdrawn" in lowered:
        return "Withdrawn"
    return "SR&MA"


def _items_to_rows(
    items: List[Dict[str, Any]],
    pmid_by_doi: Dict[str, str],
    pubmed_by_pmid: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    maximum_score = max(
        (float(item.get("_search_score") or 0.0) for item in items),
        default=0.0,
    )
    rows: List[Dict[str, Any]] = []
    for item in items:
        doi = str(item.get("DOI") or "")
        pmid = pmid_by_doi.get(doi.lower(), "")
        pubmed = pubmed_by_pmid.get(pmid, {})
        title = _first(item.get("title")) or str(pubmed.get("TI") or "")
        year = _published_year(item)
        if not year:
            year_match = re.search(r"\b(19|20)\d{2}\b", str(pubmed.get("DP") or ""))
            year = int(year_match.group(0)) if year_match else None
        journal = _first(item.get("container-title"), COCHRANE_JOURNAL)
        abstract = str(pubmed.get("AB") or "") or _clean_abstract(
            item.get("abstract")
        )
        keywords = pubmed.get("OT") or item.get("subject") or []
        link = str(item.get("URL") or f"https://doi.org/{doi}")
        hierarchy = _evidence_hierarchy(title, doi)
        score = float(item.get("_search_score") or 0.0)
        relevance = score / maximum_score if maximum_score else 0.0
        citation = (
            f"{title} {hierarchy}, {journal}, {year or '?'}, {link}, "
            f"{pmid or _cochrane_id(doi)}"
        )
        rows.append(
            {
                "Title": title,
                "Year": year,
                "Journal": journal,
                "Authors": _author_text(item.get("author")),
                "Keywords": keywords,
                "Doi": doi,
                "Cochrane_ID": _cochrane_id(doi),
                "Citation": citation,
                "Evidence_Hierarchy": hierarchy,
                "Original Link": link,
                "PMID": pmid or pd.NA,
                "Abstract": abstract,
                "Relevance": relevance,
            }
        )
    return rows


def _next_output_file() -> Path:
    table_dir = Path.cwd() / "table"
    table_dir.mkdir(parents=True, exist_ok=True)
    maximum = -1
    for path in table_dir.glob("cochrane_*.csv"):
        match = re.fullmatch(r"cochrane_(\d+)\.csv", path.name)
        if match:
            maximum = max(maximum, int(match.group(1)))
    return table_dir / f"cochrane_{maximum + 1}.csv"


def coc_search(query: str, max_results: int = 100) -> str:
    """Search Cochrane Reviews, rank them, and save a normalized CSV table.

    Search candidates come from the Crossref records registered for the
    Cochrane Database of Systematic Reviews (ISSN 1465-1858).  PubMed metadata
    is added when a DOI-to-PMID mapping exists.  The returned CSV follows the
    schema used by PM_Search and includes a normalized relevance score.
    """

    normalized_query = cochrane_query(query)
    max_results = max(1, min(int(max_results), 500))
    session = _http_session()
    items = _search_crossref(session, normalized_query, max_results)
    if not items:
        raise ValueError("No Cochrane Review records matched the query.")

    dois = [str(item.get("DOI") or "") for item in items]
    pmid_by_doi: Dict[str, str] = {}
    pubmed_by_pmid: Dict[str, Dict[str, Any]] = {}
    try:
        pmid_by_doi = _resolve_pmids(session, dois)
        pubmed_by_pmid = _fetch_pubmed_metadata(
            session, list(pmid_by_doi.values())
        )
    except (requests.RequestException, ValueError) as exc:
        # Crossref results remain useful if optional PubMed enrichment is down.
        print(f"PubMed enrichment was unavailable: {exc}")

    rows = _items_to_rows(items, pmid_by_doi, pubmed_by_pmid)
    output_file = _next_output_file()
    pd.DataFrame(rows).to_csv(output_file, index=False)
    print(f"Cochrane search result is saved in {output_file}")
    return str(output_file)
