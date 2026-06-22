"""
Created on Thu Aug 14 17:03:54 2025

@author: user
"""

import sys
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent


if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from pathlib import Path
from typing import Union
from utils.download import get_pdf_link_by_pmid,download_by_link
from utils.pmid2citation import pmid2citation
from utils.answer_rag import RAG


def answer_gen(pmid: str, question: str) -> str:
    """
    Generate an answer for a question using the paper identified by PMID.

    1. If source/PMID.pdf exists, call RAG(PMID.pdf, question) directly.
    2. If it does not exist, try to download the paper.
    3. If the download succeeds, call RAG(PMID.pdf, question) again.
    4. If the download still fails, return a message asking for manual download.
    """

    citation = pmid2citation(pmid)

    pmid_str = str(pmid).strip()
    pdf_name = f"{pmid_str}.pdf"
    source_dir = Path(__file__).resolve().parent.parent / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = source_dir / pdf_name


    if pdf_path.exists():
        output = RAG(pmid, question)

        print('\n')
        output = output + '\n\n' + 'Citation: ' + citation
        print(output)
        return output

    else:

        try:
            link = get_pdf_link_by_pmid(pmid)
            dl_result = download_by_link(link,pmid)

        except Exception:
            dl_result = None


        success = False
        final_pdf_path = pdf_path

        if isinstance(dl_result, bool):
            success = dl_result
        elif isinstance(dl_result, (str, Path)):
            p = Path(dl_result)

            if p.is_dir():
                candidate = p / pdf_name
            else:
                candidate = p
            if candidate.suffix.lower() == ".pdf" and candidate.exists():
                final_pdf_path = candidate
                success = True


        if not success and pdf_path.exists():
            final_pdf_path = pdf_path
            success = True


        if success:

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
