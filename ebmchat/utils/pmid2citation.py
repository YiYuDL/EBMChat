#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 16:35:06 2025

@author: user
"""

from Bio import Entrez, Medline
import re
import pandas as pd
from pathlib import Path


def assign_pub_type(list1):
    # 将标题转换为小写，便于匹配
    result = str(list1)[1:-1].replace('"', "'")
    title = result.lower()
    
    pattern = re.compile(r'practice guideline(s)?', re.IGNORECASE)
    
    pattern = re.compile(
        
        r'(practice\s+guidelines?|asco\s+guideline|nccn\s+guidelines|esmo\s+guideline)', 
        flags=re.IGNORECASE
    )
    
    if pattern.search(title):
    
    # # 条件判断
    # if pattern in title or "asco guideline" in title or "nccn guidelines" in title or "esmo guideline" in title:
        return 'Guideline'
    elif 'systematic review' in title or 'meta-analysis' in title:
        return 'SRMA'
    elif any(keyword in title for keyword in ['randomized trial', 'randomized controlled trial', 'phase ii trial', 'phase iii trial', 
                                               'randomized clinical trial', 'rapido trial', 'fowarc trial', 'pcar trial', 
                                               'shortrip trial', 'randomized phase ii', 'randomized phase iii', 'phase ib/ii trial']) \
         and not any(exclude in title for exclude in ['trials', 'protocol', 'pooled analysis']):
        return 'RCT'
    elif any(keyword in title for keyword in ['cohort study', 'cohort analysis', 'retrospective cohort', 'prospective cohort']):
        return 'Cohort'
    elif any(keyword in title for keyword in ['case report', 'case series', 'clinical case', 'a case']):
        return 'Case'
    else:
        return 'Others'

# pmid = 37856213
def pmid2citation(pmid):
    
    pmid = str(pmid)
    Entrez.email = 'qmming@gmail.com'
    file = 'citation.csv'
    SOURCE_DIR = Path("/home/user/source").mkdir(parents=True, exist_ok=True) or Path("/home/user/source")
    pdf_path = SOURCE_DIR / file
    df1 = pd.read_csv(pdf_path)

    mask = df1["PMID"].astype(str) == pmid

    if mask.any():
        citations = df1.loc[mask, "Citation"].dropna()
        citation = citations.values[0]

    else:
        # print('abd')
        pmid = int(pmid)
        handle = Entrez.esearch(db="pubmed", term=pmid,retmax = 1)
        record = Entrez.read(handle)
        handle.close()
        idlist = record['IdList']
        # print(idlist)
        # print(record['Count'])
        
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",retmode="text")
        records = Medline.parse(handle)     
        record = list(records)[0]
        handle.close()
        print("searched Entrez")  
        
        # for record in records:
        title = record.get("TI", "?")
        # author = record.get("AU", "?")
        journal = record.get("TA", "?")
        date = record.get("DP", "?")
        
        start = 0
        while start < len(date) and date[start] == ' ':
            start += 1
        if start < len(date):
            index = date.find(' ', start)
        else:
            index = -1  
        if index != -1:
            year = date[start:index]
        else:
            year = date[start:]
        
        # abstract = record.get("AB", "?")
        pmid = record.get("PMID","?")
        # keywords = record.get("OT", "?")
        PT = record.get("PT", "?")
        original_link = "https://doi.org/" + record.get("LID", "?").split()[0]
        # doi = record.get("LID", "?").split()[0]
        
        PT2 = assign_pub_type(PT)
        PT2 = 'Guideline' if 'asco guideline' in title.lower() else PT2
        
        citation = title + ' ' + str(PT2) + ', '  + journal + ', ' + year + ', ' + original_link + ', ' + str(pmid) + '.'
    
    return citation
# mesh_terms =record.get("MH", "?")
# articles.append({
#     'Title': title,
#     # "Type":PT,
#     'Year': year,
#     'Journal':journal,
#     'Keywords': keywords,
#     'Doi': doi,
#     'Citation':citation,
#     'Pub_type':PT2,
#     # 'Publication Types': pub_types,
#     # 'PubMed Link': pubmed_link,
#     'Original Link': original_link,
#     'Abstract': abstract,
#     # "Type":PT
# }
#     )    