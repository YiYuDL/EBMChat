#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:59:38 2024

@author: user
"""

from Bio import Entrez, Medline
import pandas as pd
from pathlib import Path
import os 
import re

'''
final_query = '"therapeutics"[Mesh] AND "rectal neoplasms"[Mesh] AND ("Intermediate"[TIAB] OR"Intermediate-stage"[TIAB] OR "Stage II"[TIAB] OR "Stage III"[TIAB] OR "Locally advanced"[TIAB] OR "Regional spread"[TIAB] OR "T3"[TIAB] OR "T4"[TIAB] OR "Invasive"[TIAB] OR "Unresectable"[TIAB]) AND 2021:3000[dp]'
final_query2 = '"therapeutics"[Mesh] AND "rectal neoplasms"[Mesh] AND ("advanced"[TIAB] OR "metastatic"[TIAB] OR "stage IV"[TIAB] OR "advanced stage"[TIAB] OR "unresectable"[TIAB] OR "terminal"[TIAB] OR "end-stage"[TIAB]) NOT "locally advanced" AND 2021:3000[dp]'

final_query3 = '"therapeutics"[Mesh] AND "rectal neoplasms"[Mesh] AND ("advanced"[TIAB] OR "metastatic"[TIAB] OR "stage IV"[TIAB] OR "advanced stage"[TIAB] OR "unresectable"[TIAB] OR "terminal"[TIAB] OR "end-stage"[TIAB]) AND 2023:3000[dp]'
final_query4 = '"therapeutics"[Mesh] AND "rectal neoplasms"[Mesh] AND ("Intermediate"[TIAB] OR"Intermediate-stage"[TIAB] OR "Stage II"[TIAB] OR "Stage III"[TIAB] OR "Locally advanced"[TIAB] OR "Regional spread"[TIAB] OR "T3"[TIAB] OR "T4"[TIAB] OR "Invasive"[TIAB] OR "Unresectable"[TIAB]) AND 2016:3000[dp]'
'''
# final_query4 = '"therapeutics"[Mesh] AND "rectal neoplasms"[Mesh] AND ("Intermediate"[TIAB] OR"Intermediate-stage"[TIAB] OR "Stage II"[TIAB] OR "Stage III"[TIAB] OR "Locally advanced"[TIAB] OR "Regional spread"[TIAB] OR "T3"[TIAB] OR "T4"[TIAB] OR "Invasive"[TIAB] OR "Unresectable"[TIAB]) AND 2016:3000[dp]'



#SOURCE_DIR = Path("/home/user/table").mkdir(parents=True, exist_ok=True) or Path("/home/user/table")

def assign_pub_type(list1):
    # 将标题转换为小写，便于匹配
    result = str(list1)[1:-1].replace('"', "'")
    title = result.lower()

    # 条件判断
    if 'practice guideline' in title:
        return 'Guideline'
    elif 'systematic review' in title or 'meta-analysis' in title:
        return 'SR&MA'
    elif any(keyword in title for keyword in ['randomized trial', 'randomized controlled trial', 'phase ii trial', 'phase iii trial', 
                                               'randomized clinical trial', 'rapido trial', 'fowarc trial', 'pcar trial', 
                                               'shortrip trial', 'randomized phase ii', 'randomized phase iii', 'phase ib/ii trial']) \
         and not any(exclude in title for exclude in ['trials', 'protocol', 'pooled analysis']):
        return 'RCT'
    elif any(keyword in title for keyword in ['cohort study', 'cohort analysis', 'retrospective cohort', 'prospective cohort']):
        return 'cohort'
    elif any(keyword in title for keyword in ['case report', 'case series', 'clinical case', 'a case']):
        return 'case'
    else:
        return 'others'

def pm_search(query):

    '''input your own email address'''
    Entrez.email = 'yourownemial@gmail.com'
    
    handle = Entrez.esearch(db="pubmed", term=query,retmax = 10000)
    record = Entrez.read(handle)
    handle.close()
    idlist = record['IdList']
    # print(idlist)
    # print(record['Count'])
    
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline",retmode="text")
    records = Medline.parse(handle)     
    records = list(records)
    handle.close()
    print("searched Entrez")  
    articles = []
    
    for record in records:
        title = record.get("TI", "?")
        # author = record.get("AU", "?")
        pmid = record.get("PMID","?")
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
        
        abstract = record.get("AB", "?")
        keywords = record.get("OT", "?")
        PT = record.get("PT", "?")
        original_link = "https://doi.org/" + record.get("LID", "?").split()[0]
        doi = record.get("LID", "?").split()[0]
        
        PT2 = assign_pub_type(PT)
        PT2 = 'Guideline' if 'asco guideline' in title.lower() else PT2
        
        citation = title + ' ' + str(PT2) + ', '  + journal + ', ' + year + ', ' + original_link +', ' + pmid
        # mesh_terms =record.get("MH", "?")
        articles.append({
            'Title': title,
            # "Type":PT,
            'Year': year,
            'Journal':journal,
            'Keywords': keywords,
            'Doi': doi,
            'Citation':citation,
            'Evidence_Hierarchy':PT2,
            # 'Publication Types': pub_types,
            # 'PubMed Link': pubmed_link,
            'Original Link': original_link,
            'PMID': pmid,
            'Abstract': abstract,
            # "Type":PT
        }
            )    
    
    
    df1 = pd.DataFrame(articles)
    
    df1['Doi'] = df1['Doi'].replace('?', pd.NA)      
    df1['Doi'] = df1['Doi'].replace(r'^\s*$', pd.NA, regex=True)  

    df1 = df1.dropna(subset=['Doi'])
    
    #parent_dir = Path.cwd().parent
    parent_dir = Path.cwd()
    source_dir = parent_dir / "table"

    if not source_dir.exists():
        source_dir.mkdir()
        print(f"table folder: {source_dir} is created")
    else:
        print(f"table folder: {source_dir} is existed")
    
    directory = str(source_dir)
    max_num = -1
    for filename in os.listdir(directory):
        match = re.match(r"table_(\d+)\.csv", filename)
        if match:
            num = int(match.group(1))
            if num > max_num:
                max_num = num
    if max_num == -1:
        new_file = "table_0.csv"
    else:
        new_file = f"table_{max_num + 1}.csv"
    
    file = directory + '/' + new_file
    # file = SOURCE_DIR / new_file
    df1.to_csv(file, index=False)
    searched_result = str(file)
    
    print(f'Searched result is saved in {searched_result}')
    return searched_result

# parent_dir = Path.cwd().parent
# source_dir = parent_dir / "table"

# if not source_dir.exists():
#     source_dir.mkdir()
#     print(f"table folder: {source_dir} is created")
# else:
#     print(f"table folder: {source_dir} is existed")


'''
out = Query2Table(final_query4,46)
'''
# out = Query2Table(final_query4,148)
