#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 08:56:24 2025

@author: user
"""

import pandas as pd
from scipy import spatial
from typing import List
from openai import OpenAI
import os
import re
from pathlib import Path

openai_api_key = os.environ.get('OPENAI_API_KEY')

client = OpenAI(api_key = openai_api_key)

def get_embedding(text):
    embeddings = []
    response = client.embeddings.create(
        model="text-embedding-3-large",
        input=text
    )
    embeddings.extend([item.embedding for item in response.data])
    return embeddings

def distances_from_embeddings(
        query_embedding: List[float],
        embeddings: List[List[float]],
        distance_metric="cosine",
    ) -> List[List]:
    """Return the distances between a query embedding and a list of embeddings."""
    distance_metrics = {
        "cosine": spatial.distance.cosine,
        "L1": spatial.distance.cityblock,
        "L2": spatial.distance.euclidean,
        "Linf": spatial.distance.chebyshev,
    }
    
    sims = [
        1-distance_metrics[distance_metric](query_embedding, embedding)
        for embedding in embeddings
    ]
    
    return sims

def get_rows_sorted_by_relevance(question, df):
    """
    Function that takes in a question string and a dataframe containing
    rows of text and associated embeddings, and returns that dataframe
    sorted from least to most relevant for that question
    """
    
    # Get embeddings for the question text
    question_embeddings = get_embedding(question)
    list1 = question_embeddings[0]
    # add a "Relevance" column containing
    # the cosine distances between each row's embeddings and the
    # embeddings of the question
    df_copy = df.copy()
    df_copy["Relevance"] = distances_from_embeddings(
        list1,
        df_copy["embedding"].values,
        distance_metric="cosine"
    )
    
    # Sort the copied dataframe by the distances and return it
    # (shorter distance = more relevant so we sort in ascending order)
    df_copy.sort_values("Relevance", ascending=True, inplace=True)
    return df_copy

def evid_eval(file, question):
    df1 = pd.read_csv(file)
    
    try:
        
        df1.loc[:, 'em'] = df1['Title'] + " " + df1['Abstract']
    except:
        print('Concat Error') 
    
    df1 = df1[(df1['Evidence_Hierarchy'] == 'Guideline') | (df1['Evidence_Hierarchy'] == 'SR&MA') | (df1['Evidence_Hierarchy'] == 'RCT')]
    df1.reset_index(drop=True, inplace=True)
    
    text1 = df1['em'].tolist()
    
    embeddings = get_embedding(text1)
    df1['embedding'] = embeddings
    
    df_copy0 = get_rows_sorted_by_relevance(question, df1)
    
    df_copy0 = df_copy0.drop(columns=['em', 'embedding'], errors='ignore')
    
    df_copy0 = df_copy0.sort_values(by='Relevance', ascending=False)
    
    df_copy0.to_csv(file, index=False)
    
    df_copy1 = df_copy0[df_copy0['Evidence_Hierarchy'] == 'Guideline'].copy()
    df_copy2 = df_copy0[df_copy0['Evidence_Hierarchy'] == 'SRMA'].copy()
    df_copy3 = df_copy0[df_copy0['Evidence_Hierarchy'] == 'RCT'].copy()
    
    non_empty_dfs = []
    Evidence_Hierarchys = []
    
    # 检查并收集非空的 DataFrame
    if not df_copy1.empty:
        non_empty_dfs.append(df_copy1)
        Evidence_Hierarchys.append('Guideline')
    if not df_copy2.empty:
        non_empty_dfs.append(df_copy2)
        Evidence_Hierarchys.append('SRMA')
    if not df_copy3.empty:
        non_empty_dfs.append(df_copy3)
        Evidence_Hierarchys.append('RCT')
    
    # 假设 df_copy1, df_copy2, df_copy3 已定义
    # dataframes = [df_copy1, df_copy2, df_copy3]
    Evidence_Hierarchys = ['Guideline', 'SRMA', 'RCT']
    selected_rows = []
    
    # 遍历每个 DataFrame 并处理
    for i, df in enumerate(non_empty_dfs):
    
        # Step 1: 按 distances 降序排序
        sorted_table = df.sort_values(by='Relevance', ascending=False).reset_index(drop=True)
        
        # Step 2: 找出最大年份
        max_year = sorted_table['Year'].max()
        
        # Step 3: 找出 distances 最大的行
        max_distances_row = sorted_table.iloc[0]
        max_distances = max_distances_row['Relevance']
        max_distances_year = max_distances_row['Year']
        
        # Step 4: 如果 distances 最大的行的 Year 也是最大年份，直接选择它
        if max_distances_year == max_year:
            selected_row = max_distances_row
        else:
            # Step 5: 找出所有 Year > max_distances_year 的行
            filtered_rows = sorted_table[sorted_table['Year'] > max_distances_year]
        
            # Step 6: 筛选出这些行中 distances 与 max_distances 差值 < 0.02 的行
            candidate_rows = filtered_rows[filtered_rows['Relevance'] >= max_distances - 0.02]
        
            if not candidate_rows.empty:
                # Step 7: 在候选行中选择 Year 最大的行
                selected_row = candidate_rows.sort_values(by='Year', ascending=False).iloc[0]
            else:
                # Step 8: 如果没有候选行，返回 distances 最大的行
                selected_row = max_distances_row
                
        selected_rows.append(selected_row)
    
    
    # 将选中的行合并为一个 DataFrame
    result_df = pd.DataFrame(selected_rows)
    
    # 显示结果
    print(result_df)
    
    if result_df.empty:
        print("No evidence found in any category.")
    else:
        # 提取 distance 和 PMID（根据存在的类型动态处理）
        # g_distance = result_df[result_df['Evidence_Hierarchy'] == 'Guideline']['Relevance'].values[0] if 'Guideline' in Evidence_Hierarchys else None
        # g_pmid = result_df[result_df['Evidence_Hierarchy'] == 'Guideline']['PMID'].values[0] if 'Guideline' in Evidence_Hierarchys else None
        
        g_df = result_df[result_df['Evidence_Hierarchy'] == 'Guideline']
        g_distance = g_df['Relevance'].values[0] if not g_df.empty else None
        g_pmid = g_df['PMID'].values[0] if not g_df.empty else None
        
        
        srma_df = result_df[result_df['Evidence_Hierarchy'] == 'SRMA']
        s_distance = srma_df['Relevance'].values[0] if not srma_df.empty else None
        s_pmid = srma_df['PMID'].values[0] if not srma_df.empty else None
        # s_pmid = result_df[result_df['Evidence_Hierarchy'] == 'SRMA']['PMID'].values[0] if 'SRMA' in Evidence_Hierarchys else None
        # s_distance = result_df[result_df['Evidence_Hierarchy'] == 'SRMA']['Relevance'].values[0] if 'SRMA' in Evidence_Hierarchys else None
        # s_pmid = result_df[result_df['Evidence_Hierarchy'] == 'SRMA']['PMID'].values[0] if 'SRMA' in Evidence_Hierarchys else None
    
        r_df = result_df[result_df['Evidence_Hierarchy'] == 'RCT']
        r_distance = r_df['Relevance'].values[0] if not r_df.empty else None
        r_pmid = r_df['PMID'].values[0] if not r_df.empty else None
    
        # r_distance = result_df[result_df['Evidence_Hierarchy'] == 'RCT']['Relevance'].values[0] if 'RCT' in Evidence_Hierarchys else None
        # r_pmid = result_df[result_df['Evidence_Hierarchy'] == 'RCT']['PMID'].values[0] if 'RCT' in Evidence_Hierarchys else None
    
        # Step 1: Guideline has priority
        if g_distance is not None and g_distance >= 0.65:
            best_pmid = g_pmid
    
        # Step 2: Compare SRMA and RCT
        else:
            high_sr = s_distance is not None and s_distance >= 0.70
            high_rct = r_distance is not None and r_distance >= 0.70
    
            # Step 2.1: Only one meets threshold
            if high_sr and not high_rct:
                best_pmid = s_pmid
            elif high_rct and not high_sr:
                best_pmid = r_pmid
    
            # Step 2.2: Both meet threshold
            elif high_sr and high_rct:
                diff = r_distance - s_distance
                if diff > 0.03:
                    best_pmid = r_pmid
                else:
                    best_pmid = s_pmid
    
            # Step 3: No evidence meets threshold
            else:
                # 构建存在的 distance 列表
                existing_distances = []
                if g_distance is not None:
                    existing_distances.append(g_distance)
                if s_distance is not None:
                    existing_distances.append(s_distance)
                if r_distance is not None:
                    existing_distances.append(r_distance)
    
                max_evidence = max(existing_distances)
    
                # Compare differences with max_evidence
                if g_distance is not None and (max_evidence - g_distance) < 0.05:
                    best_pmid = g_pmid
                elif s_distance is not None and (max_evidence - s_distance) < 0.03:
                    best_pmid = s_pmid
                else:
                    best_pmid = r_pmid
    
        # 输出最佳证据的 PMID
        print("Best Evidence PMID:", best_pmid)
        
    return best_pmid
# 提取三行的 distances 和 PMID
# g_distances = result_df.iloc[0]['Relevance']
# g_pmid = result_df.iloc[0]['PMID']

# s_distances = result_df.iloc[1]['Relevance']
# s_pmid = result_df.iloc[1]['PMID']

# r_distances = result_df.iloc[2]['Relevance']
# r_pmid = result_df.iloc[2]['PMID']

# # Step 1: Guideline has priority
# if g_distances >= 0.65:
#     best_pmid = g_pmid

# # Step 2: Compare SRMA and RCT
# else:
#     high_sr = s_distances >= 0.70
#     high_rct = r_distances >= 0.70

#     # Step 2.1: Only one meets threshold
#     if high_sr and not high_rct:
#         best_pmid = s_pmid
#     elif high_rct and not high_sr:
#         best_pmid = r_pmid

#     # Step 2.2: Both meet threshold
#     elif high_sr and high_rct:
#         diff = r_distances - s_distances
#         if diff > 0.03:
#             best_pmid = r_pmid
#         else:
#             best_pmid = s_pmid

#     # Step 3: No evidence meets threshold
#     else:
#         max_evidence = max(g_distances, s_distances, r_distances)

#         # Compare differences with max_evidence
#         if (max_evidence - g_distances) < 0.05:
#             best_pmid = g_pmid
#         elif (max_evidence - s_distances) < 0.03:
#             best_pmid = s_pmid
#         else:
#             best_pmid = r_pmid

# # 输出最佳证据的 PMID
# print("Best Evidence PMID:", best_pmid)