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


    question_embeddings = get_embedding(question)
    list1 = question_embeddings[0]


    df_copy = df.copy()
    df_copy["Relevance"] = distances_from_embeddings(
        list1,
        df_copy["embedding"].values,
        distance_metric="cosine"
    )


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


    if not df_copy1.empty:
        non_empty_dfs.append(df_copy1)
        Evidence_Hierarchys.append('Guideline')
    if not df_copy2.empty:
        non_empty_dfs.append(df_copy2)
        Evidence_Hierarchys.append('SRMA')
    if not df_copy3.empty:
        non_empty_dfs.append(df_copy3)
        Evidence_Hierarchys.append('RCT')


    Evidence_Hierarchys = ['Guideline', 'SRMA', 'RCT']
    selected_rows = []


    for i, df in enumerate(non_empty_dfs):


        sorted_table = df.sort_values(by='Relevance', ascending=False).reset_index(drop=True)


        max_year = sorted_table['Year'].max()


        max_distances_row = sorted_table.iloc[0]
        max_distances = max_distances_row['Relevance']
        max_distances_year = max_distances_row['Year']


        if max_distances_year == max_year:
            selected_row = max_distances_row
        else:

            filtered_rows = sorted_table[sorted_table['Year'] > max_distances_year]


            candidate_rows = filtered_rows[filtered_rows['Relevance'] >= max_distances - 0.02]

            if not candidate_rows.empty:

                selected_row = candidate_rows.sort_values(by='Year', ascending=False).iloc[0]
            else:

                selected_row = max_distances_row

        selected_rows.append(selected_row)


    result_df = pd.DataFrame(selected_rows)


    print(result_df)

    if result_df.empty:
        print("No evidence found in any category.")
    else:


        g_df = result_df[result_df['Evidence_Hierarchy'] == 'Guideline']
        g_distance = g_df['Relevance'].values[0] if not g_df.empty else None
        g_pmid = g_df['PMID'].values[0] if not g_df.empty else None


        srma_df = result_df[result_df['Evidence_Hierarchy'] == 'SRMA']
        s_distance = srma_df['Relevance'].values[0] if not srma_df.empty else None
        s_pmid = srma_df['PMID'].values[0] if not srma_df.empty else None


        r_df = result_df[result_df['Evidence_Hierarchy'] == 'RCT']
        r_distance = r_df['Relevance'].values[0] if not r_df.empty else None
        r_pmid = r_df['PMID'].values[0] if not r_df.empty else None


        if g_distance is not None and g_distance >= 0.65:
            best_pmid = g_pmid


        else:
            high_sr = s_distance is not None and s_distance >= 0.70
            high_rct = r_distance is not None and r_distance >= 0.70


            if high_sr and not high_rct:
                best_pmid = s_pmid
            elif high_rct and not high_sr:
                best_pmid = r_pmid


            elif high_sr and high_rct:
                diff = r_distance - s_distance
                if diff > 0.03:
                    best_pmid = r_pmid
                else:
                    best_pmid = s_pmid


            else:

                existing_distances = []
                if g_distance is not None:
                    existing_distances.append(g_distance)
                if s_distance is not None:
                    existing_distances.append(s_distance)
                if r_distance is not None:
                    existing_distances.append(r_distance)

                max_evidence = max(existing_distances)


                if g_distance is not None and (max_evidence - g_distance) < 0.05:
                    best_pmid = g_pmid
                elif s_distance is not None and (max_evidence - s_distance) < 0.03:
                    best_pmid = s_pmid
                else:
                    best_pmid = r_pmid


        print("Best Evidence PMID:", best_pmid)

    return best_pmid
