"""
Created on Fri Sep  5 17:00:27 2025

@author: user
"""

import pandas as pd


def result_show(pub_type: str, file_path: str) -> pd.DataFrame:
    """
    Filter records in a CSV file by evidence type and shorten long titles for display.

    Args:
        pub_type (str): Evidence type to filter by.
        file_path (str): Path to the CSV file.

    Returns:
        pd.DataFrame: Filtered results containing PMID, Title, Year,
        Evidence_Hierarchy, and Relevance columns.
    """

    df = pd.read_csv(file_path)


    df1 = df[df['Evidence_Hierarchy'] == pub_type].copy()


    if df1.empty:
        print(f"未找到 Pub_type 为 '{pub_type}' 的记录")
        return df1


    required_columns = ['PMID', 'Title', 'Year', 'Evidence_Hierarchy', 'Relevance']
    for col in required_columns:
        if col not in df1.columns:
            df1[col] = pd.NA


    df1 = df1[required_columns]


    if 'Title' in df1.columns:
        def truncate_title(title, max_length=70):
            """Truncate titles to max_length characters and append an ellipsis."""
            if pd.isna(title) or not isinstance(title, str):
                return title

            if len(title) <= max_length:
                return title

            return title[:max_length] + "..."


        df1['Title'] = df1['Title'].apply(truncate_title)


    print(df1.to_string(
        index=False,
        justify='left',
        max_colwidth=100,
        line_width=500
    ))

    ab = df1.to_string(
        index=False,
        justify='left',
        max_colwidth=100,
        line_width=500)

    return ab
