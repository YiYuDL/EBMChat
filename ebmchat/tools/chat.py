#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:00:27 2025

@author: user
"""

import pandas as pd
# import numpy as np

def result_show(pub_type: str, file_path: str) -> pd.DataFrame:
    """
    根据文献类型筛选CSV文件中的记录，并智能处理长标题
    
    参数:
    pub_type (str): 要筛选的文献类型
    file_path (str): CSV文件路径
    
    返回:
    pd.DataFrame: 包含PMID, Title, Year, Pub_type, distances列的筛选结果
    """
    # 读取CSV文件
    df = pd.read_csv(file_path)
    
    # 筛选Pub_type匹配的行
    df1 = df[df['Evidence_Hierarchy'] == pub_type].copy()
    
    # 检查是否找到匹配记录
    if df1.empty:
        print(f"未找到 Pub_type 为 '{pub_type}' 的记录")
        return df1
    
    # 确保包含所需列
    required_columns = ['PMID', 'Title', 'Year', 'Evidence_Hierarchy', 'Relevance']
    for col in required_columns:
        if col not in df1.columns:
            df1[col] = pd.NA
    
    # 仅保留需要的列
    df1 = df1[required_columns]
    
    # 智能截断Title列（保留首尾，中间用...替代）
    if 'Title' in df1.columns:
        def truncate_title(title, max_length=70):
            """截断标题：只保留开头max_length个字符，超出部分用...替代"""
            if pd.isna(title) or not isinstance(title, str):
                return title
            
            if len(title) <= max_length:
                return title
            
            return title[:max_length] + "..."
        
        # 应用智能截断
        df1['Title'] = df1['Title'].apply(truncate_title)
    
    # 打印所有行（不显示索引，确保单行显示）
    print(df1.to_string(
        index=False,
        justify='left',
        max_colwidth=100,  # 防止其他列过长
        line_width=500      # 确保整行在单行内显示
    ))
    
    ab = df1.to_string(
        index=False,
        justify='left',
        max_colwidth=100,  # 防止其他列过长
        line_width=500)
    
    return ab