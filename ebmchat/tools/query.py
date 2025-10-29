#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 15:03:29 2025

@author: user
"""

import openai

client = openai.OpenAI()

def generate_query(question: str) -> str:
    prompt = f"""You are an expert in evidence-based medicine, skilled at using evidence-based medicine thinking to convert clinical questions described in natural language into search queries based on the PICO principle (patient, intervention, condition, and outcome). After conversion, please place patient- or disease-related keywords first, each keyword needs to be enclosed in quotation marks.. Below are a few representative examples; please follow these examples for conversion：
1. question: What are recommended treatments for advanced non-small cell lung cancer. query: "non-small cell lung cancer " AND "advanced" AND "treatments";
2.question: What is the diagnosis of early breast cancer. query: "breast cancer" AND "early" AND "diagnosis";
3.question: What is the recommended prognosis of metastatic breast cancer. query: "breast cancer" AND "metastatic" AND "prognosis";
4.question: what is the optimal first-line regimen for  metastatic non-small-cell lung carcinomas patients with epidermal growth factor receptor mutation. query: "non-small-cell lung carcinomas" AND "metastatic" AND "first-line regimen" AND "epidermal growth factor receptor mutation";
5.question: Do combining Encorafenib and cetuximab have better effect, compared with solely drug therapy. query: "colorectal cancer" AND ("Encorafenib plus cetuximab" OR "Encorafenib, cetuximab") AND "drug therapy" AND "BRAF-mutant";
6.question: Tell me about the cost-effectiveness of nivolumab for patients with non-small cell lung cancer. query: "non-small-cell lung cancer" AND "cost-effectiveness" AND "nivolumab";
7. question: what is the proper approach among perioperative patients with T1 non-small-cell lung cancer. query: "non-small-cell lung cancer" AND "perioperative" AND "pembrolizumab" AND "T1";

Please follow the instructions and examples above and only output the final converted search query, without any additional explanations or additional content. Here the question is {question}
"""
    response = client.chat.completions.create(
        model="gpt-4.1",
        messages=[{"role": "user", "content": prompt}],
        temperature=0
    )
    return response.choices[0].message.content.strip()

import re
from Bio import Entrez

def replace(Query):
    words_and_operators = re.findall(r'\"[^\"]*\"|\S+', Query)
    
    list1 = ['"advanced"', '"metastatic"', '"stage IV""', '"advanced stage"', '"unresectable"', '"terminal"', '"end-stage"']
    replace1 = '(("advanced"[TIAB] OR "metastatic"[TIAB] OR "stage IV"[TIAB] OR "advanced stage"[TIAB] OR "unresectable"[TIAB] OR "terminal"[TIAB] OR "end-stage"[TIAB]) NOT "locally advanced"[Title])'
    # list2 = ['"Intermediate"','"Intermediate-stage"', '"Stage II"', '"Stage III"', '"Locally advanced"', '"Regional spread"', '"T3"', '"T4"', '"Invasive"', '"Unresectable"']
    # replace2 = '("Intermediate"[TIAB] OR"Intermediate-stage"[TIAB] OR "Stage II"[TIAB] OR "Stage III"[TIAB] OR "Locally advanced"[TIAB] OR "Regional spread"[TIAB] OR "T3"[TIAB] OR "T4"[TIAB] OR "Invasive"[TIAB] OR "Unresectable"[TIAB])'
    list3 = ['"Intermediate"','"Stage II"', '"Stage III"', '"T3"', '"T4"', '"Locally advanced"', '"Invasive"', '"Unresectable"',
        '"non-metastatic"','"nonmetastatic"','"early"', '"stage I"', '"incipient"', '"carcinoma in situ"', '"T1"', '"T2"', '"pre-cancerous"', '"early-stage"']
    replace3 = '("Intermediate"[TIAB] OR "Stage II"[TIAB] OR "Stage III"[TIAB] OR "Locally advanced"[TIAB] OR "Invasive"[TIAB] OR "Unresectable"[TIAB] OR "T3"[TIAB] OR "T4"[TIAB] OR "early"[TIAB] OR "stage I"[TIAB] OR "incipient"[TIAB] OR "carcinoma in situ"[TIAB] OR "T1"[TIAB] OR "T2"[TIAB] OR "pre-cancerous"[TIAB] OR "early-stage"[TIAB] OR "non-metastatic"[TIAB] OR "nonmetastatic"[TIAB])'
    
    words_and_operators = re.findall(r'\"[^\"]*\"|\S+', Query)
    
    # print(words_and_operators)
    
    # 遍历 string1 中的每个元素（单词或分隔符）
    for i, word in enumerate(words_and_operators):
        # print('i:',i)
        # print('word:',word)
        # 检查这个单词是否在 list1 中，大小写不敏感
        if any(re.fullmatch(re.escape(item),  word, flags=re.IGNORECASE) for item in list1):
            # 如果匹配，替换为 replace1
            words_and_operators[i] = replace1
            
        # elif any(re.fullmatch(re.escape(item),  word, flags=re.IGNORECASE) for item in list2):
        #     words_and_operators[i] = replace2
            
        elif any(re.fullmatch(re.escape(item),  word, flags=re.IGNORECASE) for item in list3):
            words_and_operators[i] = replace3
    # 将处理后的单词和分隔符重新合并为一个字符串
    result_string = " ".join(words_and_operators)
    return result_string

'''
mesh
'''

# 定义function1
def Mesh(input_string):
    try:
        input_string = input_string.strip('"')
        
        '''input your own email address'''
        Entrez.email = 'yourownemial@gmail.com'
       
        handle = Entrez.esearch(db="mesh", term=input_string)
        record = Entrez.read(handle)
        if len(record['TranslationSet']) == 1:
            translation = record['TranslationSet'][0]
            terms = translation['To'].split(' OR ')
            for term in terms:
                if '[MeSH Terms]' in term:
                    term1 = term.replace('[MeSH Terms]', '[Mesh]').strip()
            return term1
        else:
            term1 = input_string
                    
                    # mesh_terms.append(term.replace('[MeSH Terms]', '[Mesh]').strip())
                    # mesh_terms.append(term.replace('[MeSH Terms]', '').replace('"', '').strip())
                    # print('mesh term: ',mesh_terms)
                    # print('term1: ',term1)
                    
            return f'"{term1}"'
    except:
        term1 = input_string
        # print("***except: ", input_string)
    
    # 这个函数可以进行自定义的处理，下面只是一个例子，做小写转换
        return f'"{term1}"'

def remove_bracketed_content(text):
    cleaned = []
    bracket_depth = 0
    for char in text:
        if char == '(':
            bracket_depth += 1
        elif char == ')':
            if bracket_depth > 0:
                bracket_depth -= 1
        else:
            if bracket_depth == 0:
                cleaned.append(char)
    return ''.join(cleaned)

# 定义一个函数对引号内的部分进行处理
def query_opt(string):
   
    # 使用正则表达式找到所有引号内的内容并进行替换
    # string_without_parentheses = re.sub(r'\(.*?\)', '', string)
    # matches = re.findall(r'"[^"]*"', string_without_parentheses)
    
    cleaned_text = remove_bracketed_content(string)
    matches = re.findall(r'"([^"]*)"', cleaned_text)
    #r'"[^"]*"'
    # print("matches: ",matches)
    list1 =  [Mesh(match) for match in matches]
    
    # print('list1: ',list1)
    
    list2 = []
    start = string.find('(')
    if start == -1:
        list2 = []
    else:
        # 初始化括號計數器
        counter = 1
        end = start + 1
        # 遍歷字符直到括號完全閉合
        while end < len(string) and counter > 0:
            if string[end] == '(':
                counter += 1
            elif string[end] == ')':
                counter -= 1
            end += 1
        # 提取完整括號內容
        st2 = string[start:end] if counter == 0 else ""
        list2.append(st2)
    # print('list2: ',list2)
    
    list3 = list1 + list2
    
    # part1 = list1[0]
    # print('$$$$: ',list1[0])
    
    
    
    # 组合query3
    # query3 = f'({query1_clean}) OR ({query2})'
    
    # print(query3)
    
    par = " AND ".join(list3)
    
    par = par + " AND 2022:3000[dp]"
    
    par = replace(par)
    
    first_condition = par.split(' AND ', 1)[0].strip()   # 提取第一个AND前的条件

    # 构建query2
    query2 = f'{first_condition} AND ("nccn clinical practice guidelines" OR "nccn guidelines insights" OR "esmo clinical practice guideline*") AND 2022:3000[dp]'
    
    string = string + " AND 2022:3000[dp]"
    
    query = f'({string}) OR ({par}) OR ({query2})'
    
    # par = par + " AND 2020:3000[dp]"
    # rocessed_string = re.sub(r'"(.*?)"', replace_match, string)
    return query
