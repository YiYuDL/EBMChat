import os
from langchain_community.vectorstores import FAISS
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain.prompts import ChatPromptTemplate # 聊天提示模板
from langchain.schema.runnable import RunnablePassthrough
from langchain.schema.output_parser import StrOutputParser # 输出解析器
from langchain_community.document_loaders import PyPDFLoader
from langchain.chat_models import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from pathlib import Path
from IPython.display import display

os.environ["OPENAI_API_KEY"] = os.getenv('OPENAI_API_KEY')

# pdf_file = '/home/user/source/36307056.pdf'


def RAG(pmid,question):
    
    file = str(pmid) + '.pdf'
    
    SOURCE_DIR = Path("/home/user/source").mkdir(parents=True, exist_ok=True) or Path("/home/user/source")
    pdf_path = SOURCE_DIR / file
    
    loader = PyPDFLoader(pdf_path)
    pages = loader.load_and_split() 
    embeddings = OpenAIEmbeddings(model="text-embedding-3-large")
    
    db = FAISS.from_documents(pages, embeddings)
    retriever=db.as_retriever()
    
    template = """
        You are a dedicated AI medical assistant specializing in oncology. Based on the provided oncology literatures, provide a detailed and truthful response that addresses the specifics of the question. 
       
    Ensure your responses are:
    
    1.Directly linked to the literature, citing the chapter names of the specific sections when possible (page number is not needed), the format should be: (Section: xxx). 
    2.Presented in concise bullet points, do not bloat the answer unnecessarily. Keep it short and professional. For example, If the question is about treatment, don't include too much information about diagnosis and prognosis. If the question is about diagnosis or prognosis, follow the same principle.
    3. Honest, especially when the answer isn't available in the literatures.
    4. Comprehensive, including all relevant details about disease type, disease stage, patient age, treatment outcome.
    5. List **key points** dynamically based on the content, summarize the same treatment or diagnosis method in one point, and try to be within 4 or 5 points.
       for instance, when asked about treatment, include all the possible treatment options, like first-line, second-line,  etc; When asked about diagnosis, include imaging staging and pathological confirmation, etc; When asked about prognosis, include survival rate and risk of local recurrence and distant metastasis and survival rate, as well as their influencing factors.
       Format:  
       `Number + Colon + Brief Description + Hyphen + Detailed Explanation`.  
       Example:  
       `1. Key Point - Detailed explanation of the point.`
       the detail explanation can include its role, advantages or disadvantages briefly, strictly following the information in literatures
    6. The targeted treatments for different mutations should be illustrated in detail, they can be show as sub-point, like 1.1,1.2, 2.1, 2.2, etc.
    7. Add a **Summary** section at the end, condensing the overall content into **≤50 words**.
    
    Based on the uploaded literature, what do the authors say about {question}?
       Answer:
       """
    
    prompt = ChatPromptTemplate.from_template(template)
    
    llm = ChatOpenAI(model="gpt-4.1", temperature = 0.2)
    
    # rag_chain = (
    #         {"context": retriever, "question": RunnablePassthrough()} # 上下文信息
    #         | prompt
    #         | llm
    #         | StrOutputParser()
    # )
    
    context = retriever.invoke(question)
    
    prompt_input = {"context": context, "question": question}
    prompt_result = prompt.invoke(prompt_input)
    
    llm_result = llm.invoke(prompt_result)
    
    output = StrOutputParser().invoke(llm_result)
    
    # output = output + '\n\n' + 'Citation: ' + citation
    
    # print("output:",output)
    
    return output
# print("output:",rag_chain.invoke(question) )

# Question = 'what are recommended treatments for metastatic colorectal cancer'
# Citation = 'Metastatic colorectal cancer: ESMO Clinical Practice Guideline for diagnosis, treatment and follow-up. Guideline, Ann Oncol, 2023, https://doi.org/10.1016/j.annonc.2022.10.003, 36307056.'
# Pmid = 36307056

# aa = RAG(Pmid,Question)