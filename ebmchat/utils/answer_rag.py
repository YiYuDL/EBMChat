from langchain_community.vectorstores import FAISS
from langchain.prompts import ChatPromptTemplate
from langchain.schema.output_parser import StrOutputParser
from langchain_community.document_loaders import PyPDFLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from pathlib import Path
from utils.model_runtime import create_chat_model, get_embeddings


def RAG(pmid,question):

    file = str(pmid) + '.pdf'

    source_dir = Path(__file__).resolve().parent.parent / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = source_dir / file

    loader = PyPDFLoader(pdf_path)
    pages = loader.load()
    splitter = RecursiveCharacterTextSplitter(
        chunk_size=1200,
        chunk_overlap=200,
    )
    chunks = splitter.split_documents(pages)
    embeddings = get_embeddings()

    db = FAISS.from_documents(chunks, embeddings)
    retriever = db.as_retriever(search_kwargs={"k": 4})

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

    Literature context:
    {context}

    Based on the uploaded literature, what do the authors say about {question}?
       Answer:
       """

    prompt = ChatPromptTemplate.from_template(template)

    llm = create_chat_model()


    documents = retriever.invoke(question)
    context = "\n\n".join(document.page_content for document in documents)

    prompt_input = {"context": context, "question": question}
    prompt_result = prompt.invoke(prompt_input)

    llm_result = llm.invoke(prompt_result)

    output = StrOutputParser().invoke(llm_result)


    answer = output.strip()
    if not answer:
        raise RuntimeError("The model returned an empty answer.")
    return answer
