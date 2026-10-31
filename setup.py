from setuptools import setup, find_packages

setup(
    name='ebmchat',
    version='1.1.0',
    url='https://github.com/YiYuDL/EBMChat',
    author='Yi Yu',
    packages=find_packages(),
    author_email='yuyi689@gmail.com',
    python_requires='>=3.10',
    install_requires=[
        'biopython>=1.85,<2',
        'faiss-cpu>=1.8,<2',
        'langchain>=0.3.14,<0.4',
        'langchain-community>=0.3.14,<0.4',
        'langchain-openai>=0.3.1,<0.4',
        'langchain-text-splitters>=0.3.5,<0.4',
        'langgraph>=0.2.64,<0.3',
        'numpy>=1.26,<3',
        'openai>=1.59,<2',
        'pandas>=2.2,<3',
        'pypdf>=5.1,<6',
        'requests>=2.32,<3',
        'scipy>=1.15,<2',
        'selenium>=4.28,<5',
        'socksio>=1,<2',
        'typing-extensions>=4.12,<5',
    ],
)
