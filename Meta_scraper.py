import re
import requests
from bs4 import BeautifulSoup
import Bio
from Bio import Entrez

Entrez.email = "YourEMAIL@gmail.com"

ORGANS = ['liver', 'brain', 'heart', 'lung', 'kidney', 'intestine', 'pancreas']
KEY_TERMS = ['disease', 'syndrome']

def get_pubmed_snippets(gene_symbol, retmax=5):
    handle = Entrez.esearch(db="pubmed", term=gene_symbol, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]
    if not ids:
        return []

    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    return parse_relevant_snippets(abstracts)

def parse_relevant_snippets(text):
    snippets = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        for organ in ORGANS:
            if re.search(rf'\b{organ}\b', line, re.IGNORECASE):
                snippets.append(f"[ORGAN] {line}")
                break

        for term in KEY_TERMS:
            for match in re.finditer(rf'\b(\w+\s+){0,3}{term}(\s+\w+){0,3}\b', line, re.IGNORECASE):
                snippet = match.group(0).strip()
                snippets.append(f"[TERM] {snippet}")

    return snippets

gene_list = ["VSIG10", "KLF4", "IL9"]

for gene in gene_list:
    print(f"\n--- {gene} ---")
    results = get_pubmed_snippets(gene)
    if not results:
        print("[INFO] No relevant snippets found.")
    else:
        for res in results:
            print(res)
