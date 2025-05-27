import re
from Bio import Entrez
import pandas as pd

Entrez.email = "YOURemail@gmail.com"

ORGANS = ['liver', 'brain', 'heart', 'lung', 'kidney', 'intestine', 'pancreas']
KEY_TERMS = ['disease', 'syndrome']

def get_pubmed_snippets(gene_symbol, retmax=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=gene_symbol, retmax=retmax)
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids:
            return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
        abstracts = handle.read()
        return parse_relevant_snippets(gene_symbol, abstracts)
    except Exception as e:
        print(f"[ERROR] Failed for {gene_symbol}: {e}")
        return []

def parse_relevant_snippets(gene, text):
    results = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        # Check for organ mentions
        for organ in ORGANS:
            if re.search(rf'\b{organ}\b', line, re.IGNORECASE):
                results.append((gene, "ORGAN", line))
                break
        # Check for disease/syndrome keywords and capture snippet
        for term in KEY_TERMS:
            matches = re.finditer(rf'\b(\w+\s+){{0,3}}{term}(\s+\w+){{0,3}}\b', line, re.IGNORECASE)
            for match in matches:
                snippet = match.group(0).strip()
                results.append((gene, "TERM", snippet))
    return results

gene_list = [
    "H3-2", "MTERF1", "ICAM4", "PF4V1", "CUX1", "POT1", "CAPS2", "FTSJ1",
    "ESPL1", "CEACAM4", "TNFSF11", "SLC4A1", "MACROH2A1", "TSPO2", "DIP2A", "HIP1",
    "CTR9", "C10orf90", "ANKRD36C", "GPR6", "ICE2", "KCNC2", "FBXO24", "CD36", "PRPS2",
    "CIZ1", "IL15RA", "C22orf31", "TRIM58", "ADH7", "KCNJ10", "APLP1", "CFT1", "SLC51A", "SLC17A7"
]

all_snippets = []
for gene in gene_list:
    print(f"[INFO] Processing {gene}")
    snippets = get_pubmed_snippets(gene)
    all_snippets.extend(snippets)

df = pd.DataFrame(all_snippets, columns=["Gene", "Tag", "Snippet"])

output_file = "gene_snippets_utf7.csv"
df.to_csv(output_file, index=False, encoding='utf-7')
