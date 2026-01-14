# Portefólio de Algoritmos em Bioinformática (AASB 2025/2026)

> **Nota:** Este repositório disponibiliza uma **biblioteca Python** com os algoritmos abordados na UC, organizada por módulos temáticos e pensada para ser **facilmente importável e reutilizável por terceiros**.

## Índice
- [Enquadramento](#enquadramento)
- [Funcionalidades](#funcionalidades)
- [Estrutura do projeto](#estrutura-do-projeto)
- [Instalação](#instalação)
- [Utilização](#utilização)
  - [1) Sequências biológicas](#1-sequências-biológicas)
  - [2) Alinhamento de sequências](#2-alinhamento-de-sequências)
  - [3) Motifs e padrões](#3-motifs-e-padrões)
  - [4) BLAST simplificado](#4-blast-simplificado)
  - [5) Análise filogenética](#5-análise-filogenética)
- [Testes unitários e cobertura](#testes-unitários-e-cobertura)
- [Qualidade do código (Radon)](#qualidade-do-código-radon)
- [Documentação (Sphinx)](#documentação-sphinx)
- [Reprodutibilidade](#reprodutibilidade)
- [Autores](#autores)

---

## Enquadramento

Este trabalho prático consiste na implementação de um **portefólio de algoritmos em bioinformática** sob a forma de uma biblioteca Python, com:
- código modular e reutilizável;
- docstrings completas e documentação gerada com Sphinx;
- testes unitários com cobertura mínima;
- verificação de complexidade/qualidade com Radon.

A implementação foi consolidada a partir do notebook `Portfólio.ipynb`, que serviu como base de desenvolvimento e validação manual dos algoritmos.

---

## Funcionalidades

A biblioteca está organizada em cinco áreas principais:

1. **Sequências biológicas**
   - validação/identificação (DNA, RNA, proteína), complemento reverso, transcrição DNA→RNA, contagens simples;

2. **Alinhamento de sequências**
   - matrizes de pontos (dot plots) com janela deslizante;
   - matrizes de substituição (DNA e proteína, via BLOSUM62);
   - Needleman–Wunsch (global): cálculo + reconstrução;
   - Smith–Waterman (local): cálculo + reconstrução;
   - alinhamento progressivo (múltiplo): construção de consenso + alinhamento;

3. **Motifs e padrões**
   - procura de padrões com ambiguidades IUPAC;
   - conversão PROSITE → regex e pesquisa em sequências;
   - simulação de digestão por enzimas de restrição (sítio com `^`, ambiguidades e fragmentação);
   - PWM/PSSM: construção, scoring e subsequência mais provável;

4. **BLAST simplificado**
   - indexação por k-mers, identificação de hits (seeds), extensão sem gaps e seleção do melhor alinhamento;

5. **Análise filogenética**
   - construção de matriz de distâncias (Levenshtein);
   - UPGMA: clustering e registo das junções para posterior construção da árvore.

---

## Estrutura do projeto

A organização recomendada (e utilizada) é a seguinte:

```
bioinf_portfolio/
├── bioinf/                 # package principal
│   ├── sequencias.py
│   ├── alinhamento.py
│   ├── motifs.py
│   ├── blast.py
│   └── filogenia.py
├── tests/                  # testes unitários (pytest)
│   └── test_*.py
├── docs/                   # documentação Sphinx
│   ├── conf.py
│   ├── index.rst
│   └── _build/html/        # saída HTML (gerada)
├── exemplos/               # scripts de demonstração (opcional)
├── README.md
└── requirements.txt
```

---

## Instalação

### Pré-requisitos
- Python **3.10+** (recomendado).
- `pip` e ambiente virtual.

### 1) Criar e ativar ambiente virtual

**Linux/macOS**
```bash
python -m venv .venv
source .venv/bin/activate
```

**Windows (PowerShell)**
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

### 2) Instalar dependências

```bash
pip install -r requirements.txt
```

> Dependência principal (runtime): `blosum` (para BLOSUM62).

---

## Utilização

A biblioteca pode ser importada a partir do package `bioinf`.  
Os exemplos abaixo assumem que está a executar a partir da raiz do repositório (ou que o package está instalado/está no `PYTHONPATH`).

### 1) Sequências biológicas

```python
from bioinf.sequencias import dna_counter, reverse_complement, dna_to_rna, identificar_sequencia

seq = "ATGCCGTA"

print(identificar_sequencia(seq))      # "DNA"
print(dna_counter(seq))                # contagem A/C/G/T
print(reverse_complement(seq))         # complemento reverso
print(dna_to_rna(seq))                 # transcrição DNA → RNA
```

### 2) Alinhamento de sequências

**Dot plot (janela deslizante + stringency)**

```python
from bioinf.alinhamento import imprimir_dot_plot

imprimir_dot_plot("ACGTACGT", "ACGTTGCA", window=3, stringency=2)
```

**Needleman–Wunsch e reconstrução**

```python
from bioinf.alinhamento import needleman_wunsch
from bioinf.alinhamento import reconstruir_alinhamento  # ou de um módulo comum, se aplicável

grid, setas = needleman_wunsch("ACGT", "AGT", match=2, mismatch=-1, space=-2)
a1, a2 = reconstruir_alinhamento(setas, "ACGT", "AGT")
print(a1)
print(a2)
```

**Smith–Waterman e reconstrução**

```python
from bioinf.alinhamento import smith_waterman
from bioinf.alinhamento import reconstruir_alinhamento

grid, setas = smith_waterman("ACGT", "TACG", match=2, mismatch=-1, space=-2)
a1, a2 = reconstruir_alinhamento(setas, "ACGT", "TACG")
print(a1)
print(a2)
```

**Alinhamento múltiplo (progressivo)**

```python
from bioinf.alinhamento import alinhamento_progressivo, consenso

from bioinf.alinhamento import matriz_substituição_dna

subst = matriz_substituição_dna(match=2, mismatch=-1)
alinh = alinhamento_progressivo(["ACGT", "AGT", "ACCT"], subst, space=-2)

print(alinh)
print(consenso(alinh))
```

### 3) Motifs e padrões

**IUPAC → regex e pesquisa por sobreposições**

```python
from bioinf.motifs import iupac_to_regex, find_overlapping

regex = iupac_to_regex("ATGN")     # "ATG[ACGT]"
pos = find_overlapping("ATGACATGTTATGA", regex)
print(pos)
```

**PROSITE → regex**

```python
from bioinf.motifs import prosite_to_regex, find_prosite

pat = prosite_to_regex("C-x(2)-C")       # "C.{2}C"
pos = find_prosite("TTTCAGCAGCTTC", "C-x(2)-C")
print(pos)
```

**Enzimas de restrição: corte e fragmentação**

```python
from bioinf.motifs import digest_dna

cuts, frags = digest_dna("AAGAATTCTTGAATTC", "G^AATTC")
print(cuts)
print(frags)
```

**PWM/PSSM: scoring e melhor subsequência**

```python
from bioinf.motifs import build_pwm, pssm_from_pwm, best_subsequence

kmers = ["ATG", "ATG", "ATA", "ATC"]
pwm = build_pwm(kmers, alphabet="ACGT")
pssm = pssm_from_pwm(pwm)

best = best_subsequence(pssm, "CCCATGAAATA")
print(best)  # (melhor_subseq, posição_0based, score)
```

### 4) BLAST simplificado

```python
from bioinf.blast import alinhamento_pro

res = alinhamento_pro("ACGT", "TACGTG", k=2)
print(res)
# {'query_start': ..., 'subject_start': ..., 'tamanho': ..., 'alinhado_q': ..., 'alinhado_s': ...}
```

### 5) Análise filogenética

```python
from bioinf.filogenia import matriz_distancias, upgma

seqs = ["ACGT", "ACGG", "TCGT"]
D = matriz_distancias(seqs)
ligacoes = upgma(D, labels=seqs.copy())
print(ligacoes)  # lista de junções (label1, label2, distância)
```

---

## Testes unitários e cobertura

Os testes estão no diretório `tests/` e devem cobrir **casos normais, limites e exceções**, garantindo uma cobertura **≥ 80%**.

Executar testes e cobertura:

```bash
pytest -q --cov=bioinf --cov-report=term-missing
```

Gerar relatório HTML de cobertura (opcional):

```bash
pytest --cov=bioinf --cov-report=html
```

---

## Qualidade do código (Radon)

A qualidade/complexidade deve cumprir o requisito mínimo (Radon **≥ B**, eliminatório).  
Comandos típicos:

```bash
radon cc bioinf -a -s
radon mi bioinf -s
```

Para gerar relatórios (ex.: JSON), pode usar:

```bash
radon cc bioinf -s -j > radon_cc.json
radon mi bioinf -s -j > radon_mi.json
```

---

## Documentação (Sphinx)

A documentação é gerada com **Sphinx**, baseada em **docstrings completas** (descrição, parâmetros, retorno, exceções e exemplos).

### Gerar HTML

```bash
cd docs
make html
```

A saída fica disponível em `docs/_build/html/index.html`.

---

## Reprodutibilidade

- Fixar versões no `requirements.txt`.
- Executar todo o pipeline (instalação → testes → Radon → Sphinx) numa máquina limpa antes da entrega.
- Incluir, na entrega final, o HTML da documentação, os relatórios do Radon e os testes.

---

## Autores

- **Aluno/a:** <NOME>
- **N.º:** <NÚMERO>
- **Mestrado:** <CURSO / INSTITUIÇÃO>
- **Ano letivo:** 2025/2026
