## Portefólio - AASB 2025/2026
Este portefólio tem o intuito de apresentar uma **uma bibliotéca Python funcional, documentada e testada** com os algoritmos abordados ao longo da UC. Este codigo foi criado para ser facilmente importável e reutilizável por terceiros.

## Índice
- [Introdução](#Introdução)

- [Funcionalidades](#funcionalidades)

- [Organização](#Organização)

- [Instalação / Packages](#Instalação-/-Packages)

- [Manual do utilizador](#manual-do-utilizador)
  - [1) Sequências biológicas](#1-sequências-biológicas)
  - [2) Alinhamento de sequências](#2-alinhamento-de-sequências)
  - [3) Motifs e padrões](#3-motifs-e-padrões)
  - [4) BLAST simplificado](#4-blast-simplificado)
  - [5) Análise filogenética](#5-análise-filogenética)

- [Validação e testes](#Validação-e-testes)

- [Qualidade e complexidade (Radon)](#Qualidade-e-complexidade-(Radon))

- [Documentação (Sphinx)](#documentação-sphinx)

- [Reprodutibilidade](#reprodutibilidade)

- [Autores](#autores)


## Introdução
O presente trabalho prático tem como objetivo o desenvolvimento de um **portefólio de algoritmos em Bioinformática**, disponibilizado sob a forma de uma **biblioteca em Python**. A solução foi concebida com enfoque na **modularidade** e **reutilização**, garantindo:

- Organização por módulos
- Documentação com auxilio de **docstrings** e **Sphinx**
- Testes unitários com avaliação de cobertura
- Avaliação de qualidade via **Radon**


## Funcionalidades 
A biblioteca seguio a estrutura disponibilizada no enunciado sendo a mesma organizada em 5 modulos.

1. **Sequências biológicas**
    - Validação (DNA/RNA/proteínas), complemento/reverso, transcrição, etc

2. **Alinhamento de Sequências**
    - Matrizes de pontos e de substituição (BLOSUM/PAM)
    - Needleman-Wunsch (global): cálculo + reconstrução
    - Smith-Waterman (local): cálculo + reconstrução
    - Alinhamento progressivo (múltiplo): consenso + alinhamento

3. **Motifs e Padrõe**
    - Procura de padrões fixos com ambiguidades
    - Conversão PROSITE → regex
    - Conversão enzimas restrição → regex + posições corte + fragmentação
    - PWM e PSSM: construção, probabilidade sequência, subsequência mais provável

4. **Blast Simplificado**
    - Query map, identificação de hits, extensão, melhor alinhamento

5. **Análise Filogenética**
    - Matriz de distâncias
    - UPGMA: clustering + construção árvore


## Organização
Estrutura da biblioteca apresentada:

```text
bioinf_portfolio/
├── bioinf/                 
│   ├── sequencias.py
│   ├── alinhamento.py
│   ├── motifs.py
│   ├── blast.py
│   └── filogenia.py
├── tests/                  
│   └── test_*.py
├── docs/                   
│   ├── conf.py
│   ├── index.rst
│   └── _build/html/        
├── exemplos/              
├── README.md
└── requirements.txt 
```

## Instalação / Packages
Esta biblioteca foi concebida para ser executada num ambiente Python, com as dependências definidas em `requirements.txt`. :contentReference[oaicite:4]{index=4}

### Pré-requisitos
- **Python 3.10 +** e 'pip' disponíveis no sistema.
```python
- python -m pip install -r requirements.txt
```

## Manual do Utilizador
A biblioteca pode ser importada a partir do package `bioinf`.

### Sequências biológicas
Funções principais (ver docstrings/Sphinx):
- `identificar_sequencia(seq: str) -> str`
- `reverse_complement(seq: str) -> str`
- `dna_to_rna(seq: str) -> str`
- `dna_counter(seq: str) -> dict[str, int]`


### Alinhamento de sequencias
Funcionalidades principais:
- dot plot com janela deslizante e stringency;
- matrizes de substituição (DNA e proteína);
- alinhamento global (Needleman–Wunsch) e reconstrução;
- alinhamento local (Smith–Waterman) e reconstrução;
- alinhamento múltiplo progressivo e consenso.

### Motifs e padrões
Funcionalidades principais:
- pesquisa com ambiguidades IUPAC;
- conversão PROSITE → expressão regular e pesquisa;
- digestão por enzimas de restrição (sítio com `^`);
- construção e scoring de PWM/PSSM.

### BLAST simplificado
Funcionalidades principais:
- indexação por k-mers (seeds);
- extensão sem gaps;
- seleção do melhor alinhamento obtido.

### Análise filogenética
Funcionalidades principais:
- cálculo de matriz de distâncias (Levenshtein);
- clustering UPGMA e registo das junções para construção de árvore.


## Validação e testes
Os testes encontram-se no diretório `tests/` e são executados com o auxilio do **pytest**.

### Dependências
- `pytest>=8.0`
- `pytest-cov>=5.0`
- `coverage>=7.0`

### Documentação (Sphinx HTML)
- `Sphinx >= 7.0`
- `Sphinx-rtd-theme >= 2 .0`

### Métrica radon
- `radon >= 6 .0`

### Execução dos testes
Deve ser feita a partir da raiz do projeto
```bash
pytest
```


## Qualidade e complexidade (Radon)
Avaliação da complexidade do códgio através do **radon**

### Dependências
- `radon>=6.0`
```bash
radon cc bioinf -s -a
```
###
Nota: Trocar o `bioinf` pelo nome da pasta onde está o código, caso esta seja diferente.


## Documentação (Sphinx)
Para o projeto a documentação foi produzida através do **Sphinx**, a partir das **docstrings** do código.
### Dependências
- `sphinx>=7.0`
- `sphinx-rtd-theme>=2.0`

### Documentação HTML
```bash
cd docs
make html
```


## Reprodutibilidade 
É necessário alguns requesitos para que os resultados possam ser reproduzidos em diferentes maquinas

- **Fixar dependeências** no ficheiro `requirements.txt`
- **Executar o pipeline completo** a partir da raiz do projeto
    - Instalar dependências
    - Execução dos testes
    - Radon
    - Sphinx

## Autores

- **Alunos:**  Ana Rocha PG58813 ; Ângela Sousa PG42489 ; Bruno Ferreira PG58814 ; Fábio Mendes PG61008
- **Mestrado:** \<Bioinformática / Universidade do Minho\>  
- **Ano letivo:** 2025/2026  