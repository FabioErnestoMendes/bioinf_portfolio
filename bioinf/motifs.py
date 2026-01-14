# 3 - motifs.py
import re
import math
from collections import Counter

# Padrões com ambiguidades

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]",
}


def _clean_seq(seq):
    if seq is None:
        raise TypeError("seq não pode ser None")
    return seq.upper().strip()


def _find_overlapping_positions(seq, rgx):
    return [m.start() for m in re.finditer("(?=(" + rgx + "))", seq)]


def iupac_to_regex(pat):
    """
    Converte um padrão IUPAC (DNA) para uma expressão regular (regex).

    Args:
        pat (str): Padrão IUPAC (ex.: "CCWGG", "ATGN").

    Returns:
        str: Regex equivalente (ex.: "CC[AT]GG").

    Raises:
        ValueError: Se o padrão estiver vazio ou contiver símbolos IUPAC inválidos.

    Example:
        iupac_to_regex("CCWGG") -> "CC[AT]GG"
    """
    pat = _clean_seq(pat)
    if not pat:
        raise ValueError("Padrão vazio")

    out = []
    for c in pat:
        if c not in IUPAC:
            raise ValueError("IUPAC inválido: " + c)
        out.append(IUPAC[c])
    return "".join(out)


def find_overlapping(seq, pat_iupac):
    """
    Encontra todas as ocorrências (com sobreposição) de um padrão IUPAC numa sequência.

    Args:
        seq (str): Sequência alvo (DNA).
        pat_iupac (str): Padrão IUPAC (DNA).

    Returns:
        list[int]: Lista de posições (base 0) onde começa cada correspondência.

    Raises:
        ValueError: Se o padrão IUPAC estiver vazio ou inválido.

    Example:
        find_overlapping("AAAA", "AA") -> [0, 1, 2]
    """
    seq = _clean_seq(seq)
    rgx = iupac_to_regex(pat_iupac)
    return _find_overlapping_positions(seq, rgx)


# PROSITE -> regex


def prosite_to_regex(prosite):
    """
    Converte um padrão PROSITE (versão mínima) para regex.

    Suporta o básico:
    - Remove '-' (separadores)
    - 'x' ou 'X' -> '.' (qualquer aminoácido)
    - '{ABC}' -> '[^ABC]' (qualquer exceto ABC)
    - '(n)' -> '{n}'
    - '(n,m)' -> '{n,m}'

    Nota:
        Isto NÃO implementa PROSITE completo (ex.: âncoras '<' '>', alguns detalhes avançados).

    Args:
        prosite (str): Padrão PROSITE (ex.: "C-x(2)-C").

    Returns:
        str: Regex equivalente.

    Raises:
        ValueError: Se o padrão estiver vazio.

    Example:
        prosite_to_regex("C-x(2)-C") -> "C.{2}C"
    """
    if prosite is None:
        raise TypeError("prosite não pode ser None")

    p = prosite.strip().replace("-", "")
    if not p:
        raise ValueError("PROSITE vazio")

    p = p.replace("x", ".").replace("X", ".")
    p = re.sub(r"\{([A-Za-z]+)\}", r"[^\1]", p)
    p = re.sub(r"\((\d+)\)", r"{\1}", p)
    p = re.sub(r"\((\d+),(\d+)\)", r"{\1,\2}", p)
    return p


def find_prosite(seq, prosite):
    """
    Encontra todas as ocorrências (com sobreposição) de um padrão PROSITE numa sequência.

    Args:
        seq (str): Sequência alvo (proteína, tipicamente).
        prosite (str): Padrão PROSITE (mínimo).

    Returns:
        list[int]: Lista de posições (base 0) onde começa cada correspondência.

    Raises:
        ValueError: Se o PROSITE estiver vazio.
    """
    seq = _clean_seq(seq)
    rgx = prosite_to_regex(prosite)
    return _find_overlapping_positions(seq, rgx)


# Enzimas de restrição -> regex + cortes + fragmentos


def digest_dna(seq, restriction_site):
    """
    Simula a digestão de uma sequência de DNA por uma enzima de restrição.

    A posição de corte é indicada por '^' dentro do sítio.
    O sítio pode conter ambiguidades IUPAC (ex.: W, N, etc.).

    Args:
        seq (str): Sequência de DNA.
        restriction_site (str): Sítio com '^' (ex.: "G^AATTC", "CC^WGG").

    Returns:
        tuple[list[int], list[str]]:
            - cuts: posições de corte (base 0) na sequência
            - frags: lista de fragmentos resultantes (em ordem)

    Raises:
        ValueError: Se faltar '^', se o motivo ficar vazio, ou se a posição de corte for inválida.

    Example:
        digest_dna("GAATTCC", "G^AATTC") -> ([1], ["G", "AATTCC"])
    """
    seq = _clean_seq(seq)
    site = _clean_seq(restriction_site)

    if "^" not in site:
        raise ValueError("O sítio de restrição deve conter '^'")

    cut_i = site.index("^")
    motif = site.replace("^", "")

    if not motif:
        raise ValueError("Motivo de restrição vazio")

    if cut_i < 0 or cut_i > len(motif):
        raise ValueError("Posição de corte fora do intervalo")

    rgx = iupac_to_regex(motif)
    starts = _find_overlapping_positions(seq, rgx)

    cuts = sorted(set(st + cut_i for st in starts))

    frags = []
    last = 0
    for c in cuts:
        frags.append(seq[last:c])
        last = c
    frags.append(seq[last:])

    return cuts, frags


# PWM e PSSM


def _validate_kmers(kmers, alphabet):
    if kmers is None:
        raise TypeError("kmers não pode ser None")
    if not kmers:
        raise ValueError("Nenhum kmer fornecido")

    k = len(kmers[0])
    if k == 0:
        raise ValueError("kmer vazio")

    if any(len(x) != k for x in kmers):
        raise ValueError("Todos os kmers devem ter o mesmo comprimento")

    allowed = set(alphabet)
    cleaned = []
    for x in kmers:
        x2 = _clean_seq(x)
        if any(c not in allowed for c in x2):
            raise ValueError("Kmer contém símbolos fora do alfabeto: " + x)
        cleaned.append(x2)

    return cleaned, k


def build_pwm(kmers, alphabet="ACGT"):
    """
    Constrói uma PWM (Position Weight Matrix) a partir de uma lista de kmers.

    Args:
        kmers (list[str]): Lista de sequências de igual comprimento (k).
        alphabet (str): Alfabeto permitido (padrão "ACGT").

    Returns:
        list[dict[str, float]]: PWM como lista de colunas, cada coluna {base: prob}.

    Raises:
        ValueError: Se kmers estiver vazio, se k=0, se comprimentos diferirem,
            ou se existirem símbolos fora do alfabeto.
    """
    kmers2, k = _validate_kmers(kmers, alphabet)

    pwm = []
    for i in range(k):
        col = [x[i] for x in kmers2]
        counts = Counter(col)
        total = len(col)
        pwm.append({b: counts.get(b, 0) / total for b in alphabet})

    return pwm


def pssm_from_pwm(pwm, alphabet="ACGT"):
    """
    Converte PWM em PSSM usando log2(p / bg), com bg uniforme (`1/|alphabet|`).

    Args:
        pwm (list[dict[str, float]]): PWM.
        alphabet (str): Alfabeto.

    Returns:
        list[dict[str, float]]: PSSM (scores log2).

    Notes:
        Se p=0, o score fica -inf.
    """
    bg = 1.0 / len(alphabet)
    pssm = []
    for col in pwm:
        d = {}
        for b in alphabet:
            p = col.get(b, 0)
            if p == 0:
                d[b] = float("-inf")
            else:
                d[b] = math.log2(p / bg)
        pssm.append(d)
    return pssm


def score_kmer(pssm, kmer):
    """
    Calcula o score de um kmer dado uma PSSM.

    Args:
        pssm (list[dict[str, float]]): PSSM.
        kmer (str): Sequência com comprimento igual ao número de colunas da PSSM.

    Returns:
        float: Score total.

    Raises:
        ValueError: Se o comprimento do kmer não corresponder à PSSM.
    """
    kmer = _clean_seq(kmer)
    if len(kmer) != len(pssm):
        raise ValueError("O comprimento do kmer deve corresponder ao comprimento da PSSM")

    s = 0.0
    for i, b in enumerate(kmer):
        s += pssm[i].get(b, float("-inf"))
    return s


def best_subsequence(pssm, seq):
    """
    Encontra a subsequência (window) com maior score segundo a PSSM.

    Args:
        pssm (list[dict[str, float]]): PSSM.
        seq (str): Sequência onde procurar.

    Returns:
        tuple[str, int, float]: (melhor_subseq, posição_inicial_0based, melhor_score)

    Raises:
        ValueError: Se a sequência for mais curta do que o comprimento do motivo (len(pssm)).
    """
    seq = _clean_seq(seq)
    k = len(pssm)

    if k == 0:
        raise ValueError("PSSM vazia")
    if len(seq) < k:
        raise ValueError("Sequência menor que o comprimento do motivo")

    best = ""
    best_pos = -1
    best_score = float("-inf")

    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i + k]
        sc = score_kmer(pssm, kmer)
        if sc > best_score:
            best_score = sc
            best = kmer
            best_pos = i

    return best, best_pos, best_score
