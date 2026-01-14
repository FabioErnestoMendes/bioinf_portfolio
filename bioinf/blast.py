from collections import defaultdict

def novo_indice(sequencia, k=3):
    """
    Cria um índice de k-mers (substrings de tamanho k) de uma sequência.

    O índice devolve, para cada k-mer, a lista das posições (0-based) onde esse
    k-mer começa na sequência.

    Args:
        sequencia (str): Sequência (por exemplo, DNA/proteína) onde será criado o índice.
        k (int, optional): Tamanho do k-mer. Tem de ser > 0. Por omissão é 3.

    Returns:
        collections.defaultdict[list[int]]: Dicionário (defaultdict) que mapeia cada k-mer
        para uma lista de posições onde ocorre.

    Raises:
        ValueError: Se `k` for menor ou igual a 0.

    Example:
        >>> dict(novo_indice("ATAT", k=2))
        {'AT': [0, 2], 'TA': [1]}
    """

    if k <= 0:  #Verificar se k menor que zero. Caso não parar.
        raise ValueError("k tem de ser > 0")
    indice = defaultdict(list)
    for i in range(len(sequencia) - k + 1):
        indice[sequencia[i:i+k]].append(i)
    return indice



def busca_pares(query, subject, k=3):   #Função para encontramos os hits
    """
    Encontra todos os hits (seeds) entre `query` e `subject` com base em k-mers.

    Um hit é um par (i, j) onde:
    - `i` é a posição inicial do k-mer na `query`
    - `j` é a posição inicial do mesmo k-mer no `subject`

    Args:
        query (str): Sequência query (onde é construído o índice de k-mers).
        subject (str): Sequência subject (onde são procurados os k-mers).
        k (int, optional): Tamanho do k-mer (seed). Por omissão é 3.

    Returns:
        list[tuple[int, int]]: Lista de pares (i, j) com todos os hits encontrados.

    Example:
        >>> busca_pares("ATAT", "GATAT", k=2)
        [(0, 1), (2, 1), (1, 2), (0, 3), (2, 3)]
    """

    indice = novo_indice(query, k)
    pares = []      #Guardar todos os hits como pares de coordenadas
    for j in range(len(subject) - k + 1):       #Garantir que a ultima sbstring tem k comrpimento
        seccao = subject[j:j+k]
        for i in indice.get(seccao, ()):
            pares.append((i, j))
    return pares


def estender_alem(query, subject, i, j, k):
    """
    Estende um seed (i, j) para a esquerda e para a direita (sem gaps, match exato).

    A extensão tenta aumentar o alinhamento enquanto os caracteres forem iguais
    em ambas as sequências. O seed inicial tem comprimento `k`.

    Args:
        query (str): Sequência query.
        subject (str): Sequência subject.
        i (int): Posição inicial do seed na query.
        j (int): Posição inicial do seed no subject.
        k (int): Comprimento do seed.

    Returns:
        tuple[int, int, int]: Triplo (start_i, start_j, tamanho), onde:
        - start_i (int): Início do alinhamento na query
        - start_j (int): Início do alinhamento no subject
        - tamanho (int): Comprimento total do alinhamento estendido

    Example:
        >>> estender_alem("ACGT", "TACGTG", 1, 2, 2)
        (0, 1, 4)
    """

    left = 0    #Seguir para a esquerda e serve para contar quantos caracteres iguais vamos conseguir estender para trás
    while (i - left - 1 >= 0 and j - left - 1 >= 0 and      #Garantir que vamos estar dentro dos limites
           query[i - left - 1] == subject[j - left - 1]):
        left += 1

    # direita (começa após o seed)
    right = k      #Usar o tamanho k e nao 0 porque já temos garantidos k caracteres alinhados a partir de (i, j)
    while (i + right < len(query) and j + right < len(subject) and          #Usar o while para estendermos para a direita e garantir que não excedemos a sequencia
           query[i + right] == subject[j + right]):
        right += 1

    start_i = i - left
    start_j = j - left
    tamanho = left + right
    return start_i, start_j, tamanho


def alinhamento_pro(query, subject, k=3):       #Serve para executar o blast e devolver o melhor alinhamento (sem gaps) que encontramos entre a query e subject

    """
    Executa um BLAST simplificado e devolve o melhor alinhamento sem gaps.

    Processo:
    1) Gera seeds (k-mers) comuns entre query e subject.
    2) Para cada seed (i, j), estende para a esquerda/direita com match exato.
    3) Seleciona o alinhamento com maior comprimento.

    Se não houver seeds possíveis (sequências curtas) ou não existirem hits,
    devolve `None`.

    Args:
        query (str): Sequência query.
        subject (str): Sequência subject.
        k (int, optional): Tamanho do seed (k-mer). Por omissão é 3.

    Returns:
        dict[str, object] | None: Dicionário com o melhor alinhamento, ou `None` se
        não houver alinhamento. O dicionário contém:
        - "query_start" (int): índice inicial na query
        - "subject_start" (int): índice inicial no subject
        - "tamanho" (int): comprimento do alinhamento
        - "alinhado_q" (str): segmento alinhado da query
        - "alinhado_s" (str): segmento alinhado do subject

    Example:
        >>> alinhamento_pro("ACGT", "TACGTG", k=2)
        {'query_start': 0, 'subject_start': 1, 'tamanho': 4, 'alinhado_q': 'ACGT', 'alinhado_s': 'ACGT'}
    """

    if len(query) < k or len(subject) < k:
        return None  # não há seeds possiveis

    pares = busca_pares(query, subject, k)  #Obter todos os hits entre duas sequências
    if not pares:
        return None  # sem hits

    melhor = None  # (tamanho, start_i, start_j)

    for i, j in pares:
        start_i, start_j, tamanho = estender_alem(query, subject, i, j, k)
        if melhor is None or tamanho > melhor[0]:
            melhor = (tamanho, start_i, start_j)

    tamanho, start_i, start_j = melhor
    alinhado_q = query[start_i:start_i + tamanho]
    alinhado_s = subject[start_j:start_j + tamanho]

    return {
        "query_start": start_i,
        "subject_start": start_j,
        "tamanho": tamanho,
        "alinhado_q": alinhado_q,
        "alinhado_s": alinhado_s
    }
