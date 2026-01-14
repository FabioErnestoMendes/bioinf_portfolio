def distancia(s1, s2):
    """Calcula a distância de edição (Levenshtein) entre duas strings.

    Operações permitidas (custo 1 cada):
        1) Inserir um carácter
        2) Remover um carácter
        3) Substituir um carácter

    Args:
        s1 (str): Primeira string (ex.: sequência de DNA). Pode ser vazia "".
        s2 (str): Segunda string (ex.: sequência de DNA). Pode ser vazia "".

    Returns:
        int: Número mínimo de operações necessárias para transformar `s1` em `s2`.

    Raises:
        TypeError: Se `s1` ou `s2` não forem strings.

    Example:
        >>> distancia("ATCG", "ATCG")
        0
        >>> distancia("AAA", "TTT")
        3
        >>> distancia("", "A")
        1
    """
    if not isinstance(s1, str) or not isinstance(s2, str):
        raise TypeError("distancia: s1 e s2 devem ser strings.")

    mat = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]

    for i in range(len(s1) + 1):
        for j in range(len(s2) + 1):
            if i == 0:
                mat[i][j] = j
            elif j == 0:
                mat[i][j] = i
            else:
                custo_subst = 0 if s1[i - 1] == s2[j - 1] else 1
                mat[i][j] = min(
                    mat[i - 1][j] + 1,
                    mat[i][j - 1] + 1,
                    mat[i - 1][j - 1] + custo_subst
                )

    return mat[len(s1)][len(s2)]


def matriz_distancias(lista_seqs):
    """Constrói uma matriz de distâncias (simétrica) entre todas as sequências.

    Dada uma lista de sequências, calcula a distância de Levenshtein entre cada par
    e devolve uma matriz NxN onde:

        - dist[i][j] = distância entre lista_seqs[i] e lista_seqs[j]
        - dist é simétrica: dist[i][j] == dist[j][i]
        - diagonal é zero: dist[i][i] == 0

    Args:
        lista_seqs (list[str]): Lista de sequências (strings).

    Returns:
        list[list[int]]: Matriz NxN de distâncias.

    Raises:
        TypeError: Se `lista_seqs` não for uma lista, ou se contiver elementos não-string.

    Example:
        >>> matriz_distancias(["AAA", "AAT", "TTT"])
        [[0, 1, 3],
         [1, 0, 2],
         [3, 2, 0]]
    """
    if not isinstance(lista_seqs, list):
        raise TypeError("matriz_distancias: lista_seqs deve ser uma lista.")
    for s in lista_seqs:
        if not isinstance(s, str):
            raise TypeError("matriz_distancias: todas as sequências devem ser strings.")

    n = len(lista_seqs)
    dist = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = distancia(lista_seqs[i], lista_seqs[j])
            dist[i][j] = d
            dist[j][i] = d

    return dist




def _validar_lista(obj, nome):
    """Garante que obj é list."""
    if not isinstance(obj, list):
        raise TypeError(f"{nome}: deve ser uma lista.")


def _validar_labels(labels):
    """Garante que labels é list[str] e tem pelo menos 2 elementos."""
    _validar_lista(labels, "upgma: labels")
    for l in labels:
        if not isinstance(l, str):
            raise TypeError("upgma: todos os labels devem ser strings.")
    if len(labels) < 2:
        raise ValueError("upgma: são necessárias pelo menos 2 sequências.")


def _validar_matriz_lista(distancias):
    """Garante que distancias é uma lista de listas."""
    _validar_lista(distancias, "upgma: distancias")
    for row in distancias:
        if not isinstance(row, list):
            raise TypeError("upgma: distancias deve ser uma matriz (lista de listas).")


def _validar_matriz_quadrada(distancias, n):
    """Garante que a matriz é NxN."""
    if len(distancias) != n:
        raise ValueError("upgma: matriz de distâncias deve ser NxN.")
    for row in distancias:
        if len(row) != n:
            raise ValueError("upgma: matriz de distâncias deve ser NxN.")


def _validar_matriz_numerica(distancias):
    """Garante que todos os valores da matriz são int ou float."""
    for row in distancias:
        for v in row:
            if not isinstance(v, (int, float)):
                raise TypeError("upgma: distancias deve conter apenas valores numéricos.")


def _encontrar_par_mais_proximo(distancias):
    """Encontra (x, y) com menor distância para x < y.

    Args:
        distancias (list[list[int|float]]): Matriz de distâncias.

    Returns:
        tuple[int, int, float]: (x, y, menor_distância).
    """
    n = len(distancias)
    x, y = 0, 1
    min_dist = float(distancias[0][1])

    i = 0
    while i < n:
        row = distancias[i]
        j = i + 1
        while j < n:
            d = row[j]
            if d < min_dist:
                min_dist = float(d)
                x, y = i, j
            j += 1
        i += 1

    return x, y, min_dist


def _nova_linha_media(distancias, x, y):
    """Calcula a nova linha de distâncias por média aritmética simples.

    Args:
        distancias (list[list[int|float]]): Matriz atual.
        x (int): Índice do cluster x.
        y (int): Índice do cluster y.

    Returns:
        list[float]: Distâncias do novo cluster para os restantes clusters.
    """
    n = len(distancias)
    nova = []
    k = 0
    while k < n:
        if k != x and k != y:
            nova.append((float(distancias[x][k]) + float(distancias[y][k])) / 2.0)
        k += 1
    return nova


def _remover_indices(distancias, clusters, labels, x, y):
    """Remove linhas/colunas e entradas x,y (mutando as estruturas)."""
    for idx in sorted([x, y], reverse=True):
        distancias.pop(idx)
        clusters.pop(idx)
        labels.pop(idx)

    col_to_pop = x if x < y else y
    i = 0
    while i < len(distancias):
        distancias[i].pop(col_to_pop)
        i += 1


def _adicionar_cluster(distancias, clusters, labels, novo_cluster, nova_linha, novo_label):
    """Adiciona novo cluster e atualiza a matriz."""
    distancias.append(nova_linha + [0.0])

    i = 0
    while i < len(distancias) - 1:
        distancias[i].append(nova_linha[i])
        i += 1

    clusters.append(novo_cluster)
    labels.append(novo_label)


def upgma(distancias, labels):
    """Executa o algoritmo UPGMA e devolve as junções efetuadas.

    Passos do algoritmo::

        1) Começa com um cluster por sequência (cada label é um cluster).
        2) Procura o par de clusters com menor distância.
        3) Regista a junção (label1, label2, distância).
        4) Junta clusters e atualiza as distâncias usando média aritmética simples.
        5) Repete até restar um único cluster.

    Args:
        distancias (list[list[int|float]]):
            Matriz NxN de distâncias. distancias[i][j] é a distância entre labels[i] e labels[j].
        labels (list[str]):
            Lista de labels (tamanho N), na mesma ordem das linhas/colunas da matriz.

    Returns:
        list[tuple[str, str, float]]:
            Lista de ligações/junções. Cada elemento é:
                (labelA, labelB, distancia)
            Existem exatamente N-1 junções para N sequências.

    Raises:
        TypeError: Se `labels` não for list[str], ou `distancias` não for matriz numérica.
        ValueError: Se N < 2, ou se `distancias` não for NxN.

    Example:
        >>> seqs = ["AAA", "AAT", "TTT"]
        >>> D = matriz_distancias(seqs)
        >>> lig = upgma(D, seqs.copy())
        >>> len(lig)
        2

    """


    _validar_labels(labels)
    _validar_matriz_lista(distancias)
    _validar_matriz_quadrada(distancias, len(labels))
    _validar_matriz_numerica(distancias)

    clusters = [[i] for i in range(len(labels))]
    ligacoes = []

    while len(clusters) > 1:
        x, y, min_dist = _encontrar_par_mais_proximo(distancias)
        ligacoes.append((labels[x], labels[y], float(min_dist)))

        novo_cluster = clusters[x] + clusters[y]
        nova_linha = _nova_linha_media(distancias, x, y)

        _remover_indices(distancias, clusters, labels, x, y)

        novo_label = "(" + ",".join(ligacoes[-1][:2]) + ")"
        _adicionar_cluster(distancias, clusters, labels, novo_cluster, nova_linha, novo_label)

    return ligacoes
