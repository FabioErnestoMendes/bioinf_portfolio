#Alinhamento de Sequências

#Matrizes de Pontos:
def identificar_sequencia(seq):
  """
    Identifica se a sequência é DNA, RNA ou proteína.

    Args:
        seq (str): Sequência biológica.

    Returns:
        str: "DNA", "RNA", "AMINO" ou "ERRO" se não for possível identificar.

    Example:
        >>> identificar_sequencia("ACGT")
        'DNA'
  """
  seq = seq.upper()
  dna = set("ACTG")
  rna = set("ACUG")
  proteina = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

  if set(seq).issubset({'A','C','G'}):
      return "Introduza mais informação da sequência para poder determinar com mais exatidão"
  elif set(seq).issubset(dna):
      return "DNA"
  elif set(seq).issubset(rna):
      return "RNA"
  elif set(seq).issubset(proteina):
      return "AMINO"
  else:
      return "ERRO"

def matriz_de_zeros(seq1: str, seq2: str):
    """
    Inicia uma matriz de zeros com tamanho len(seq1)+1 x len(seq2)+1.

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.

    Returns:
        list[list[int]]: Matriz preenchida com zeros.

    Example:
        >>> matriz_de_zeros("AC", "GT")
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    """
    return [[0 for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]


def janela_match(seq1: str, seq2: str, pos1: int, pos2: int, window: int):
    """
    Compara duas janelas de comprimento `window` nas sequências.

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        pos1 (int): Posição inicial na primeira sequência.
        pos2 (int): Posição inicial na segunda sequência.
        window (int): Tamanho da janela a comparar.

    Returns:
        int: Número de posições iguais dentro da janela.

    Example:
        >>> janela_match("ACGT", "ACGA", 0, 0, 3)
        3
    """
    iguais = 0
    for i in range(window):
        if seq1[pos1 + i] == seq2[pos2 + i]:
            iguais += 1
    return iguais


def dot_plot_janela(seq1: str, seq2: str, window: int, stringency: int):
    """
    Cria dot plot usando janela deslizante e stringency mínima.

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        window (int): Tamanho da janela.
        stringency (int): Número mínimo de matches na janela.

    Returns:
        list[list[int]]: Matriz de 0s e 1s indicando matches.

    Example:
        >>> dot_plot_janela("ACGT", "ACGA", 3, 2)
        [[1, 0], [0, 0]]
    """
    matriz = matriz_de_zeros(seq1, seq2)
    limite1, limite2 = len(seq1) - window + 1, len(seq2) - window + 1

    for p1 in range(limite1):
        for p2 in range(limite2):
            if janela_match(seq1, seq2, p1, p2, window) >= stringency:
                matriz[p1][p2] = 1
    return matriz


def imprimir_dot_plot(seq1: str, seq2: str, window: int, stringency: int):
    """
    Imprime um dot plot com sliding window e stringency,
    mostrando explicitamente as janelas de ambas as
    sequências nos eixos e alinhando os pontos corretamente.
    """

    matriz = dot_plot_janela(seq1, seq2, window, stringency)

    limite1 = len(seq1) - window + 1
    limite2 = len(seq2) - window + 1

    janelas1 = [seq1[i:i + window] for i in range(limite1)]
    janelas2 = [seq2[j:j + window] for j in range(limite2)]

    col_width = window

    print(f"Sequência 1: {seq1}")
    print(f"Sequência 2: {seq2}\n")

    header = " " * (col_width + 2)
    for j in janelas2:
        header += j.ljust(col_width + 1)
    print(header)

    for janela, linha in zip(janelas1, matriz):
        row = janela.ljust(col_width) + "  "
        for v in linha:
            row += ("•".ljust(col_width + 1) if v == 1 else " " * (col_width + 1))
        print(row)


#Needleman-Wunsch (global) e Smith-Waterman (local):
# cálculo + reconstrução

def matriz_substituição_proteína():
    """
    Cria a matriz de substituição BLOSUM62 para proteínas.

    Returns:
        dict[tuple[str, str], int]: Matriz de substituição para cada par de aminoácidos.

    Example:
        >>> subst = matriz_substituição_proteína()
        >>> subst[("A", "A")]
        4
    """
    import blosum as bl
    matriz = bl.BLOSUM(62)
    return {(a, b): matriz[a][b] for a in matriz for b in matriz[a]}


def matriz_substituição_dna(match: int, mismatch: int):
    """
    Cria uma matriz de substituição para DNA.

    Args:
        match (int): Pontuação para match.
        mismatch (int): Pontuação para mismatch.

    Returns:
        dict[str, dict[str, int]]: Matriz de substituição.

    Example:
        >>> matriz_substituição_dna(2, -1)["A"]["C"]
        -1
    """
    return {
        "A": {"A": match, "G": mismatch, "C": mismatch, "T": mismatch},
        "C": {"A": mismatch, "G": mismatch, "C": match, "T": mismatch},
        "G": {"A": mismatch, "G": match, "C": mismatch, "T": mismatch},
        "T": {"A": mismatch, "G": mismatch, "C": mismatch, "T": match}
  }


def primeira_linha_e_coluna(seq1: str, seq2: str, space: int):
    """
    Atribui os valores à primeira linha e à primeira coluna do alinhamento global.

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        space (int): Penalidade de gap.

    Returns:
        tuple: (matriz_scores, matriz_setas)
            - matriz_scores: lista de listas com scores iniciais.
            - matriz_setas: lista de listas com setas ("↑", "←") na primeira linha e coluna.

    Example:
        >>> matriz, setas = primeira_linha_e_coluna("AC", "GT", -2)
        >>> matriz
        [[0, -2, -4], [-2, 0, 0], [-4, 0, 0]]
    """
    matriz = matriz_de_zeros(seq1, seq2)
    setas = [["" for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]

    for p2 in range(len(seq2)):
        matriz[0][p2 + 1] = matriz[0][p2] + space
        setas[0][p2 + 1] = "←"

    for p1 in range(len(seq1)):
        matriz[p1 + 1][0] = matriz[p1][0] + space
        setas[p1 + 1][0] = "↑"

    return matriz, setas


def escolha_de_matriz(seq1: str, seq2: str, match: int, mismatch: int):
    """
    Determina a matriz de substituição mais adequada (DNA ou proteína).

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        match (int): Pontuação para match (DNA).
        mismatch (int): Pontuação para mismatch (DNA).

    Returns:
        dict: Matriz de substituição adequada.

    Raises:
        ValueError: Se os tipos das sequências forem incompatíveis.

    Example:
        >>> escolha_de_matriz("ACG", "AGT", 2, -1)
        {'A': {'A': 2, 'G': -1, 'C': -1, 'T': -1}, ...}
    """
    if identificar_sequencia(seq1) == "DNA" and identificar_sequencia(seq2) == "DNA":
        return matriz_substituição_dna(match, mismatch)
    if identificar_sequencia(seq1) == "AMINO" and identificar_sequencia(seq2) == "AMINO":
        return matriz_substituição_proteína()
    raise ValueError("Tipos de sequência incompatíveis")


def melhor_movimento(diag: int, cima: int, esquerda: int) -> tuple[int, str]:
    """
    Determina o melhor movimento entre diagonal, cima ou esquerda.

    Args:
        diag (int): Score da diagonal.
        cima (int): Score de cima.
        esquerda (int): Score da esquerda.

    Returns:
        tuple[int, str]: Melhor score e seta correspondente ("↖", "↑", "←").

    Example:
        >>> melhor_movimento(3, 2, 1)
        (3, '↖')
    """
    melhor = max(diag, cima, esquerda)
    if melhor == diag:
        return melhor, "↖"
    if melhor == cima:
        return melhor, "↑"
    return melhor, "←"


def montar_grid(matriz, seq1, seq2):
    """
    Monta um grid de scores com a primeira linha e coluna mostrando as sequências.

    Args:
        matriz (list[list[int]]): Matriz de scores (inicial ou calculada).
        seq1 (str): Primeira sequência (linha).
        seq2 (str): Segunda sequência (coluna).

    Returns:
        list[list[str|int]]: Grid de scores com seq2 na primeira linha e seq1 na primeira coluna.

    Example:
        >>> matriz = [[0, -2, -4], [-2, 0, 0], [-4, 0, 0]]
        >>> montar_grid(matriz, "AC", "GT")
        [[' ', 'G', 'T'], ['A', -2, -4], ['C', 0, 0]]
    """
    grid = [[" "] + list(seq2)]
    for i, a in enumerate(seq1):
        linha = [a] + matriz[i + 1]
        grid.append(linha)
    return grid


def imprimir_grid(grid, setas):
    """
    Imprime o grid de scores e a matriz de setas de forma legível.

    Args:
        grid (list[list[str|int]]): Grid de scores com seq2 na primeira linha e seq1 na primeira coluna.
        setas (list[list[str]]): Matriz de setas ("↖", "↑", "←", "STOP") correspondente ao grid.

    Returns:
        None

    Example:
        >>> grid = montar_grid(matriz, "AC", "GT")
        >>> imprimir_grid(grid, setas)
        Matriz de Scores:
        [' ', 'G', 'T']
        ['A', -2, -4]
        ['C', 0, 0]

        Matriz de Setas:
        ['', '←', '←']
        ['↑', '↖', '↖']
        ['↑', '↑', '↖']
    """
    print("Matriz de Scores:")
    for linha in grid:
        print(linha)
    print("\nMatriz de Setas:")
    for linha in setas:
        print(linha)


def needleman_wunsch(seq1: str, seq2: str, match: int = 2, mismatch: int = -3, space: int = -4):
    """
    Realiza um alinhamento global usando o algoritmo Needleman-Wunsch.

    Imprime a matriz de scores e a matriz de setas em formato legível (grid).

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        match (int, optional): Pontuação para match. Default é 2.
        mismatch (int, optional): Pontuação para mismatch. Default é -3.
        space (int, optional): Penalidade de gap. Default é -4.

    Returns:
        tuple[list[list[str|int]], list[list[str]]]:
            - grid_scores: Matriz de scores com a primeira linha e coluna mostrando seq2 e seq1.
            - setas: Matriz de setas ("↖", "↑", "←") indicando o caminho do alinhamento.

    Example:
        >>> grid, setas = needleman_wunsch("ACG", "AG")
    """
    matriz, setas = primeira_linha_e_coluna(seq1, seq2, space)
    subst = escolha_de_matriz(seq1, seq2, match, mismatch)

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            diag = matriz[i][j] + (subst.get((a, b), subst.get(a, {}).get(b, 0)))
            cima = matriz[i][j + 1] + space
            esquerda = matriz[i + 1][j] + space
            score, seta = melhor_movimento(diag, cima, esquerda)
            matriz[i + 1][j + 1] = score
            setas[i + 1][j + 1] = seta

    grid = montar_grid(matriz, seq1, seq2)
    imprimir_grid(grid, setas)

    return grid,setas


def smith_waterman(seq1: str, seq2: str, match: int = 2, mismatch: int = -3, space: int = -4):
    """
    Realiza um alinhamento local usando o algoritmo Smith-Waterman.

    Imprime a matriz de scores e a matriz de setas em formato legível (grid).

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        match (int, optional): Pontuação para match. Default é 2.
        mismatch (int, optional): Pontuação para mismatch. Default é -3.
        space (int, optional): Penalidade de gap. Default é -4.

    Returns:
        tuple[list[list[str|int]], list[list[str]]]:
            - grid_scores: Matriz de scores com a primeira linha e coluna mostrando seq2 e seq1.
            - setas: Matriz de setas ("↖", "↑", "←", "STOP") indicando o caminho do alinhamento.

    Example:
        >>> grid, setas = smith_waterman("ACG", "AG")
    """
    matriz = matriz_de_zeros(seq1, seq2)
    setas = [["" for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    subst = escolha_de_matriz(seq1, seq2, match, mismatch)

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            diag = matriz[i][j] + (subst.get((a, b), subst.get(a, {}).get(b, 0)))
            cima = matriz[i][j + 1] + space
            esquerda = matriz[i + 1][j] + space
            melhor_score = max(0, diag, cima, esquerda)

            if melhor_score == 0:
                setas[i + 1][j + 1] = "STOP"
            else:
                setas[i + 1][j + 1] = "↖" if melhor_score == diag else "↑" if melhor_score == cima else "←"

            matriz[i + 1][j + 1] = melhor_score

    grid = montar_grid(matriz, seq1, seq2)
    imprimir_grid(grid, setas)

    return grid,setas


def reconstruir_alinhamento(setas: list[list[str]], seq1: str, seq2: str):
    """
    Reconstrói o alinhamento final a partir da matriz de setas.

    Funciona tanto para Needleman-Wunsch quanto para Smith-Waterman.
    Para Smith-Waterman, o caminho termina ao encontrar "STOP".

    Args:
        setas (list[list[str]]): Matriz de setas ("↖", "↑", "←", "STOP").
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.

    Returns:
        tuple[str, str]: Sequências alinhadas.

    Example:
        >>> grid, setas = needleman_wunsch("ACG", "AG")
        >>> reconstruir_alinhamento(setas, "ACG", "AG")
        ('ACG', 'A-G')
    """
    i, j = len(seq1), len(seq2)
    a1, a2 = [], []

    while i > 0 or j > 0:
        direcao = setas[i][j]

        if direcao == "STOP":
            break
        elif direcao == "↖":
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direcao == "↑":
            a1.append(seq1[i - 1])
            a2.append("-")
            i -= 1
        elif direcao == "←":
            a1.append("-")
            a2.append(seq2[j - 1])
            j -= 1
        else:
            break

    return "".join(reversed(a1)), "".join(reversed(a2))


#Alinhamento progressivo (múltiplo):
# consenso + alinhamento

def consenso(alinhamento: list[str]):
    """
    Calcula a sequência consenso de um alinhamento múltiplo.

    Args:
        alinhamento (list[str]): Lista de sequências alinhadas.

    Returns:
        str: Sequência consenso.

    Example:
        >>> consenso(["ACG", "A-G", "ACG"])
        'ACG'
    """
    seq_consenso = ""
    for coluna in zip(*alinhamento):
        bases = {}
        for simbolo in coluna:
            if simbolo != "-":
                bases[simbolo] = bases.get(simbolo, 0) + 1
        seq_consenso += max(bases, key=bases.get)
    return seq_consenso


def melhor_movimento(diag: int, cima: int, esquerda: int) -> tuple[int, str]:
    """
    Determina o melhor movimento entre diagonal, cima ou esquerda.

    Args:
        diag (int): Score da diagonal.
        cima (int): Score de cima.
        esquerda (int): Score da esquerda.

    Returns:
        tuple[int, str]: Melhor score e seta correspondente ("↖", "↑", "←").

    Example:
        >>> melhor_movimento(3, 2, 1)
        (3, '↖')
    """
    melhor = max(diag, cima, esquerda)
    if melhor == diag:
        return melhor, "↖"
    if melhor == cima:
        return melhor, "↑"
    return melhor, "←"


def alinhar_par(seq1: str, seq2: str, subst: dict, space: int):
    """
    Realiza um alinhamento simples entre duas sequências.

    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        subst (dict): Matriz de substituição.
        space (int): Penalidade de gap.

    Returns:
        tuple[str, str]: Sequências alinhadas.

    Example:
        >>> subst = matriz_substituição_dna(2, -1)
        >>> alinhar_par("ACG", "AG", subst, -2)
        ('ACG', 'A-G')
    """
    matriz = [[0 for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    setas = [["" for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]

    for i in range(len(seq1)):
        matriz[i + 1][0] = matriz[i][0] + space
        setas[i + 1][0] = "↑"

    for j in range(len(seq2)):
        matriz[0][j + 1] = matriz[0][j] + space
        setas[0][j + 1] = "←"

    for i, a in enumerate(seq1):
        for j, b in enumerate(seq2):
            diag = matriz[i][j] + subst[a][b]
            cima = matriz[i][j + 1] + space
            esquerda = matriz[i + 1][j] + space

            valor, seta = melhor_movimento(diag, cima, esquerda)
            matriz[i + 1][j + 1] = valor
            setas[i + 1][j + 1] = seta

    return reconstruir_alinhamento(setas, seq1, seq2)


def alinhar_consenso(alinhamento: list[str], nova_seq: str, subst: dict, space: int):
    """
    Alinha uma nova sequência ao consenso existente.

    Args:
        alinhamento (list[str]): Lista de sequências alinhadas.
        nova_seq (str): Nova sequência a ser alinhada.
        subst (dict): Matriz de substituição.
        space (int): Penalidade de gap.

    Returns:
        list[str]: Novo alinhamento múltiplo.

    Example:
        >>> subst = matriz_substituição_dna(2, -1)
        >>> alinhar_consenso(["ACG", "A-C"], "AG", subst, -2)
        ['ACG', 'A-C', 'A-G']
    """
    cons = consenso(alinhamento)
    cons_al, nova_al = alinhar_par(cons, nova_seq, subst, space)

    novo_alinhamento = []
    idx = 0
    for simbolo in cons_al:
        if simbolo == "-":
            novo_alinhamento.append("-" * len(alinhamento))
        else:
            novo_alinhamento.append("".join(seq[idx] for seq in alinhamento))
            idx += 1
    alinhamento_final = ["".join(coluna) for coluna in zip(*novo_alinhamento)]
    alinhamento_final.append(nova_al)
    return alinhamento_final


def alinhamento_progressivo(seqs: list[str], subst: dict, space: int = -4):
    """
    Realiza um alinhamento múltiplo progressivo de sequências.

    Args:
        seqs (list[str]): Lista de sequências.
        subst (dict): Matriz de substituição.
        space (int, optional): Penalidade de gap. Default é -4.

    Returns:
        list[str]: Lista de sequências alinhadas.

    Example:
        >>> seqs = ["ACG", "AC", "AG"]
        >>> subst = matriz_substituição_dna(2, -1)
        >>> alinhamento_progressivo(seqs, subst)
        ['ACG', 'A-C', 'A-G']
    """
    alinhamento = [seqs[0]]
    for seq in seqs[1:]:
        alinhamento = alinhar_consenso(alinhamento, seq, subst, space)
    return alinhamento