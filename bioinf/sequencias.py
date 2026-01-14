def is_dna(dna):
    """
    Verifica se uma sequência é DNA válida (contém apenas A, G, C, T).

    Args:
        dna (str): Sequência a ser verificada.

    Returns:
        None

    Prints:
        "Sequência válida" se todos os caracteres forem A, G, C ou T.
        "Sequência inválida" caso contrário.

    Example:
        >>> is_dna("AGCTAG")
        Sequência válida
        >>> is_dna("AGCTX")
        Sequência inválida
    """
    dna = dna.upper()
    seq = set("AGCT")
    if set(dna).issubset(seq):
        print("Sequência válida")
    else:
        print("Sequência inválida")
    return



def dna_counter(dna):
    """
    Conta a ocorrência de cada nucleotídeo em uma sequência de DNA.

    Args:
        dna (str): Sequência de DNA a ser analisada.

    Returns:
        dict[str, int]: Dicionário com o número de ocorrências de cada nucleotídeo.
                        As chaves são 'A', 'G', 'C' e 'T'.

    Example:
        >>> dna_counter("AGCTAGC")
        {'A': 2, 'G': 2, 'C': 2, 'T': 1}
        >>> dna_counter("TTTAAA")
        {'A': 3, 'G': 0, 'C': 0, 'T': 3}
    """
    dna = dna.upper()
    dic = {
        "A": dna.count("A"),
        "G": dna.count("G"),
        "C": dna.count("C"),
        "T": dna.count("T")
        }
    return dic

def reverse_complement(dna):
    """
    Retorna o complemento reverso de uma sequência de DNA.

    O complemento reverso é obtido invertendo a sequência
    e trocando cada nucleotídeo pelo seu complementar::

        A <-> T
        G <-> C

    Args:
        dna (str): Sequência de DNA a ser convertida.

    Returns:
        str: Sequência complementar reversa do DNA.

    Examples:
        >>> reverse_complement("ATGC")
        'GCAT'
        >>> reverse_complement("GGCTA")
        'TAGCC'
    """
    dna = dna.upper()
    complementar = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    return "".join(complementar[n] for n in dna[::-1])


def dna_to_rna(dna):
    """
    Converte uma sequência de DNA em RNA complementar.

    A conversão segue as regras:
        A -> U
        T -> A
        G -> C
        C -> G

    Args:
        dna (str): Sequência de DNA a ser transcrita.

    Returns:
        str: Sequência de RNA resultante da transcrição.

    Example:
        >>> dna_to_rna("ATGC")
        'UACG'
        >>> dna_to_rna("GGCTA")
        'CCGAU'
    """
    dna = dna.upper()
    transcricao = {
      "A":"U",
      "T":"A",
      "G":"C",
      "C":"G"
      }
    return "".join(transcricao[letra] for letra in dna)

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