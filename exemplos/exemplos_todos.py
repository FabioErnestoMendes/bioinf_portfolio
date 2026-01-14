
##Exemplos de alinhamento
from bioinf.alinhamento import *

# Matriz de pontos

seq1 = "ACGT"
seq2 = "ACGA"

print("DOT PLOT:")
imprimir_dot_plot(seq1, seq2, window=3, stringency=2)

# Alinhamento global (Needleman-Wunsch)

print("\nNEEDLEMAN-WUNSCH:")
grid, setas = needleman_wunsch("ACG", "AG")
a1, a2 = reconstruir_alinhamento(setas, "ACG", "AG")
print(a1)
print(a2)

# Alinhamento local (Smith-Waterman)

print("\nSMITH-WATERMAN:")
grid, setas = smith_waterman("ACGT", "CG")
a1, a2 = reconstruir_alinhamento(setas, "ACGT", "CG")
print(a1)
print(a2)

# Matriz de substituição DNA

subst = matriz_substituição_dna(2, -1)
print("\nMATRIZ DNA A vs C:", subst["A"]["C"])

# Alinhamento par a par

print("\nALINHAR PAR:")
a1, a2 = alinhar_par("ACG", "AG", subst, -2)
print(a1)
print(a2)

# Consenso

print("\nCONSENSO:")
print(consenso(["ACG", "A-G", "ACG"]))

# Alinhamento progressivo

print("\nALINHAMENTO MÚLTIPLO:")
seqs = ["ACG", "AC", "AG"]
alinhamento = alinhamento_progressivo(seqs, subst)
for s in alinhamento:
    print(s)

## Exemplos de blast 

from bioinf.blast import *

def teste(query, subject, k=3):
    print("=" * 50)
    print(f"query   = {query}")
    print(f"subject = {subject}")
    print(f"k = {k}")
    print("resultado:", alinhamento_pro(query, subject, k))

# Funcionamento basico
print("INDICE (k=2) para 'ACGT':", novo_indice("ACGT", k=2))
print("PARES (k=2) entre 'ACGT' e 'TTACGTAA':", busca_pares("ACGT", "TTACGTAA", k=2))

# 1) Match perfeito 
teste("ACGT", "TTACGTAA", k=3)

# 2) Query maior, mas o melhor match é parcial 
teste("ACGTT", "GGACGTA", k=3)

# 3) Repetições: Escolher o maior
teste("AAAAC", "TAAAAGAAAAC", k=3)

# 4) Sem hits nenhuns logo vai devolver None
teste("ACGT", "TTTTTT", k=3)

# 5) Sequência menor do que k vai devolver None
teste("AC", "ACGT", k=3)

# 6) Exemplo com k=2 
teste("GATTACA", "TTAC", k=2)

# 7) k inválido (<=0) vai dar erro (ValueError)
try:
    teste("ACGT", "TTACGTAA", k=0)
except ValueError as e:
    print("Erro esperado com k=0:", e)




##Exemplos de filogenia 

from bioinf.filogenia import *

if __name__ == "__main__":
    seqs = ["AAA", "AAT", "ATT", "TTT"]
    labels = seqs.copy()

    print("Distâncias individuais:")
    print("AAA vs AAT =", distancia("AAA", "AAT"))
    print("AAA vs ATT =", distancia("AAA", "ATT"))
    print("AAA vs TTT =", distancia("AAA", "TTT"))
    print()

    print("Matriz de distâncias:")
    D = matriz_distancias(seqs)
    for linha in D:
        print(linha)
    print()

    print("UPGMA (ligações):")
    ligacoes = upgma(D, labels)
    for lig in ligacoes:
        print(lig)


##Exemplos de motifs
from bioinf.motifs import *


print("IUPAC → REGEX:")
padrao = "CCWGG"
regex = iupac_to_regex(padrao)
print("Padrão IUPAC:", padrao)
print("Regex:", regex)

print("\nBUSCA IUPAC:")
seq = "ACCCAGGCCAGG"
posicoes = find_overlapping(seq, "CCWGG")
print("Sequência:", seq)
print("Ocorrências (0-based):", posicoes)

print("\nPROSITE → REGEX:")
prosite = "C-x(2)-C"
regex = prosite_to_regex(prosite)
print("Padrão PROSITE:", prosite)
print("Regex:", regex)

print("\nBUSCA PROSITE:")
seq_prot = "MCKAACFGCQQC"
posicoes = find_prosite(seq_prot, "C-x(2)-C")
print("Sequência:", seq_prot)
print("Ocorrências (0-based):", posicoes)

print("\nDIGESTÃO POR ENZIMA DE RESTRIÇÃO:")
dna = "GAATTCCGAATTC"
enzima = "G^AATTC"
cuts, frags = digest_dna(dna, enzima)
print("DNA:", dna)
print("Enzima:", enzima)
print("Cortes (0-based):", cuts)
print("Fragmentos:")
for f in frags:
    print(f)

print("\nPWM:")
kmers = ["ACG", "AAG", "ATG", "ACG"]
pwm = build_pwm(kmers)
for i, col in enumerate(pwm):
    print(f"Posição {i}:", col)

print("\nPSSM:")
pssm = pssm_from_pwm(pwm)
for i, col in enumerate(pssm):
    print(f"Posição {i}:", col)

print("\nSCORE DE K-MER:")
kmer = "ACG"
print("K-mer:", kmer)
print("Score:", score_kmer(pssm, kmer))

print("\nMELHOR SUBSEQUÊNCIA (PSSM SCAN):")
seq = "TTACGATGACGTT"
melhor, pos, score = best_subsequence(pssm, seq)
print("Sequência:", seq)
print("Melhor subsequência:", melhor)
print("Posição inicial (0-based):", pos)
print("Score:", score)


## Exemplos de sequencias

from bioinf.sequencias import *

# Verificação de DNA

print("TESTE is_dna:")
is_dna("AGCTAG")
is_dna("AGCTX")

# Contagem de nucleótidos

print("\nTESTE dna_counter:")
print(dna_counter("AGCTAGC"))
print(dna_counter("TTTAAA"))

# Complemento reverso

print("\nTESTE reverse_complement:")
print(reverse_complement("ATGC"))
print(reverse_complement("GGCTA"))

# Transcrição DNA → RNA

print("\nTESTE dna_to_rna:")
print(dna_to_rna("ATGC"))
print(dna_to_rna("GGCTA"))

# Identificação de sequência

print("\nTESTE identificar_sequencia:")
print(identificar_sequencia("ACGT"))
print(identificar_sequencia("ACGU"))
print(identificar_sequencia("MKWVTFISLL"))
print(identificar_sequencia("ACG"))
print(identificar_sequencia("123XYZ"))
