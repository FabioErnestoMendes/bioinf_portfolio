import pytest
from bioinf.sequencias import is_dna,dna_counter,dna_to_rna,reverse_complement,identificar_sequencia

def test_is_dna_valida(capsys):
  is_dna("AGCTAG")
  captured = capsys.readouterr()
  assert "Sequência válida" in captured.out

def test_is_dna_invalida(capsys):
  is_dna("AGCTX")
  captured = capsys.readouterr()
  assert "Sequência inválida" in captured.out

def test_is_dna_sequencia_vazia():
  assert is_dna("") is None

def test_dna_counter_normal():
  assert dna_counter("AGCTAGC") == {"A": 2, "G": 2, "C": 2, "T": 1}

def test_dna_counter_tamanho_um():
  assert dna_counter("A") == {"A": 1, "G": 0, "C": 0, "T": 0}

def test_dna_counter_vazio():
  assert dna_counter("") == {"A": 0, "G": 0, "C": 0, "T": 0}

def test_dna_to_rna():
  assert dna_to_rna("ATGC") == "UACG"

def test_reverse_complement():
  assert reverse_complement("ATGC") == "GCAT"

def test_reverse_complement_tamanho_um():
  assert reverse_complement("A") == "T"

def test_identificar_dna():
  assert identificar_sequencia("ACGT") == "DNA"

def test_identificar_rna():
  assert identificar_sequencia("ACGU") == "RNA"

def test_identificar_amino():
  assert identificar_sequencia("MKTLL") == "AMINO"

def test_identificar_insuficiente():
  msg = identificar_sequencia("ACG")
  assert "Introduza mais informação" in msg

def test_identificar_erro():
  assert identificar_sequencia("1234") == "ERRO"