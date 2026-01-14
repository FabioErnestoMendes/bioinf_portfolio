import unittest

from bioinf.alinhamento import (
    matriz_de_zeros,
    janela_match,
    dot_plot_janela,
    matriz_substituição_dna,
    primeira_linha_e_coluna,
    melhor_movimento,
    needleman_wunsch,
    smith_waterman,
    reconstruir_alinhamento,
    consenso,
    alinhar_par,
    alinhar_consenso,
    alinhamento_progressivo,
    escolha_de_matriz
)


class TestMatrizesBasicas(unittest.TestCase):

    def test_matriz_de_zeros_dimensoes(self):
        m = matriz_de_zeros("AC", "GT")
        self.assertEqual(len(m), 3)
        self.assertEqual(len(m[0]), 3)

    def test_janela_match_total(self):
        self.assertEqual(janela_match("ACGT", "ACGT", 0, 0, 4), 4)

    def test_janela_match_parcial(self):
        self.assertEqual(janela_match("ACGT", "AGGT", 0, 0, 4), 3)


class TestDotPlot(unittest.TestCase):

    def test_dot_plot_janela_basico(self):
        m = dot_plot_janela("ACGT", "ACGA", 3, 2)
        self.assertEqual(m[0][0], 1)

    def test_dot_plot_sem_matches(self):
        m = dot_plot_janela("AAAA", "TTTT", 2, 2)
        self.assertTrue(all(all(v == 0 for v in linha) for linha in m))


class TestSubstituicao(unittest.TestCase):

    def test_matriz_substituicao_dna(self):
        subst = matriz_substituição_dna(2, -1)
        self.assertEqual(subst["A"]["A"], 2)
        self.assertEqual(subst["A"]["C"], -1)


class TestInicializacao(unittest.TestCase):

    def test_primeira_linha_e_coluna(self):
        matriz, setas = primeira_linha_e_coluna("AC", "GT", -2)
        self.assertEqual(matriz[0][1], -2)
        self.assertEqual(setas[1][0], "↑")


class TestMovimentos(unittest.TestCase):

    def test_melhor_movimento_diagonal(self):
        score, seta = melhor_movimento(5, 3, 1)
        self.assertEqual(seta, "↖")

    def test_melhor_movimento_cima(self):
        score, seta = melhor_movimento(1, 4, 2)
        self.assertEqual(seta, "↑")


class TestAlinhamentoGlobal(unittest.TestCase):

    def test_needleman_sequencias_identicas(self):
        grid, setas = needleman_wunsch("ATGC", "ATGC")
        alin1, alin2 = reconstruir_alinhamento(setas, "ATGC", "ATGC")
        self.assertEqual(alin1, "ATGC")
        self.assertEqual(alin2, "ATGC")

    def test_needleman_tamanho_1(self):
        grid, setas = needleman_wunsch("A", "A")
        alin1, alin2 = reconstruir_alinhamento(setas, "A", "A")
        self.assertEqual(alin1, "A")

    def test_escolha_matriz_invalida(self):
        with self.assertRaises(ValueError):
            escolha_de_matriz("ACG", "XYZ", 2, -1)


class TestAlinhamentoLocal(unittest.TestCase):

    def test_smith_waterman_basico(self):
        grid, setas = smith_waterman("ACGT", "CG")
        alin1, alin2 = reconstruir_alinhamento(setas, "ACGT", "CG")
        self.assertTrue(len(alin1) > 0)
        self.assertEqual(len(alin1), len(alin2))


class TestReconstrucao(unittest.TestCase):

    def test_reconstruir_com_gap(self):
        grid, setas = needleman_wunsch("ACG", "AG")
        a1, a2 = reconstruir_alinhamento(setas, "ACG", "AG")
        self.assertEqual(len(a1), len(a2))


class TestAlinhamentoMultiplo(unittest.TestCase):

    def test_consenso(self):
        self.assertEqual(consenso(["ACG", "A-G", "ACG"]), "ACG")

    def test_alinhar_par(self):
        subst = matriz_substituição_dna(2, -1)
        a1, a2 = alinhar_par("ACG", "AG", subst, -2)
        self.assertEqual(len(a1), len(a2))

    def test_alinhar_consenso(self):
        subst = matriz_substituição_dna(2, -1)
        alinh = alinhar_consenso(["ACG", "A-C"], "AG", subst, -2)
        self.assertEqual(len(alinh), 3)

    def test_alinhamento_progressivo(self):
        subst = matriz_substituição_dna(2, -1)
        seqs = ["ACG", "AC", "AG"]
        alinh = alinhamento_progressivo(seqs, subst)
        self.assertEqual(len(alinh), 3)


if __name__ == "__main__":
    unittest.main()