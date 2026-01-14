import unittest
from bioinf.blast import novo_indice, busca_pares, estender_alem, alinhamento_pro


import unittest
from bioinf.blast import novo_indice, busca_pares, estender_alem, alinhamento_pro

class TestBlast(unittest.TestCase):

    def test_novo_indice_k_invalido(self):
        with self.assertRaises(ValueError):
            novo_indice("ABCDE", k=0)
        with self.assertRaises(ValueError):
            novo_indice("ABCDE", k=-1)

    def test_novo_indice_normal(self):
        idx = novo_indice("ATATAT", k=3)
        self.assertEqual(idx["ATA"], [0, 2])
        self.assertEqual(idx["TAT"], [1, 3])

    def test_novo_indice_sequencia_menor_que_k(self):
        idx = novo_indice("AB", k=3)
        self.assertEqual(len(idx), 0)

    def test_busca_pares_normal(self):
        pares = busca_pares("ATATAT", "TAT", k=3)
        self.assertEqual(pares, [(1, 0), (3, 0)])

    def test_busca_pares_sem_hits(self):
        self.assertEqual(busca_pares("AAAA", "TTTT", k=2), [])

    def test_busca_pares_subject_menor_que_k(self):
        self.assertEqual(busca_pares("ABCDE", "AB", k=3), [])

    def test_estender_alem_extende_esquerda_e_direita(self):
        start_i, start_j, tamanho = estender_alem("ABCDE", "ABCDE", i=1, j=1, k=3)
        self.assertEqual((start_i, start_j, tamanho), (0, 0, 5))

    def test_estender_alem_sem_extensao(self):
        start_i, start_j, tamanho = estender_alem("ABCXX", "ABCYY", i=0, j=0, k=3)
        self.assertEqual((start_i, start_j, tamanho), (0, 0, 3))

    def test_alinhamento_pro_sem_seeds(self):
        self.assertIsNone(alinhamento_pro("AB", "ABCDE", k=3))
        self.assertIsNone(alinhamento_pro("ABCDE", "AB", k=3))

    def test_alinhamento_pro_sem_hits(self):
        self.assertIsNone(alinhamento_pro("AAAA", "TTTT", k=2))

    def test_alinhamento_pro_melhor_alinhamento(self):
        query = "AACCTT"
        subject = "GGACCTTA"
        res = alinhamento_pro(query, subject, k=3)

        self.assertIsNotNone(res)
        self.assertEqual(res["alinhado_q"], "ACCTT")
        self.assertEqual(res["alinhado_s"], "ACCTT")
        self.assertEqual(res["tamanho"], 5)
        self.assertEqual(res["query_start"], 1)
        self.assertEqual(res["subject_start"], 2)


if __name__ == "__main__":
    unittest.main()
