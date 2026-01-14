import unittest
from bioinf.filogenia import distancia, matriz_distancias, upgma

class TestFuncoes(unittest.TestCase):

    def test_distancia(self):
        # Sequências iguais
        self.assertEqual(distancia("ATCG", "ATCG"), 0)

        # Sequências completamente diferentes
        self.assertEqual(distancia("AAA", "TTT"), 3)

        # Casos limite
        self.assertEqual(distancia("", ""), 0)
        self.assertEqual(distancia("A", ""), 1)
        self.assertEqual(distancia("", "A"), 1)

    def test_distancia_typeerror(self):
        with self.assertRaises(TypeError):
            distancia(None, "A")
        with self.assertRaises(TypeError):
            distancia("A", 123)

    def test_matriz_distancias(self):
        sequencias = ["AAA", "AAT", "ATT", "TTT"]
        matriz = matriz_distancias(sequencias)
        n = len(matriz)

        # Verifica se a matriz é simétrica
        for i in range(n):
            for j in range(n):
                self.assertEqual(matriz[i][j], matriz[j][i])

        # Verifica se a diagonal principal é zero
        for i in range(n):
            self.assertEqual(matriz[i][i], 0)

        # Valores conhecidos
        self.assertEqual(matriz[0][1], 1)
        self.assertEqual(matriz[0][3], 3)

    def test_matriz_distancias_typeerror(self):
        with self.assertRaises(TypeError):
            matriz_distancias("AAA")  # não é lista
        with self.assertRaises(TypeError):
            matriz_distancias(["AAA", 1])  # elemento não-string

    def test_upgma(self):
        sequencias = ["CCG", "GT", "GTA", "AAT", "AT", "ACG", "ACGT"]
        D = matriz_distancias(sequencias)

        ligacoes = upgma(D, sequencias.copy())

        self.assertEqual(len(ligacoes), len(sequencias) - 1)

        for lig in ligacoes:
            self.assertEqual(len(lig), 3)
            self.assertIsInstance(lig[0], str)
            self.assertIsInstance(lig[1], str)
            self.assertIsInstance(lig[2], (int, float))

    def test_upgma_muta_distancias_e_labels(self):
        labels = ["A", "B", "C"]
        dist = [
            [0, 2, 4],
            [2, 0, 3],
            [4, 3, 0],
        ]
        upgma(dist, labels)
        self.assertEqual(len(labels), 1)
        self.assertEqual(len(dist), 1)
        self.assertEqual(len(dist[0]), 1)

    def test_upgma_erros_labels(self):
        with self.assertRaises(TypeError):
            upgma([[0, 1], [1, 0]], "AB")
        with self.assertRaises(TypeError):
            upgma([[0, 1], [1, 0]], ["A", 1])
        with self.assertRaises(ValueError):
            upgma([[0]], ["A"])  # N < 2

    def test_upgma_erros_matriz(self):
        with self.assertRaises(TypeError):
            upgma("matriz", ["A", "B"])
        with self.assertRaises(TypeError):
            upgma([0, 1], ["A", "B"])
        with self.assertRaises(ValueError):
            upgma([[0, 1, 2], [1, 0, 3]], ["A", "B", "C"])
        with self.assertRaises(TypeError):
            upgma([[0, "x"], [1, 0]], ["A", "B"])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFuncoes)
    unittest.TextTestRunner(verbosity=2).run(suite)