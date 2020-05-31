import unittest
from cuadripolos import Cuadripolo

a = Cuadripolo()
b = Cuadripolo()
a.cuadri = [[1, 2], [3, 4]]
b.cuadri = [[4, 3], [2, 1]]

class TestSum(unittest.TestCase):
    def test_matrix_mul(self):
        d = a * b
        self.assertEqual(d.cuadri, [[8, 5], [20, 13]])

    def test_matrix_rmul(self):
        d = b * a
        self.assertEqual(d.cuadri, [[13, 20], [5, 8]])

if __name__ == '__main__':
    unittest.main()
