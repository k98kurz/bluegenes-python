from context import genes, algorithm
from random import random
import unittest


class TestAlgorithmForGenes(unittest.TestCase):
    def population(self) -> list[genes.Gene]:
        return [
            genes.Gene("test", [1.0, 2.0, 3.0, 4.0]),
            genes.Gene("test", [2.1, 1.1, 4.1, 3.1]),
            genes.Gene("test", [0.2, 3.2, 2.2, 1.2]),
            genes.Gene("test", [3.3, 4.3, 1.3, 2.3]),
        ]

    def setUp(self) -> None:
        self.parents = self.population()
        return super().setUp()

    def test_gene_child_from_parents_returns_recombination_of_two_random_parents(self):
        child = algorithm.gene_child_from_parents(self.parents)
        print(child)

        parent_bases = {
            i: self.parents[i].bases
            for i in range(len(self.parents))
        }
        parents = set()
        for base in child.bases:
            for k, v in parent_bases.items():
                if base in v:
                    parents.add(k)
        assert len(parents) <= 2


if __name__ == '__main__':
    unittest.main()
