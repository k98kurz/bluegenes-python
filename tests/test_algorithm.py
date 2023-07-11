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

    def test_optimize_gene_exits_after_max_iterations_or_achieving_fitness_target(self):
        target = 123.456

        def measure_fitness(gene: genes.Gene, target: float) -> float:
            return 1 / abs(sum(gene.bases) - target)

        def mutate_gene(gene: genes.Gene) -> genes.Gene:
            for i in range(len(gene.bases)):
                val = random()
                if val <= 0.1:
                    gene.bases[i] *= random()
                elif val <= 0.2:
                    gene.bases[i] /= random()
                elif val <= 0.6:
                    gene.bases[i] += random()
                else:
                    gene.bases[i] -= random()
            return gene

        count, population = algorithm.optimize_gene(
            measure_fitness=measure_fitness,
            mutate_gene=mutate_gene,
            initial_population=self.population(),
            fitness_target=target,
            max_iterations=1000,
            population_size=100
        )

        assert type(count) is int and count <= 1000
        assert type(population) is list
        assert all(type(p) is genes.Gene for p in population)
        assert len(population) == 100
        best = (sum(population[0].bases), population[0])

        assert count < 1000 or (best[0] - target)/target < 0.01


class TestAlgorithmForAlleles(unittest.TestCase):
    def population(self) -> list[genes.Allele]:
        return [
            genes.Allele("test", [
                genes.Gene("gn1", [1.0, 2.0, 3.0, 4.0]),
                genes.Gene("gn2", [2.1, 1.1, 4.1, 3.1]),
            ]),
            genes.Allele("test", [
                genes.Gene("gn1", [1.2, 2.2, 3.2, 4.2]),
                genes.Gene("gn2", [2.3, 1.3, 4.3, 3.3]),
            ]),
            genes.Allele("test", [
                genes.Gene("gn1", [1.4, 2.4, 3.4, 4.4]),
                genes.Gene("gn2", [2.5, 1.5, 4.5, 3.5]),
            ]),
            genes.Allele("test", [
                genes.Gene("gn1", [1.6, 2.6, 3.6, 4.6]),
                genes.Gene("gn2", [2.7, 1.7, 4.7, 3.7]),
            ]),
        ]

    def setUp(self) -> None:
        self.parents = self.population()
        return super().setUp()

    def test_allele_child_from_parents_returns_recombination_of_two_random_parents(self):
        child = algorithm.allele_child_from_parents(self.parents)

        parent_genes = {
            i: self.parents[i].genes
            for i in range(len(self.parents))
        }
        parent_bases = {
            i: [b for g in pg for b in g.bases]
            for i, pg in parent_genes.items()
        }
        parents = set()
        for gene in child.genes:
            for base in gene.bases:
                for k, v in parent_bases.items():
                    if base in v:
                        parents.add(k)
        assert len(parents) <= 2

    def test_optimize_allele_exits_after_max_iterations_or_achieving_fitness_target(self):
        target = 123.456

        def measure_fitness(allele: genes.Allele, target: float) -> float:
            sums = [sum(g.bases) for g in allele.genes]
            return 1 / abs(sum(sums) - target)

        def mutate_gene(gene: genes.Gene) -> genes.Gene:
            for i in range(len(gene.bases)):
                val = random()
                if val <= 0.1:
                    gene.bases[i] *= random()
                elif val <= 0.2:
                    gene.bases[i] /= random()
                elif val <= 0.6:
                    gene.bases[i] += random()
                else:
                    gene.bases[i] -= random()
            return gene

        def mutate_allele(allele: genes.Allele) -> genes.Allele:
            allele.genes = [mutate_gene(g) for g in allele.genes]
            return allele

        count, population = algorithm.optimize_allele(
            measure_fitness=measure_fitness,
            mutate_allele=mutate_allele,
            initial_population=self.population(),
            fitness_target=target,
            max_iterations=1000,
            population_size=100
        )

        assert type(count) is int and count <= 1000
        assert type(population) is list
        assert all(type(p) is genes.Allele for p in population)
        assert len(population) == 100
        sums = [sum(g.bases) for g in population[0].genes]
        best = (sum(sums), population[0])

        assert count < 1000 or (best[0] - target)/target < 0.01


if __name__ == '__main__':
    unittest.main()
