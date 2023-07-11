from context import genes
from random import random
import unittest


class TestGene(unittest.TestCase):
    def setUp(self) -> None:
        self.gene = genes.Gene("test", [1,2,3])
        return super().setUp()

    def test_Gene_is_hashable(self):
        hash(self.gene)

    def test_copy_returns_identical_Gene(self):
        copy = self.gene.copy()
        assert type(copy) is genes.Gene
        assert copy.name == self.gene.name
        assert copy.bases == self.gene.bases
        assert copy.bases is not self.gene.bases

    def test_insert_inserts_base_at_specified_or_random_index(self):
        self.gene.insert(0, 5)
        assert self.gene.bases[0] == 5
        self.gene.insert(3, 7)
        assert self.gene.bases[3] == 7

        indices = set()
        for _ in range(10):
            gene = genes.Gene("test", [1,2,3])
            gene.insert(base=4)
            indices.add(gene.bases.index(4))
        assert len(indices) > 1

    def test_append_adds_base_to_end_of_bases(self):
        self.gene.append(9)
        assert self.gene.bases[-1] == 9
        self.gene.append(10)
        assert self.gene.bases[-1] == 10

    def test_insert_sequence_at_specified_or_random_index(self):
        self.gene.insert_sequence(1, [8,9])
        assert self.gene.bases[1:][:2] == [8,9]

        indices = set()
        for _ in range(10):
            gene = genes.Gene('test', [1,2,3])
            gene.insert_sequence(sequence=[4,5])
            assert 4 in gene.bases
            assert 5 in gene.bases
            assert gene.bases.index(4) + 1 == gene.bases.index(5)
            indices.add(gene.bases.index(4))
        assert len(indices) > 1

    def test_delete_removes_base_at_specified_or_random_index(self):
        self.gene.delete(1)
        assert self.gene.bases == [1,3]

        results = set()
        for _ in range(10):
            gene = genes.Gene("test", [1,2,3])
            gene.delete()
            results.add(tuple(gene.bases))
        assert len(results) > 1

    def test_delete_sequence_removes_bases_at_specified_or_random_index_with_specified_or_random_size(self):
        gene = genes.Gene("big", list(range(10)))
        gene.delete_sequence(3, 5)
        assert gene.bases == [0, 1, 2, 8, 9]

        results = set()
        for _ in range(10):
            gene = genes.Gene("big", list(range(10)))
            gene.delete_sequence()
            results.add(tuple(gene.bases))
        assert len(results) > 1
        assert len([len(r) for r in results]) > 1

    def test_substitute_replaces_base_with_specified_or_random_base_at_specified_or_random_index(self):
        self.gene.substitute(1, 5)
        assert self.gene.bases[1] == 5

        results = set()
        for _ in range(10):
            gene = genes.Gene("big", list(range(10)))
            gene.substitute()
            results.add(tuple(gene.bases))
        assert len(results) > 1
        assert all(len(r) == 10 for r in results)

    def test_recombine_recombines_genes_at_specified_or_random_indices(self):
        gene1 = genes.Gene("base", list(range(5)))
        gene2 = genes.Gene("diff", list(range(5, 10)))
        gene3 = gene1.recombine(gene2, [1, 3])
        assert gene3.bases == [0, 6, 7, 3, 4]
        assert gene3.name != gene1.name
        assert gene3.name[0] == gene1.name[0]
        assert gene3.name != gene2.name
        assert gene3.name[-1] == gene2.name[-1]

        results = set()
        names = set()
        for _ in range(10):
            gene1 = genes.Gene("base", list(range(10)))
            gene2 = genes.Gene("diff", list(range(10, 20)))
            gene3 = gene1.recombine(gene2)
            results.add(tuple(gene3.bases))
            names.add(gene3.name)
        assert len(results) > 1
        assert len(names) > 1

    def test_make_returns_random_Gene(self):
        gene = genes.Gene.make(5)
        assert type(gene) is genes.Gene
        assert len(gene.bases) == 5

        results = set()
        for _ in range(10):
            results.add(genes.Gene.make(10))
        assert len(results) > 8
        assert len(set([g.name for g in results])) >= 9
        assert len(set([tuple(g.bases) for g in results])) >= 9

    def test_make_with_base_factory_sets_bases_with_correct_type(self):
        float_factory = lambda: random()
        gene = genes.Gene.make(5, base_factory=float_factory)
        assert all(type(b) is float for b in gene.bases)
        assert len(set(gene.bases)) >= 4

        float_factory = lambda max: random()*max
        gene = genes.Gene.make(5, base_factory=float_factory, factory_args=[-2])
        assert all(-2 <= b <= 0 for b in gene.bases)

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.gene.to_dict()
        assert type(encoded) is dict
        decoded = genes.Gene.from_dict(encoded)
        assert type(decoded) is genes.Gene
        assert hash(self.gene) == hash(decoded)


if __name__ == '__main__':
    unittest.main()
