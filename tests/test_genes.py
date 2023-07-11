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
        for g in results:
            for base in g:
                assert base in list(range(20))

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

        float_factory = lambda mult: random()*mult
        gene = genes.Gene.make(5, base_factory=float_factory, factory_args=[-2])
        assert all(-2 <= b <= 0 for b in gene.bases)

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.gene.to_dict()
        assert type(encoded) is dict
        decoded = genes.Gene.from_dict(encoded)
        assert type(decoded) is genes.Gene
        assert self.gene == decoded
        assert self.gene is not decoded

        float_factory = lambda: random()
        gene = genes.Gene.make(5, base_factory=float_factory)
        encoded = gene.to_dict()
        assert type(encoded) is dict
        decoded = genes.Gene.from_dict(encoded)
        assert type(decoded) is genes.Gene
        assert decoded == gene
        assert decoded is not gene

        str_factory = lambda: genes.random_str(5)
        gene = genes.Gene.make(5, base_factory=str_factory)
        encoded = gene.to_dict()
        assert type(encoded) is dict
        decoded = genes.Gene.from_dict(encoded)
        assert type(decoded) is genes.Gene
        assert decoded == gene
        assert decoded is not gene


class TestAllele(unittest.TestCase):
    def first_allele(self) -> genes.Allele:
        return genes.Allele("test", [
            genes.Gene("g1", [1, 2, 3]),
            genes.Gene("g2", [4, 5, 6]),
            genes.Gene("g3", [7, 8, 9]),
        ])

    def alt_allele(self) -> genes.Allele:
        return genes.Allele("test", [
            genes.Gene("g1", [10, 11, 12]),
            genes.Gene("g2", [13, 14, 15]),
            genes.Gene("g3", [16, 17, 18]),
        ])

    def second_allele(self) -> genes.Allele:
        return genes.Allele("diff", [
            genes.Gene("g4", [10, 11, 12]),
            genes.Gene("g5", [13, 14, 15]),
            genes.Gene("g6", [16, 17, 18]),
        ])

    def first_gene(self) -> genes.Gene:
        return genes.Gene("test", [1,2,3])

    def second_gene(self) -> genes.Gene:
        return genes.Gene("diff", [4,5,6])

    def setUp(self) -> None:
        self.allele = self.first_allele()
        self.alternate = self.alt_allele()
        self.other = self.second_allele()
        return super().setUp()

    def test_Allele_is_hashable(self):
        hash(self.allele)

    def test_copy_returns_identical_Allele_containing_identical_genes(self):
        copy = self.allele.copy()
        assert type(copy) is genes.Allele
        assert copy is not self.allele
        assert copy == self.allele

        for a, b in zip(copy.genes, self.allele.genes):
            assert a == b
            assert a.name == b.name
            assert a.bases == b.bases
            assert a is not b

    def test_insert_inserts_specified_or_random_gene_at_specified_or_random_index(self):
        g4 = genes.Gene("g4", [10, 11, 12])
        g5 = genes.Gene("g5", [13, 14, 15])
        self.allele.insert(0, g4)
        assert self.allele.genes[0] == g4
        self.allele.insert(2, g5)
        assert self.allele.genes[2] == g5

        indices = set()
        ref_gene = self.first_gene()
        for _ in range(5):
            allele = self.first_allele()
            allele.insert(gene=ref_gene)
            indices.add(allele.genes.index(ref_gene))
        assert len(indices) > 1

        gs = set()
        for _ in range(5):
            allele = self.first_allele()
            allele.insert(0)
            gs.add(allele.genes[0])
        assert len(gs) >= 4

    def test_duplicate_duplicates_Gene_at_specified_or_random_index(self):
        self.allele.duplicate(0)
        assert self.allele.genes[0] == self.allele.genes[1]
        assert self.allele.genes[0] is not self.allele.genes[1]

        indices = set()
        for _ in range(5):
            allele = self.first_allele()
            allele.duplicate()
            for i in range(len(allele.genes)):
                if allele.genes[i] == allele.genes[i+1]:
                    indices.add(i)
                    break
        assert len(indices) >= 2

    def test_delete_deletes_Gene_at_specified_or_random_index(self):
        assert len(self.allele.genes) == 3
        self.allele.delete(0)
        assert len(self.allele.genes) == 2
        assert self.allele.genes[0].name == "g2"
        assert self.allele.genes[1].name == "g3"

        indices = set()
        for _ in range(5):
            allele = self.first_allele()
            allele.delete()
            if allele.genes[0].name != "g1":
                indices.add(0)
            elif allele.genes[1].name != "g2":
                indices.add(1)
            else:
                indices.add(2)
        assert len(indices) >= 2

    def test_substitute_replaces_Gene_with_specified_or_random_Gene_at_specified_or_random_index(self):
        ref_gene = self.first_gene()
        self.allele.substitute(1, ref_gene)
        assert self.allele.genes[1] == ref_gene

        gs = set()
        indices = set()
        for _ in range(10):
            allele = self.first_allele()
            g_at_1 = allele.genes[1]
            allele.substitute(1)
            assert allele.genes[1] != g_at_1
            gs.add(allele.genes[1])
            allele.substitute(gene=g_at_1)
            indices.add(allele.genes.index(g_at_1))
        assert len(gs) >= 9
        assert len(indices) >= 2

    def test_recombine_swaps_Genes_between_Alleles_at_specified_or_random_indices(self):
        al3 = self.allele.recombine(self.other, [1], False)
        assert al3.genes == [self.allele.genes[0], *self.other.genes[1:]]

        swapset = set()
        for _ in range(5):
            al1 = self.first_allele()
            al2 = self.second_allele()
            al3 = al1.recombine(al2, recombine_genes=False)
            swapped = False
            swaps = []
            for i in range(len(al3.genes)):
                if swapped and al3.genes[i] == al1.genes[i]:
                    swapped = False
                    swaps.append(i)
                elif al3.genes[i] == al2.genes[i]:
                    swapped = True
                    swaps.append(i)
            swapset.add(tuple(swaps))
        assert len(swapset) >= 2

    def test_recombine_recombines_underlying_Genes_with_matching_name_and_index_by_default(self):
        allele = self.allele.recombine(self.alternate, [1])
        recombined_genes = 0
        allele_bases = [b for g in self.allele.genes for b in g.bases]
        other_bases = [b for g in self.other.genes for b in g.bases]
        for g in allele.genes:
            if g not in self.allele.genes and g not in self.alternate.genes:
                recombined_genes += 1
            assert all(b in allele_bases or b in other_bases for b in g.bases)
        assert recombined_genes >= 1

    def test_recombine_does_not_recombine_underlying_Genes_with_different_names_by_default(self):
        allele = self.allele.recombine(self.other, [1])
        recombined_genes = 0
        for g in allele.genes:
            if g not in self.allele.genes and g not in self.other.genes:
                recombined_genes += 1
        assert recombined_genes == 0

    def test_recombine_can_recombine_underlying_Genes_with_different_names(self):
        allele = self.allele.recombine(self.other, [1], match_genes=False)
        recombined_genes = 0
        for g in allele.genes:
            if g not in self.allele.genes and g not in self.other.genes:
                recombined_genes += 1
        assert recombined_genes >= 1

    def test_make_returns_random_Allele_with_random_Genes(self):
        allele = genes.Allele.make(3, 5)
        assert type(allele) is genes.Allele
        assert len(allele.genes) == 3
        assert all(len(g.bases) == 5 for g in allele.genes)

        results = set()
        for _ in range(10):
            results.add(genes.Allele.make(5, 10))
        assert len(results) >= 8
        assert len(set([a.name for a in results])) >= 8
        assert len(set([tuple(a.genes) for a in results])) >= 8

    def test_make_with_base_factory_sets_bases_with_correct_type(self):
        float_factory = lambda: random()
        allele = genes.Allele.make(3, 5, base_factory=float_factory)
        assert all(type(b) is float for g in allele.genes for b in g.bases)
        assert all(len(set(tuple(g.bases))) >= 4 for g in allele.genes)

        float_factory = lambda mult: random()*mult
        allele = genes.Allele.make(3, 5, base_factory=float_factory, factory_args=[-2])
        assert all(-2 <= b <= 0 for g in allele.genes for b in g.bases)

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.allele.to_dict()
        assert type(encoded) is dict
        decoded = genes.Allele.from_dict(encoded)
        assert decoded == self.allele
        assert decoded is not self.allele

        float_factory = lambda: random()
        allele = genes.Allele.make(3, 5, base_factory=float_factory)
        encoded = allele.to_dict()
        assert type(encoded) is dict
        decoded = genes.Allele.from_dict(encoded)
        assert decoded == allele
        assert decoded is not allele

        str_factory = lambda: genes.random_str(5)
        allele = genes.Allele.make(3, 5, base_factory=str_factory)
        encoded = allele.to_dict()
        assert type(encoded) is dict
        decoded = genes.Allele.from_dict(encoded)
        assert decoded == allele
        assert decoded is not allele


if __name__ == '__main__':
    unittest.main()
