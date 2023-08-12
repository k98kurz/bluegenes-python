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


class TestNucleosome(unittest.TestCase):
    def first_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("test", [
            genes.Gene("g1", [1, 2, 3]),
            genes.Gene("g2", [4, 5, 6]),
            genes.Gene("g3", [7, 8, 9]),
            genes.Gene("g4", [-7, -8, -9]),
        ])

    def alt_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("test", [
            genes.Gene("g1", [10, 11, 12]),
            genes.Gene("g2", [13, 14, 15]),
            genes.Gene("g3", [16, 17, 18]),
            genes.Gene("g4", [-16, -17, -18]),
        ])

    def second_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("diff", [
            genes.Gene("g4", [10, 11, 12]),
            genes.Gene("g5", [13, 14, 15]),
            genes.Gene("g6", [16, 17, 18]),
            genes.Gene("g6", [116, 117, 118]),
        ])

    def first_gene(self) -> genes.Gene:
        return genes.Gene("test", [1,2,3])

    def second_gene(self) -> genes.Gene:
        return genes.Gene("diff", [4,5,6])

    def setUp(self) -> None:
        self.nucleosome = self.first_nucleosome()
        self.alternate = self.alt_nucleosome()
        self.other = self.second_nucleosome()
        return super().setUp()

    def test_Nucleosome_is_hashable(self):
        hash(self.nucleosome)

    def test_copy_returns_identical_Nucleosome_containing_identical_genes(self):
        copy = self.nucleosome.copy()
        assert type(copy) is genes.Nucleosome
        assert copy is not self.nucleosome
        assert copy == self.nucleosome

        for a, b in zip(copy.genes, self.nucleosome.genes):
            assert a == b
            assert a.name == b.name
            assert a.bases == b.bases
            assert a is not b

    def test_insert_inserts_specified_or_random_gene_at_specified_or_random_index(self):
        g4 = genes.Gene("g4", [10, 11, 12])
        g5 = genes.Gene("g5", [13, 14, 15])
        self.nucleosome.insert(0, g4)
        assert self.nucleosome.genes[0] == g4
        self.nucleosome.insert(2, g5)
        assert self.nucleosome.genes[2] == g5

        indices = set()
        ref_gene = self.first_gene()
        for _ in range(10):
            nucleosome = self.first_nucleosome()
            nucleosome.insert(gene=ref_gene)
            indices.add(nucleosome.genes.index(ref_gene))
        assert len(indices) > 1

        gs = set()
        for _ in range(10):
            nucleosome = self.first_nucleosome()
            nucleosome.insert(0)
            gs.add(nucleosome.genes[0])
        assert len(gs) >= 4

    def test_append_adds_Gene_to_end_of_genes(self):
        gg1 = genes.Gene("gg1", [1,2,3])
        gg2 = genes.Gene("gg2", [1,2,3])
        self.nucleosome.append(gg1)
        assert self.nucleosome.genes[-1] == gg1
        self.nucleosome.append(gg2)
        assert self.nucleosome.genes[-1] == gg2

    def test_duplicate_duplicates_Gene_at_specified_or_random_index(self):
        self.nucleosome.duplicate(0)
        assert self.nucleosome.genes[0] == self.nucleosome.genes[1]
        assert self.nucleosome.genes[0] is not self.nucleosome.genes[1]

        indices = set()
        for _ in range(10):
            nucleosome = self.first_nucleosome()
            nucleosome.duplicate()
            for i in range(len(nucleosome.genes)):
                if nucleosome.genes[i] == nucleosome.genes[i+1]:
                    indices.add(i)
                    break
        assert len(indices) >= 2

    def test_delete_deletes_Gene_at_specified_or_random_index(self):
        assert len(self.nucleosome.genes) == 4
        self.nucleosome.delete(0)
        assert len(self.nucleosome.genes) == 3
        assert self.nucleosome.genes[0].name == "g2"
        assert self.nucleosome.genes[1].name == "g3"
        assert self.nucleosome.genes[2].name == "g4"

        indices = set()
        for _ in range(10):
            nucleosome = self.first_nucleosome()
            nucleosome.delete()
            if nucleosome.genes[0].name != "g1":
                indices.add(0)
            elif nucleosome.genes[1].name != "g2":
                indices.add(1)
            elif nucleosome.genes[2].name != "g3":
                indices.add(2)
            else:
                indices.add(3)
        assert len(indices) >= 2

    def test_substitute_replaces_Gene_with_specified_or_random_Gene_at_specified_or_random_index(self):
        ref_gene = self.first_gene()
        self.nucleosome.substitute(1, ref_gene)
        assert self.nucleosome.genes[1] == ref_gene

        gs = set()
        indices = set()
        for _ in range(10):
            nucleosome = self.first_nucleosome()
            g_at_1 = nucleosome.genes[1]
            nucleosome.substitute(1)
            assert nucleosome.genes[1] != g_at_1
            gs.add(nucleosome.genes[1])
            nucleosome.substitute(gene=g_at_1)
            indices.add(nucleosome.genes.index(g_at_1))
        assert len(gs) >= 9
        assert len(indices) >= 2

    def test_recombine_swaps_Genes_between_Nucleosomes_at_specified_or_random_indices(self):
        al3 = self.nucleosome.recombine(self.other, [1], False)
        assert al3.genes == [self.nucleosome.genes[0], *self.other.genes[1:]]

        swapset = set()
        for _ in range(10):
            al1 = self.first_nucleosome()
            al2 = self.second_nucleosome()
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

    def test_recombine_recombines_underlying_Genes_with_different_name_and_index_by_default(self):
        nucleosome = self.nucleosome.recombine(self.other, [1])
        recombined_genes = 0
        nucleosome_bases = [b for g in self.nucleosome.genes for b in g.bases]
        other_bases = [b for g in self.other.genes for b in g.bases]
        for g in nucleosome.genes:
            if g not in self.nucleosome.genes and g not in self.alternate.genes:
                recombined_genes += 1
            assert all(b in nucleosome_bases or b in other_bases for b in g.bases)
        assert recombined_genes >= 1

    def test_recombine_does_not_recombine_underlying_Genes_with_different_names_when_specified(self):
        nucleosome = self.nucleosome.recombine(self.other, [1], match_genes=True)
        recombined_genes = 0
        for g in nucleosome.genes:
            if g not in self.nucleosome.genes and g not in self.other.genes:
                recombined_genes += 1
        assert recombined_genes == 0

    def test_recombine_can_recombine_underlying_Genes_with_different_names(self):
        nucleosome = self.nucleosome.recombine(self.other, [1], match_genes=False)
        recombined_genes = 0
        for g in nucleosome.genes:
            if g not in self.nucleosome.genes and g not in self.other.genes:
                recombined_genes += 1
        assert recombined_genes >= 1

    def test_make_returns_random_Nucleosome_with_random_Genes(self):
        nucleosome = genes.Nucleosome.make(3, 5)
        assert type(nucleosome) is genes.Nucleosome
        assert len(nucleosome.genes) == 3
        assert all(len(g.bases) == 5 for g in nucleosome.genes)

        results = set()
        for _ in range(10):
            results.add(genes.Nucleosome.make(5, 10))
        assert len(results) >= 8
        assert len(set([a.name for a in results])) >= 8
        assert len(set([tuple(a.genes) for a in results])) >= 8

    def test_make_with_base_factory_sets_bases_with_correct_type(self):
        float_factory = lambda: random()
        nucleosome = genes.Nucleosome.make(3, 5, base_factory=float_factory)
        assert all(type(b) is float for g in nucleosome.genes for b in g.bases)
        assert all(len(set(tuple(g.bases))) >= 4 for g in nucleosome.genes)

        float_factory = lambda mult: random()*mult
        nucleosome = genes.Nucleosome.make(3, 5, base_factory=float_factory, factory_args=[-2])
        assert all(-2 <= b <= 0 for g in nucleosome.genes for b in g.bases)

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.nucleosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Nucleosome.from_dict(encoded)
        assert decoded == self.nucleosome
        assert decoded is not self.nucleosome

        float_factory = lambda: random()
        nucleosome = genes.Nucleosome.make(3, 5, base_factory=float_factory)
        encoded = nucleosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Nucleosome.from_dict(encoded)
        assert decoded == nucleosome
        assert decoded is not nucleosome

        str_factory = lambda: genes.random_str(5)
        nucleosome = genes.Nucleosome.make(3, 5, base_factory=str_factory)
        encoded = nucleosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Nucleosome.from_dict(encoded)
        assert decoded == nucleosome
        assert decoded is not nucleosome


class TestChromosome(unittest.TestCase):
    def first_gene(self) -> genes.Gene:
        return genes.Gene("test", [1,2,3])

    def second_gene(self) -> genes.Gene:
        return genes.Gene("diff", [4,5,6])

    def first_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("test", [
            genes.Gene("g1", [1, 2, 3]),
            genes.Gene("g2", [4, 5, 6]),
            genes.Gene("g3", [7, 8, 9]),
            genes.Gene("g4", [-7, -8, -9]),
        ])

    def alt_nucleosome(self, name: str = "test") -> genes.Nucleosome:
        return genes.Nucleosome(name, [
            genes.Gene("g1", [10, 11, 12]),
            genes.Gene("g2", [13, 14, 15]),
            genes.Gene("g3", [16, 17, 18]),
            genes.Gene("g4", [-16, -17, -18]),
        ])

    def second_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("diff", [
            genes.Gene("g4", [10, 11, 12]),
            genes.Gene("g5", [13, 14, 15]),
            genes.Gene("g6", [16, 17, 18]),
            genes.Gene("g6", [116, 117, 118]),
        ])

    def third_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("fizz", [
            genes.Gene("bz1", [20, 21, 22]),
            genes.Gene("bz2", [30, 31, 32]),
            genes.Gene("bz3", [40, 41, 42]),
            genes.Gene("bz4", [140, 141, 142]),
        ])

    def first_chromosome(self) -> genes.Chromosome:
        return genes.Chromosome("test", [
            self.first_nucleosome(),
            self.second_nucleosome()
        ])

    def alt_chromosome(self) -> genes.Chromosome:
        return genes.Chromosome("test", [
            self.alt_nucleosome(),
            self.second_nucleosome()
        ])

    def second_chromosome(self) -> genes.Chromosome:
        return genes.Chromosome("diff", [
            self.second_nucleosome(),
            self.third_nucleosome(),
        ])

    def setUp(self) -> None:
        self.chromosome = self.first_chromosome()
        self.alternate = self.alt_chromosome()
        self.other = self.second_chromosome()
        return super().setUp()

    def test_Chromosome_is_hashable(self):
        hash(self.chromosome)

    def test_copy_returns_identical_Chromosome_containing_identical_Nucleosomes(self):
        copy = self.chromosome.copy()
        assert type(copy) is genes.Chromosome
        assert copy is not self.chromosome
        assert copy == self.chromosome

        for a, b in zip(copy.nucleosomes, self.chromosome.nucleosomes):
            assert a == b
            assert a.name == b.name
            assert a.genes == b.genes
            assert a is not b

    def test_insert_inserts_specified_or_random_nucleosome_at_specified_or_random_index(self):
        al4 = genes.Nucleosome.make(3, 5, "al4")
        al5 = genes.Nucleosome.make(3, 5, "al5")
        self.chromosome.insert(0, al4)
        assert self.chromosome.nucleosomes[0] == al4
        self.chromosome.insert(2, al5)
        assert self.chromosome.nucleosomes[2] == al5

        indices = set()
        ref_nucleosome = self.first_nucleosome().append(self.second_gene())
        for _ in range(10):
            chromosome = self.alt_chromosome()
            chromosome.insert(nucleosome=ref_nucleosome)
            indices.add(chromosome.nucleosomes.index(ref_nucleosome))
        assert len(indices) > 1

        gs = set()
        for _ in range(10):
            chromosome = self.first_chromosome()
            chromosome.insert(0)
            gs.add(chromosome.nucleosomes[0])
        assert len(gs) >= 4

    def test_append_adds_Nucleosome_to_end_of_nucleosomes(self):
        gl1 = genes.Nucleosome("gl1", [self.first_gene()])
        gl2 = genes.Nucleosome("gl2", [self.second_gene()])
        self.chromosome.append(gl1)
        assert self.chromosome.nucleosomes[-1] == gl1
        self.chromosome.append(gl2)
        assert self.chromosome.nucleosomes[-1] == gl2

    def test_duplicate_duplicates_Nucleosome_at_specified_or_random_index(self):
        self.chromosome.duplicate(0)
        assert self.chromosome.nucleosomes[0] == self.chromosome.nucleosomes[1]
        assert self.chromosome.nucleosomes[0] is not self.chromosome.nucleosomes[1]

        indices = set()
        for _ in range(10):
            chromosome = self.first_chromosome().append(self.alt_nucleosome())
            chromosome.duplicate()
            for i in range(len(chromosome.nucleosomes)):
                if chromosome.nucleosomes[i] == chromosome.nucleosomes[i+1]:
                    indices.add(i)
                    break
        assert len(indices) >= 2

    def test_delete_deletes_Nucleosome_at_specified_or_random_index(self):
        assert len(self.chromosome.nucleosomes) == 2
        self.chromosome.delete(0)
        assert len(self.chromosome.nucleosomes) == 1
        assert self.chromosome.nucleosomes[0].name == "diff"

        indices = set()
        for _ in range(10):
            chromosome = self.first_chromosome().append(self.alt_nucleosome("fizz"))
            chromosome.delete()
            if chromosome.nucleosomes[0].name != "test":
                indices.add(0)
            elif chromosome.nucleosomes[1].name != "diff":
                indices.add(1)
            else:
                indices.add(2)
        assert len(indices) >= 2

    def test_substitute_replaces_Nucleosome_with_specified_or_random_Nucleosome_at_specified_or_random_index(self):
        ref_nucleosome = self.first_nucleosome()
        assert self.chromosome.nucleosomes[1] != ref_nucleosome
        self.chromosome.substitute(1, ref_nucleosome)
        assert self.chromosome.nucleosomes[1] == ref_nucleosome

        gs = set()
        indices = set()
        for _ in range(20):
            chromosome = self.first_chromosome()
            a_at_1 = chromosome.nucleosomes[1]
            chromosome.substitute(1)
            assert chromosome.nucleosomes[1] != a_at_1
            gs.add(chromosome.nucleosomes[1])
            chromosome.substitute(nucleosome=a_at_1)
            indices.add(chromosome.nucleosomes.index(a_at_1))
        assert len(gs) >= 9
        assert len(indices) >= 2

    def test_recombine_swaps_Nucleosomes_between_Chromosomes_at_specified_or_random_indices(self):
        al3 = self.chromosome.recombine(self.other, [1], False)
        assert al3.nucleosomes == [self.chromosome.nucleosomes[0], *self.other.nucleosomes[1:]]

        swapset = set()
        for _ in range(10):
            ch1 = self.first_chromosome()
            ch2 = self.second_chromosome()
            ch3 = ch1.recombine(ch2, match_nucleosomes=True, recombine_genes=False)
            swapped = False
            swaps = []
            for i in range(len(ch3.nucleosomes)):
                if swapped and ch3.nucleosomes[i] == ch1.nucleosomes[i]:
                    swapped = False
                    swaps.append(i)
                elif ch3.nucleosomes[i] == ch2.nucleosomes[i]:
                    swapped = True
                    swaps.append(i)
            swapset.add(tuple(swaps))
        assert len(swapset) >= 2

    def test_recombine_recombines_underlying_Nucleosomes_and_Genes_with_matching_name_and_index_by_default(self):
        chromosome = self.chromosome.recombine(self.alternate, [1])
        main_genes = [g for a in self.chromosome.nucleosomes for g in a.genes]
        main_gene_names = [g.name for g in main_genes]
        main_bases = [b for g in main_genes for b in g.bases]
        other_genes = [b for g in self.alternate.nucleosomes for b in g.genes]
        other_gene_names = [g.name for g in other_genes]
        other_bases = [b for g in other_genes for b in g.bases]
        recombined_nucleosomes = 0
        recombined_genes = 0
        for c in chromosome.nucleosomes:
            if c not in self.chromosome.nucleosomes and c not in self.alternate.nucleosomes:
                recombined_nucleosomes += 1
            for g in c.genes:
                if g not in main_genes and g not in other_genes:
                    recombined_genes += 1
                assert g.name in main_gene_names or g.name in other_gene_names
                assert all(b in main_bases or b in other_bases for b in g.bases)
        assert recombined_nucleosomes >= 1
        assert recombined_genes >= 1

    def test_recombine_does_not_recombine_underlying_Nucleosomes_or_Genes_with_different_names_when_specified(self):
        chromosome = self.chromosome.recombine(self.other, [1], match_nucleosomes=True)
        main_genes = [g for a in self.chromosome.nucleosomes for g in a.genes]
        other_genes = [b for g in self.other.nucleosomes for b in g.genes]
        recombined_nucleosomes = 0
        recombined_genes = 0

        for a in chromosome.nucleosomes:
            if a not in self.chromosome.nucleosomes and a not in self.other.nucleosomes:
                recombined_nucleosomes += 1
            for g in a.genes:
                if g not in main_genes and g not in other_genes:
                    recombined_genes += 1
        assert recombined_nucleosomes == 0
        assert recombined_genes == 0

    def test_recombine_can_recombine_underlying_Nucleosomes_and_Genes_with_different_names(self):
        chromosome = self.chromosome.recombine(
            self.other, [1], match_nucleosomes=False, match_genes=False
        )
        main_genes = [g for a in self.chromosome.nucleosomes for g in a.genes]
        main_gene_names = [g.name for g in main_genes]
        other_genes = [b for g in self.other.nucleosomes for b in g.genes]
        other_gene_names = [g.name for g in other_genes]
        recombined_nucleosomes = 0
        recombined_genes = 0
        for a in chromosome.nucleosomes:
            if a not in self.chromosome.nucleosomes and a not in self.other.nucleosomes:
                recombined_nucleosomes += 1
            for g in a.genes:
                if g.name not in main_gene_names and g.name not in other_gene_names:
                    recombined_genes += 1
        assert recombined_nucleosomes >= 1
        assert recombined_genes >= 1

    def test_make_returns_random_Chromosome_with_random_Nucleosomes_with_random_Genes(self):
        chromosome = genes.Chromosome.make(3, 5, 7)
        assert type(chromosome) is genes.Chromosome
        assert len(chromosome.nucleosomes) == 3
        assert all(len(a.genes) == 5 for a in chromosome.nucleosomes)
        assert all(len(g.bases) == 7 for a in chromosome.nucleosomes for g in a.genes)

        results: set[genes.Chromosome] = set()
        for _ in range(10):
            results.add(genes.Chromosome.make(5, 10, 20))
        assert len(results) >= 8
        assert len(set([c.name for c in results])) >= 8
        assert len(set([tuple(c.nucleosomes) for c in results])) >= 8
        assert len(set([tuple(a.genes) for c in results for a in c.nucleosomes])) >= 8

    def test_make_with_base_factory_sets_bases_with_correct_type(self):
        float_factory = lambda: random()
        chromosome = genes.Chromosome.make(3, 5, 7, base_factory=float_factory)
        assert all(
            type(b) is float
            for a in chromosome.nucleosomes
            for g in a.genes
            for b in g.bases
        )
        assert all(
            len(set(tuple(g.bases))) >= 6
            for a in chromosome.nucleosomes
            for g in a.genes
        )

        float_factory = lambda mult: random()*mult
        chromosome = genes.Chromosome.make(3, 5, 7, base_factory=float_factory, factory_args=[-2])
        assert all(
            -2 <= b <= 0
            for a in chromosome.nucleosomes
            for g in a.genes
            for b in g.bases
        )

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.chromosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Chromosome.from_dict(encoded)
        assert decoded == self.chromosome
        assert decoded is not self.chromosome

        float_factory = lambda: random()
        chromosome = genes.Chromosome.make(3, 5, 7, base_factory=float_factory)
        encoded = chromosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Chromosome.from_dict(encoded)
        assert decoded == chromosome
        assert decoded is not chromosome

        str_factory = lambda: genes.random_str(5)
        chromosome = genes.Chromosome.make(3, 5, 7, base_factory=str_factory)
        encoded = chromosome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Chromosome.from_dict(encoded)
        assert decoded == chromosome
        assert decoded is not chromosome


class TestGnome(unittest.TestCase):
    def first_gene(self) -> genes.Gene:
        return genes.Gene("test", [1,2,3])

    def second_gene(self) -> genes.Gene:
        return genes.Gene("diff", [4,5,6])

    def first_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("test", [
            genes.Gene("g1", [1, 2, 3]),
            genes.Gene("g2", [4, 5, 6]),
            genes.Gene("g3", [7, 8, 9]),
            genes.Gene("g4", [-7, -8, -9]),
        ])

    def alt_nucleosome(self, name: str = "test") -> genes.Nucleosome:
        return genes.Nucleosome(name, [
            genes.Gene("g1", [10, 11, 12]),
            genes.Gene("g2", [13, 14, 15]),
            genes.Gene("g3", [16, 17, 18]),
            genes.Gene("g4", [-16, -17, -18]),
        ])

    def second_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("diff", [
            genes.Gene("g4", [10, 11, 12]),
            genes.Gene("g5", [13, 14, 15]),
            genes.Gene("g6", [16, 17, 18]),
            genes.Gene("g6", [116, 117, 118]),
        ])

    def third_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("fizz", [
            genes.Gene("bz1", [20, 21, 22]),
            genes.Gene("bz2", [30, 31, 32]),
            genes.Gene("bz3", [40, 41, 42]),
            genes.Gene("bz4", [140, 141, 142]),
        ])

    def fourth_nucleosome(self) -> genes.Nucleosome:
        return genes.Nucleosome("buzz", [
            genes.Gene("bz4", [220, 221, 222]),
            genes.Gene("bz5", [230, 231, 232]),
            genes.Gene("bz6", [240, 241, 242]),
            genes.Gene("bz7", [250, 251, 252]),
        ])

    def first_chromosome(self) -> genes.Chromosome:
        return genes.Chromosome("test", [
            self.first_nucleosome(),
            self.second_nucleosome()
        ])

    def alt_chromosome(self, name: str = "test") -> genes.Chromosome:
        return genes.Chromosome(name, [
            self.alt_nucleosome(),
            self.second_nucleosome()
        ])

    def second_chromosome(self, name: str = "diff") -> genes.Chromosome:
        return genes.Chromosome(name, [
            self.second_nucleosome(),
            self.third_nucleosome(),
        ])

    def third_chromosome(self) -> genes.Chromosome:
        return genes.Chromosome("fizz", [
            self.third_nucleosome(),
            self.fourth_nucleosome(),
        ])

    def first_genome(self) -> genes.Genome:
        return genes.Genome("test", [
            self.first_chromosome(),
            self.second_chromosome(),
        ])

    def alt_genome(self) -> genes.Genome:
        return genes.Genome("test", [
            self.alt_chromosome(),
            self.second_chromosome(),
        ])

    def second_genome(self) -> genes.Genome:
        return genes.Genome("diff", [
            self.second_chromosome("asdf"),
            self.third_chromosome(),
        ])

    def setUp(self) -> None:
        self.genome = self.first_genome()
        self.alternate = self.alt_genome()
        self.other = self.second_genome()
        return super().setUp()

    def test_Genome_is_hashable(self):
        hash(self.genome)

    def test_copy_returns_identical_Genome_containing_identical_Chromosomes(self):
        copy = self.genome.copy()
        assert type(copy) is genes.Genome
        assert copy is not self.genome
        assert copy == self.genome

        for a, b in zip(copy.chromosomes, self.genome.chromosomes):
            assert a == b
            assert a.name == b.name
            assert a.nucleosomes == b.nucleosomes
            assert a is not b

    def test_insert_inserts_specified_or_random_Chromosome_at_specified_or_random_index(self):
        ch4 = genes.Chromosome.make(3, 5, 7, "ch4")
        ch5 = genes.Chromosome.make(3, 5, 7, "ch5")
        self.genome.insert(0, ch4)
        assert self.genome.chromosomes[0] == ch4
        self.genome.insert(2, ch5)
        assert self.genome.chromosomes[2] == ch5

        indices = set()
        ref_chromosome = self.first_chromosome()
        for _ in range(10):
            genome = self.alt_genome().append(self.third_chromosome())
            genome.insert(chromosome=ref_chromosome)
            indices.add(genome.chromosomes.index(ref_chromosome))
        assert len(indices) > 1

        gs = set()
        for _ in range(10):
            genome = self.first_genome()
            genome.insert(0)
            gs.add(genome.chromosomes[0])
        assert len(gs) >= 4

    def test_append_adds_Chromosome_to_end_of_chromosomes(self):
        ch1 = genes.Chromosome("ch1", [self.first_chromosome()])
        ch2 = genes.Chromosome("ch2", [self.second_chromosome()])
        self.genome.append(ch1)
        assert self.genome.chromosomes[-1] == ch1
        self.genome.append(ch2)
        assert self.genome.chromosomes[-1] == ch2

    def test_duplicate_duplicates_Chromosome_at_specified_or_random_index(self):
        self.genome.duplicate(0)
        assert self.genome.chromosomes[0] == self.genome.chromosomes[1]
        assert self.genome.chromosomes[0] is not self.genome.chromosomes[1]

        indices = set()
        for _ in range(10):
            genome = self.first_genome().append(self.alt_chromosome())
            genome.duplicate()
            for i in range(len(genome.chromosomes)):
                if genome.chromosomes[i] == genome.chromosomes[i+1]:
                    indices.add(i)
                    break
        assert len(indices) >= 2

    def test_delete_deletes_Chromosome_at_specified_or_random_index(self):
        assert len(self.genome.chromosomes) == 2
        self.genome.delete(0)
        assert len(self.genome.chromosomes) == 1
        assert self.genome.chromosomes[0].name == "diff"

        indices = set()
        for _ in range(10):
            genome = self.first_genome().append(self.alt_chromosome("fizz"))
            genome.delete()
            if genome.chromosomes[0].name != "test":
                indices.add(0)
            elif genome.chromosomes[1].name != "diff":
                indices.add(1)
            else:
                indices.add(2)
        assert len(indices) >= 2

    def test_substitute_replaces_Chromosome_with_specified_or_random_Chromosome_at_specified_or_random_index(self):
        ref_chromosome = self.first_chromosome()
        assert self.genome.chromosomes[1] != ref_chromosome
        self.genome.substitute(1, ref_chromosome)
        assert self.genome.chromosomes[1] == ref_chromosome

        gs = set()
        indices = set()
        for _ in range(10):
            genome = self.first_genome().append(self.third_chromosome())
            c_at_1 = genome.chromosomes[1]
            genome.substitute(1)
            assert genome.chromosomes[1] != c_at_1
            gs.add(genome.chromosomes[1])
            genome.substitute(chromosome=c_at_1)
            indices.add(genome.chromosomes.index(c_at_1))
        assert len(gs) >= 9
        assert len(indices) >= 2

    def test_recombine_swaps_Chromosomes_between_Genomes_at_specified_or_random_indices(self):
        gm3 = self.genome.recombine(self.other, [1], False)
        assert gm3.chromosomes == [self.genome.chromosomes[0], *self.other.chromosomes[1:]]

        swapset = set()
        for _ in range(10):
            gm1 = self.first_genome()
            gm2 = self.second_genome()
            gm3 = gm1.recombine(gm2, recombine_chromosomes=False)
            swapped = False
            swaps = []
            for i in range(len(gm3.chromosomes)):
                if swapped and gm3.chromosomes[i] == gm1.chromosomes[i]:
                    swapped = False
                    swaps.append(i)
                elif gm3.chromosomes[i] == gm2.chromosomes[i]:
                    swapped = True
                    swaps.append(i)
            swapset.add(tuple(swaps))
        assert len(swapset) >= 2

    def test_recombine_recombines_underlying_Chromosomes_Nucleosomes_and_Genes_with_matching_name_and_index_by_default(self):
        genome = self.genome.recombine(self.alternate, [1])
        main_nucleosomes = [a for c in self.genome.chromosomes for a in c.nucleosomes]
        main_genes = [g for a in main_nucleosomes for g in a.genes]
        main_gene_names = [g.name for g in main_genes]
        main_bases = [b for g in main_genes for b in g.bases]
        other_nucleosomes = [a for c in self.alternate.chromosomes for a in c.nucleosomes]
        other_genes = [g for a in other_nucleosomes for g in a.genes]
        other_gene_names = [g.name for g in other_genes]
        other_bases = [b for g in other_genes for b in g.bases]
        recombined_chromosomes = 0
        recombined_nucleosomes = 0
        recombined_genes = 0
        for c in genome.chromosomes:
            if c not in self.genome.chromosomes and c not in self.alternate.chromosomes:
                recombined_chromosomes += 1
            for a in c.nucleosomes:
                if a not in main_nucleosomes and a not in other_nucleosomes:
                    recombined_nucleosomes += 1
                for g in a.genes:
                    if g not in main_genes and g not in other_genes:
                        recombined_genes += 1
                    assert g.name in main_gene_names or g.name in other_gene_names
                    assert all(b in main_bases or b in other_bases for b in g.bases)
        assert recombined_chromosomes >= 1
        assert recombined_nucleosomes >= 1
        assert recombined_genes >= 1

    def test_recombine_does_not_recombine_underlying_Chromosomes_or_Genes_with_different_names_when_specified(self):
        genome = self.genome.recombine(self.other, [1], match_chromosomes=True, match_nucleosomes=True, match_genes=True)
        main_nucleosomes = [a for c in self.genome.chromosomes for a in c.nucleosomes]
        main_genes = [g for a in main_nucleosomes for g in a.genes]
        other_nucleosomes = [a for c in self.other.chromosomes for a in c.nucleosomes]
        other_genes = [g for a in other_nucleosomes for g in a.genes]
        recombined_chromosomes = 0
        recombined_nucleosomes = 0
        recombined_genes = 0

        for c in genome.chromosomes:
            if c not in self.genome.chromosomes and c not in self.other.chromosomes:
                recombined_chromosomes += 1
            for a in c.nucleosomes:
                if a not in main_nucleosomes and a not in other_nucleosomes:
                    recombined_nucleosomes += 1
                for g in a.genes:
                    if g not in main_genes and g not in other_genes:
                        recombined_genes += 1
        assert recombined_chromosomes == 0
        assert recombined_nucleosomes == 0
        assert recombined_genes == 0

    def test_recombine_can_recombine_underlying_Chromosomes_Nucleosomes_and_Genes_with_different_names(self):
        genome = self.genome.recombine(
            self.other, [1], match_chromosomes=False, match_nucleosomes=False,
            match_genes=False
        )
        main_nucleosomes = [a for c in self.genome.chromosomes for a in c.nucleosomes]
        main_genes = [g for a in main_nucleosomes for g in a.genes]
        other_nucleosomes = [a for c in self.other.chromosomes for a in c.nucleosomes]
        other_genes = [g for a in other_nucleosomes for g in a.genes]
        recombined_chromosomes = 0
        recombined_nucleosomes = 0
        recombined_genes = 0
        for c in genome.chromosomes:
            if c not in self.genome.chromosomes and c not in self.other.chromosomes:
                recombined_chromosomes += 1
            for a in c.nucleosomes:
                if a not in main_nucleosomes and a not in other_nucleosomes:
                    recombined_nucleosomes += 1
                for g in a.genes:
                    if g not in main_genes and g not in other_genes:
                        recombined_genes += 1
        assert recombined_chromosomes >= 1
        assert recombined_nucleosomes >= 1
        assert recombined_genes >= 1

    def test_make_returns_random_Genome_with_random_Chromosomes_Nucleosomes_and_Genes(self):
        genome = genes.Genome.make(2, 3, 5, 7)
        assert type(genome) is genes.Genome
        assert len(genome.chromosomes) == 2
        assert all(len(c.nucleosomes) == 3 for c in genome.chromosomes)
        assert all(len(a.genes) == 5 for c in genome.chromosomes for a in c.nucleosomes)
        assert all(
            len(g.bases) == 7
            for c in genome.chromosomes
            for a in c.nucleosomes
            for g in a.genes
        )

        results: set[genes.Genome] = set()
        for _ in range(10):
            results.add(genes.Genome.make(3, 5, 10, 20))
        assert len(results) >= 8
        assert len(set([n.name for n in results])) >= 8
        assert len(set([tuple(n.chromosomes) for n in results])) >= 8
        assert len(set([tuple(c.nucleosomes) for n in results for c in n.chromosomes])) >= 8

    def test_make_with_base_factory_sets_bases_with_correct_type(self):
        float_factory = lambda: random()
        genome = genes.Genome.make(2, 3, 5, 7, base_factory=float_factory)
        assert all(
            type(b) is float
            for c in genome.chromosomes
            for a in c.nucleosomes
            for g in a.genes
            for b in g.bases
        )
        assert all(
            len(set(tuple(g.bases))) >= 6
            for c in genome.chromosomes
            for a in c.nucleosomes
            for g in a.genes
        )

        float_factory = lambda mult: random()*mult
        genome = genes.Genome.make(2, 3, 5, 7, base_factory=float_factory, factory_args=[-2])
        assert all(
            -2 <= b <= 0
            for c in genome.chromosomes
            for a in c.nucleosomes
            for g in a.genes
            for b in g.bases
        )

    def test_to_dict_and_from_dict_e2e(self):
        encoded = self.genome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Genome.from_dict(encoded)
        assert decoded == self.genome
        assert decoded is not self.genome

        float_factory = lambda: random()
        genome = genes.Genome.make(2, 3, 5, 7, base_factory=float_factory)
        encoded = genome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Genome.from_dict(encoded)
        assert decoded == genome
        assert decoded is not genome

        str_factory = lambda: genes.random_str(5)
        genome = genes.Genome.make(2, 3, 5, 7, base_factory=str_factory)
        encoded = genome.to_dict()
        assert type(encoded) is dict
        decoded = genes.Genome.from_dict(encoded)
        assert decoded == genome
        assert decoded is not genome


if __name__ == '__main__':
    unittest.main()
