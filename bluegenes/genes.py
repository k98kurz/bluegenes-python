from __future__ import annotations
from .errors import tert, typert, vert
from dataclasses import dataclass, field
from math import ceil, log
from random import randint
from typing import Callable


alphanumerics = [
    *[chr(i) for i in range(48, 58)],
    *[chr(i) for i in range(65, 91)],
    *[chr(i) for i in range(97, 123)],
]


def random_str(size: int) -> str:
    """Returns a str of random alphanumeric chars."""
    l = len(alphanumerics)-1
    return "".join([
        alphanumerics[randint(0, l)] for _ in range(size)
    ])


@dataclass
class Gene:
    """Represents a gene comprised of coding bases with a name of some type."""
    name: str = field(default_factory=lambda: random_str(4))
    bases: list[int|float|str] = field(default_factory=list)

    def copy(self) -> Gene:
        """Returns an exact copy of the Gene."""
        return Gene(name=self.name, bases=[*self.bases])

    def insert(self, index: int = None, base: int|float|str = None) -> Gene:
        """Inserts the base at the index. If index is None, the base is
            inserted at a random index. If base is None, adds a random
            int base. Returns self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.bases)-1)
        base = base if base is not None else randint(0, max(self.bases))
        typert(index, int, "index")
        vert(index < len(self.bases), "index must be < len(bases)")
        typert(base, [int, float, str], "base")
        self.bases.insert(index, base)
        return self

    def append(self, base: int|float|str = None) -> Gene:
        """Adds a base to the end of the gene. If base is None, adds a
            random int base. Returns self for chaining operations.
        """
        base = base if base is not None else randint(0, max(self.bases))
        typert(base, [int, float, str], "base")
        self.bases.append(base)
        return self

    def insert_sequence(self, index: int = None,
                        sequence: list[int|float|str] = None) -> Gene:
        """Inserts the sequence at the index. If index is None, the
            sequence is inserted at a random index. If sequence is None,
            a random sequence in size between 1 and the current len of
            the gene bases will be inserted.
        """
        index = index if index is not None else randint(0, len(self.bases)-1)
        if sequence is None:
            size = randint(1, len(self.bases))
            max_base = max(self.bases)
            sequence = [randint(0, max_base) for _ in range(size)]
        typert(index, int, "index")
        vert(index < len(self.bases), "index must be < len(bases)")
        tert(type(sequence) is list, "sequence must be list[int|float|str]")
        tert(all(type(s) in (int, float, str) for s in sequence),
             "sequence must be list[int|float|str]")
        self.bases = [*self.bases[:index], *sequence, *self.bases[index:]]
        return self

    def delete(self, index: int = None) -> Gene:
        """Deletes the base at the index. If index is None, a random
            base is deleted. Returns self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.bases)-1)
        typert(index, int, "index")
        vert(index < len(self.bases), "index must be < len(bases)")
        del self.bases[index]
        return self

    def delete_sequence(self, index: int = None, size: int = None) -> Gene:
        """Deletes size bases beggining at the index. If index is None,
            a random index is used. If size is None, a random size is
            used. Returns self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.bases)-1)
        size = size if size is not None else randint(1, len(self.bases)-index)
        typert(index, int, "index")
        vert(index < len(self.bases), "index must be < len(bases)")
        typert(size, int, "size")
        vert(size > 0, "size must be > 0")
        del self.bases[index:index+size]
        return self

    def substitute(self, index: int = None, base: int|float|str = None) -> Gene:
        """Substitutes the base at the index with the given base. If
            index is None, a random index will be used. If base is None,
            a random int base will be used. Returns self for chaining
            operations.
        """
        index = index if index is not None else randint(0, len(self.bases)-1)
        base = base if base is not None else randint(0, max(self.bases))
        typert(index, int, "index")
        vert(index < len(self.bases), "index must be < len(bases)")
        typert(base, [int, float, str], "base")
        self.bases[index] = base
        return self

    def recombine(self, other: Gene, indices: list[int] = None) -> Gene:
        """Recombines with another gene at the given indexes. If indices
            is None, between 1 and ceil(log(len(self.bases))) random
            indices will be chosen. Returns the resultant Gene.
        """
        typert(other, Gene, "other")
        vert(len(other.bases) > 0, "other must have bases")
        max_size = min(len(self.bases), len(other.bases))
        max_swaps = ceil(log(max_size)) or 1
        tert(indices is None or type(indices) is list,
             "indices must be list[int] or None")
        if type(indices) is list:
            tert(all(type(i) is int for i in indices),
                 "indices must be list[int] or None")
            vert(len(indices) <= max_size,
                 f"can have at most {max_size} indices")
        else:
            swaps = randint(1, max_swaps)
            indices = list(set([randint(0, max_size-1) for _ in range(swaps)]))
            indices.sort()

        name = self.name
        if self.name != other.name:
            name_size = min(len(self.name), len(other.name))
            name_swap = randint(1, name_size-1)
            name = self.name[:name_swap] + other.name[name_swap:]

        bases = [*self.bases]
        swapped = False
        for i in indices:
            bases[i:] = self.bases[i:] if swapped else other.bases[i:]
            swapped = not swapped

        return Gene(name=name, bases=bases)

    @classmethod
    def make(cls, n_bases: int, name: str = None, *, max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {},) -> Gene:
        """Makes and returns a randomized Gene."""
        if base_factory is not None and callable(base_factory):
            bases = [
                base_factory(*factory_args, **factory_kwargs)
                for _ in range(n_bases)
            ]
        else:
            bases = [randint(0, max_base_size) for _ in range(n_bases)]
        if name:
            return cls(name=name, bases=bases)
        return cls(bases=bases)

    def to_dict(self) -> dict:
        """Serialize the Gene to a dict."""
        return {self.name: [*self.bases]}

    @classmethod
    def from_dict(cls, data: dict) -> Gene:
        """Deserialize a Gene from a dict."""
        for name, bases in data.items():
            return cls(name=name, bases=bases)

    def __hash__(self) -> int:
        """Make Gene hashable."""
        return hash((self.name, tuple(self.bases)))


@dataclass
class Nucleosome:
    name: str = field(default_factory=lambda: random_str(3))
    genes: list[Gene] = field(default_factory=list)

    def copy(self) -> Nucleosome:
        """Returns an exact copy of the Nucleosome, copying the underlying
            Genes as well.
        """
        return Nucleosome(name=self.name, genes=[g.copy() for g in self.genes])

    def insert(self, index: int = None, gene: Gene = None, *,
               n_bases: int = None, max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Nucleosome:
        """Inserts the gene at the index. If index is None, the gene is
            inserted at a random index. If gene is None, adds a random
            gene using Gene.make, passing the kwargs. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        typert(index, int, "index")
        vert(0 <= index <= len(self.genes), "index out of range")

        if gene is None:
            if n_bases is None:
                n_bases = randint(0, max(len(g.bases) for g in self.genes))
            typert(n_bases, int, "n_bases")
            gene = Gene.make(
                n_bases, max_base_size=max_base_size, base_factory=base_factory,
                factory_args=factory_args,factory_kwargs=factory_kwargs
            )
        typert(gene, Gene, "gene")

        self.genes.insert(index, gene)
        return self

    def append(self, gene: Gene = None, /, *, n_bases: int = None,
               max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Nucleosome:
        """Adds a gene to the end of the nucleosome. If gene is None, adds a
            random Gene using Gene.make and passing kwargs. Returns self
            for chaining operations.
        """
        return self.insert(
            len(self.genes), gene, n_bases=n_bases, max_base_size=max_base_size,
            base_factory=base_factory, factory_args=factory_args,
            factory_kwargs=factory_kwargs
        )

    def duplicate(self, index: int = None) -> Nucleosome:
        """Duplicates the Gene at the index. If index is None, the Gene
            at a random index is duplicated. Returns self for chaining
            operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        typert(index, int, "index")
        vert(0 <= index < len(self.genes), "index out of range")
        return self.insert(index, self.genes[index].copy())

    def delete(self, index: int = None) -> Nucleosome:
        """Deletes the Gene at the index. If index is None, a Gene is
            deleted at a random index. Returns self for chaining
            operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        del self.genes[index]
        return self

    def substitute(self, index: int = None, gene: Gene = None, *,
                   n_bases: int = None, max_base_size: int = 10,
                   base_factory: Callable[[None], int|float|str] = None,
                   factory_args: list = [], factory_kwargs: dict = {}) -> Nucleosome:
        """Substitutes the Gene at the index with the given gene. If
            index is None, a random index will be used. If gene is None,
            adds a random gene using Gene.make, passing the kwargs.
            Returns self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        typert(index, int, "index")
        vert(index < len(self.genes), "index out of range")

        if gene is None:
            n_bases = n_bases or randint(0, max(len(g.bases) for g in self.genes))
            typert(n_bases, int, "n_bases")
            gene = Gene.make(
                n_bases, max_base_size=max_base_size, base_factory=base_factory,
                factory_args=factory_args,factory_kwargs=factory_kwargs
            )
        typert(gene, Gene, "gene")

        self.genes[index] = gene
        return self

    def recombine(self, other: Nucleosome, indices: list[int] = None,
                  recombine_genes: bool = True,
                  match_genes: bool = False) -> Nucleosome:
        """Recombines with the other Nucleosome, swapping at the given
            indices. If indices is None, between 1 and ceil(log(len(self.genes)))
            random indices will be chosen. Recombines individual Genes
            if recombine_genes is True. Recombines only Genes with
            matching names if match_genes is True. Returns the new
            Nucleosome.
        """
        typert(other, Nucleosome, "other")
        vert(len(other.genes) > 0, "other must have genes")
        max_size = min(len(self.genes), len(other.genes))
        max_swaps = ceil(log(max_size)) or 1
        tert(indices is None or type(indices) is list,
             "indices must be list[int] or None")
        if type(indices) is list:
            tert(all(type(i) is int for i in indices),
                 "indices must be list[int] or None")
            vert(len(indices) <= max_size,
                 f"can have at most {max_size} indices")
        else:
            swaps = randint(0, max_swaps)
            indices = list(set([randint(0, max_size-1) for _ in range(swaps)]))
            indices.sort()

        name = self.name
        if self.name != other.name:
            name_size = min(len(self.name), len(other.name))
            name_swap = randint(1, name_size-1)
            name = self.name[:name_swap] + other.name[name_swap:]

        genes = [*self.genes]
        other_genes = [*other.genes]
        swapped = False
        for i in indices:
            genes[i:] = self.genes[i:] if swapped else other.genes[i:]
            other_genes[i:] = other.genes[i:] if swapped else self.genes[i:]
        genes = [g.copy() for g in genes]
        other_genes = [g.copy() for g in other_genes]

        if recombine_genes:
            for i in range(max_size):
                if genes[i].name == other_genes[i].name or not match_genes:
                    genes[i] = genes[i].recombine(other_genes[i])

        return Nucleosome(name=name, genes=genes)

    @classmethod
    def make(cls, n_genes: int, n_bases: int, name: str = None, *,
             max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {}) -> Nucleosome:
        """Makes and returns a Nucleosome of randomized Genes."""
        genes = [
            Gene.make(
                n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
            for _ in range(n_genes)
        ]
        if name:
            return cls(name=name, genes=genes)
        return cls(genes=genes)

    def to_dict(self) -> dict:
        """Serialize the Nucleosome to a dict."""
        return {
            self.name: [gene.to_dict() for gene in self.genes],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Nucleosome:
        """Deserialize a Nucleosome from a dict."""
        for name, genes in data.items():
            unpacked = [Gene.from_dict(d) for d in genes]
            return cls(name=name, genes=unpacked)

    def __hash__(self) -> int:
        """Make Nucleosome hashable."""
        return hash((self.name, tuple(self.genes)))


@dataclass
class Chromosome:
    name: str = field(default_factory=lambda: random_str(2))
    nucleosomes: list[Nucleosome] = field(default_factory=list)

    def copy(self) -> Chromosome:
        """Returns an exact copy of the Chromosome, copying the
            underlying Nucleosomes as well.
        """
        return Chromosome(name=self.name, nucleosomes=[a.copy() for a in self.nucleosomes])

    def insert(self, index: int = None, nucleosome: Nucleosome = None, *,
               n_genes: int = None, n_bases: int = None, max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Chromosome:
        """Inserts the nucleosome at the index. If index is None, the nucleosome
            is inserted at a random index. If nucleosome is None, adds a
            random nucleosome using Nucleosome.make, passing the kwargs. Returns
            self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.nucleosomes)-1)
        typert(index, int, "index")
        vert(0 <= index <= len(self.nucleosomes), "index out of range")

        if nucleosome is None:
            if n_genes is None:
                n_genes = randint(0, max(len(a.genes) for a in self.nucleosomes))
            if n_bases is None:
                n_bases = randint(0, max(len(g.bases) for a in self.nucleosomes for g in a.genes))
            typert(n_genes, int, "n_genes")
            typert(n_bases, int, "n_bases")
            nucleosome = Nucleosome.make(
                n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
        typert(nucleosome, Nucleosome, "nucleosome")

        self.nucleosomes.insert(index, nucleosome)
        return self

    def append(self, nucleosome: Nucleosome = None, *, n_genes: int = None,
               n_bases: int = None, max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Chromosome:
        """Adds a nucleosome to the end of the chromosome. If nucleosome is None,
            adds a random Nucleosome using Nucleosome.make and passing kwargs.
            Returns self for chaining operations.
        """
        return self.insert(
            len(self.nucleosomes), nucleosome, n_genes=n_genes, n_bases=n_bases,
            max_base_size=max_base_size, base_factory=base_factory,
            factory_args=factory_args, factory_kwargs=factory_kwargs
        )

    def duplicate(self, index: int = None) -> Chromosome:
        """Duplicates the Nucleosome at the index. If index is None, the
            Nucleosome at a random index is duplicated. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.nucleosomes)-1)
        typert(index, int, "index")
        vert(0 <= index < len(self.nucleosomes), "index out of range")
        return self.insert(index, self.nucleosomes[index].copy())

    def delete(self, index: int = None) -> Chromosome:
        """Deletes the Nucleosome at the index. If index is None, a Nucleosome
            is deleted at a random index. Returns self for chaining
            operations.
        """
        index = index if index is not None else randint(0, len(self.nucleosomes)-1)
        del self.nucleosomes[index]
        return self

    def substitute(self, index: int = None, nucleosome: Nucleosome = None, *,
                   n_genes: int = None, n_bases: int = None,
                   max_base_size: int = 10,
                   base_factory: Callable[[None], int|float|str] = None,
                   factory_args: list = [], factory_kwargs: dict = {}) -> Chromosome:
        """Substitutes the Nucleosome at the index with the given nucleosome. If
            index is None, a random index will be used. If nucleosome is
            None, adds a random nucleosome using Nucleosome.make, passing the
            kwargs. Returns self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.nucleosomes)-1)
        typert(index, int, "index")
        vert(index < len(self.nucleosomes), "index out of range")

        if nucleosome is None:
            if n_genes is None:
                n_genes = randint(0, max(len(a.genes) for a in self.nucleosomes))
            if n_bases is None:
                n_bases = randint(0, max(len(g.bases) for a in self.nucleosomes for g in a.genes))
            typert(n_genes, int, "n_genes")
            typert(n_bases, int, "n_bases")
            nucleosome = Nucleosome.make(
                n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
        typert(nucleosome, Nucleosome, "nucleosome")

        self.nucleosomes[index] = nucleosome
        return self

    def recombine(self, other: Chromosome, indices: list[int] = None,
                  recombine_nucleosomes: bool = True, match_nucleosomes: bool = False,
                  recombine_genes: bool = True, match_genes: bool = False
                  ) -> Chromosome:
        """Recombines with the other Chromosome, swapping at the given
            indices. If indices is None, between 1 and ceil(log(len(self.nucleosomes)))
            random indices will be chosen. Recombines individual Nucleosomes
            if recombine_allels is True. Recombines only Nucleosomes with
            matching names if match_nucleosomes is True. Recombines
            individual Genes if recombine_genes is True. Recombines only
            Genes with matching names if match_genes is True. Returns
            the new Chromosome.
        """
        typert(other, Chromosome, "other")
        vert(len(other.nucleosomes) > 0, "other must have nucleosomes")
        max_size = min(len(self.nucleosomes), len(other.nucleosomes))
        max_swaps = ceil(log(max_size)) or 1
        tert(indices is None or type(indices) is list,
             "indices must be list[int] or None")
        if type(indices) is list:
            tert(all(type(i) is int for i in indices),
                 "indices must be list[int] or None")
            vert(len(indices) <= max_size,
                 f"can have at most {max_size} indices")
        else:
            swaps = randint(0, max_swaps)
            indices = list(set([randint(0, max_size-1) for _ in range(swaps)]))
            indices.sort()

        name = self.name
        if self.name != other.name:
            name_size = min(len(self.name), len(other.name))
            name_swap = randint(1, name_size-1)
            name = self.name[:name_swap] + other.name[name_swap:]

        nucleosomes = [*self.nucleosomes]
        other_nucleosomes = [*other.nucleosomes]
        swapped = False
        for i in indices:
            nucleosomes[i:] = self.nucleosomes[i:] if swapped else other.nucleosomes[i:]
            other_nucleosomes[i:] = other.nucleosomes[i:] if swapped else self.nucleosomes[i:]
        nucleosomes = [g.copy() for g in nucleosomes]
        other_nucleosomes = [g.copy() for g in other_nucleosomes]

        if recombine_nucleosomes:
            for i in range(max_size):
                if nucleosomes[i].name == other_nucleosomes[i].name or not match_nucleosomes:
                    nucleosomes[i] = nucleosomes[i].recombine(
                        other_nucleosomes[i],
                        recombine_genes=recombine_genes,
                        match_genes=match_genes
                    )

        return Chromosome(name=name, nucleosomes=nucleosomes)

    @classmethod
    def make(cls, n_nucleosomes: int, n_genes: int, n_bases: int, name: str = None,
             *, max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {}) -> Chromosome:
        """Makes and returns a Chromosome of randomized Nucleosomes."""
        nucleosomes = [
            Nucleosome.make(
                n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
            for _ in range(n_nucleosomes)
        ]
        if name:
            return cls(name=name, nucleosomes=nucleosomes)
        return cls(nucleosomes=nucleosomes)

    def to_dict(self) -> dict:
        """Serialize the Chromosome to a dict."""
        return {
            self.name: [nucleosome.to_dict() for nucleosome in self.nucleosomes],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Chromosome:
        """Deserialize a Chromosome from a dict."""
        for name, nucleosomes in data.items():
            unpacked = [Nucleosome.from_dict(d) for d in nucleosomes]
            return cls(name=name, nucleosomes=unpacked)

    def __hash__(self) -> int:
        """Make Chromosome hashable."""
        return hash((self.name, tuple(self.nucleosomes)))


@dataclass
class Genome:
    name: str = field(default_factory=lambda: random_str(6))
    chromosomes: list[Chromosome] = field(default_factory=list)

    def copy(self) -> Genome:
        """Returns an exact copy of the Genome, copying the
            underlying Chromosome as well.
        """
        return Genome(name=self.name, chromosomes=[a.copy() for a in self.chromosomes])

    def insert(self, index: int = None, chromosome: Chromosome = None, *,
               n_nucleosomes: int = None, n_genes: int = None, n_bases: int = None,
               max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Chromosome:
        """Inserts the chromosome at the index. If index is None, the chromosome
            is inserted at a random index. If chromosome is None, adds a
            random chromosome using Chromosome.make, passing the kwargs. Returns
            self for chaining operations.
        """
        index = index if index is not None else randint(0, len(self.chromosomes)-1)
        typert(index, int, "index")
        vert(0 <= index <= len(self.chromosomes), "index out of range")

        if chromosome is None:
            if n_nucleosomes is None:
                n_nucleosomes = randint(0, max(len(c.nucleosomes) for c in self.chromosomes))
            if n_genes is None:
                n_genes = randint(0, max(
                    len(a.genes)
                    for c in self.chromosomes
                    for a in c.nucleosomes
                ))
            if n_bases is None:
                n_bases = randint(0, max(
                    len(g.bases)
                    for c in self.chromosomes
                    for a in c.nucleosomes
                    for g in a.genes
                ))
            typert(n_nucleosomes, int, "n_nucleosomes")
            typert(n_genes, int, "n_genes")
            typert(n_bases, int, "n_bases")
            chromosome = Chromosome.make(
                n_nucleosomes, n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
        typert(chromosome, Chromosome, "chromosome")

        self.chromosomes.insert(index, chromosome)
        return self

    def append(self, chromosome: Chromosome = None, *, n_nucleosomes: int = None,
               n_genes: int = None, n_bases: int = None, max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Genome:
        """Adds a chromosome to the end of the genome. If chromosome is
            None, adds a random Chromosome using Chromosome.make and
            passing kwargs. Returns self for chaining operations.
        """
        return self.insert(
            len(self.chromosomes), chromosome, n_nucleosomes=n_nucleosomes,
            n_genes=n_genes, n_bases=n_bases, max_base_size=max_base_size,
            base_factory=base_factory, factory_args=factory_args,
            factory_kwargs=factory_kwargs
        )

    def duplicate(self, index: int = None) -> Genome:
        """Duplicates the Chromosome at the index. If index is None, the
            Chromosome at a random index is duplicated. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.chromosomes)-1)
        typert(index, int, "index")
        vert(0 <= index < len(self.chromosomes), "index out of range")
        return self.insert(index, self.chromosomes[index].copy())

    def delete(self, index: int = None) -> Genome:
        """Deletes the Chromosome at the index. If index is None, a
            Chromosome is deleted at a random index. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.chromosomes)-1)
        del self.chromosomes[index]
        return self

    def substitute(self, index: int = None, chromosome: Chromosome = None, *,
                   n_nucleosomes: int = None, n_genes: int = None,
                   n_bases: int = None, max_base_size: int = 10,
                   base_factory: Callable[[None], int|float|str] = None,
                   factory_args: list = [], factory_kwargs: dict = {}) -> Genome:
        """Substitutes the Chromosome at the index with the given
            chromosome. If index is None, a random index will be used.
            If chromosome is None, adds a random chromosome using
            Chromosome.make, passing the kwargs. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.chromosomes)-1)
        typert(index, int, "index")
        vert(index < len(self.chromosomes), "index out of range")

        if chromosome is None:
            if n_nucleosomes is None:
                n_nucleosomes = randint(0, max(len(c.nucleosomes) for c in self.chromosomes))
            if n_genes is None:
                n_genes = randint(0, max(
                    len(a.genes)
                    for c in self.chromosomes
                    for a in c.nucleosomes
                ))
            if n_bases is None:
                n_bases = randint(0, max(
                    len(g.bases)
                    for c in self.chromosomes
                    for a in c.nucleosomes
                    for g in a.genes
                ))
            typert(n_nucleosomes, int, "n_nucleosomes")
            typert(n_genes, int, "n_genes")
            typert(n_bases, int, "n_bases")
            chromosome = Chromosome.make(
                n_nucleosomes, n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
        typert(chromosome, Chromosome, "chromosome")

        self.chromosomes[index] = chromosome
        return self

    def recombine(self, other: Genome, indices: list[int] = None,
                  recombine_chromosomes: bool = True,
                  match_chromosomes: bool = False,
                  recombine_nucleosomes: bool = True, match_nucleosomes: bool = False,
                  recombine_genes: bool = True, match_genes: bool = False
                  ) -> Genome:
        """Recombines with the other Genome, swapping at the given
            indices. If indices is None, between 1 and
            ceil(log(len(self.chromosomes))) random indices will be
            chosen. Recombines individual Chromosomes
            if recombine_chromosomes is True. Recombines only
            Chromosomes with matching names if match_chromosomes is True.
            Recombines individual Nucleosomes if recombine_nucleosomes is True.
            Recombines only Nucleosomes with matching names if match_nucleosomes
            is True. Recombines individual Genes if recombine_genes is
            True. Recombines only Genes with matching names if
            match_genes is True. Returns the new Genome.
        """
        typert(other, Genome, "other")
        vert(len(other.chromosomes) > 0, "other must have chromosomes")
        max_size = min(len(self.chromosomes), len(other.chromosomes))
        max_swaps = ceil(log(max_size)) or 1
        tert(indices is None or type(indices) is list,
             "indices must be list[int] or None")
        if type(indices) is list:
            tert(all(type(i) is int for i in indices),
                 "indices must be list[int] or None")
            vert(len(indices) <= max_size,
                 f"can have at most {max_size} indices")
        else:
            swaps = randint(0, max_swaps)
            indices = list(set([randint(0, max_size-1) for _ in range(swaps)]))
            indices.sort()

        name = self.name
        if self.name != other.name:
            name_size = min(len(self.name), len(other.name))
            name_swap = randint(1, name_size-1)
            name = self.name[:name_swap] + other.name[name_swap:]

        chromosomes = [*self.chromosomes]
        other_chromosomes = [*other.chromosomes]
        swapped = False
        for i in indices:
            chromosomes[i:] = self.chromosomes[i:] if swapped else other.chromosomes[i:]
            other_chromosomes[i:] = other.chromosomes[i:] if swapped else self.chromosomes[i:]
        chromosomes = [g.copy() for g in chromosomes]
        other_chromosomes = [g.copy() for g in other_chromosomes]

        if recombine_chromosomes:
            for i in range(max_size):
                if chromosomes[i].name == other_chromosomes[i].name or not match_chromosomes:
                    chromosomes[i] = chromosomes[i].recombine(
                        other_chromosomes[i],
                        recombine_nucleosomes=recombine_nucleosomes,
                        recombine_genes=recombine_genes,
                        match_nucleosomes=match_nucleosomes,
                        match_genes=match_genes
                    )

        return Genome(name=name, chromosomes=chromosomes)

    @classmethod
    def make(cls, n_chromosomes: int, n_nucleosomes: int, n_genes: int,
             n_bases: int, name: str = None, *, max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {}) -> Genome:
        """Makes and returns a Genome of randomized Chromosomes."""
        chromosomes = [
            Chromosome.make(
                n_nucleosomes, n_genes, n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
            for _ in range(n_chromosomes)
        ]
        if name:
            return cls(name=name, chromosomes=chromosomes)
        return cls(chromosomes=chromosomes)

    def to_dict(self) -> dict:
        """Serialize the Genome to a dict."""
        return {
            self.name: [chromosome.to_dict() for chromosome in self.chromosomes],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Genome:
        """Deserialize a Genome from a dict."""
        for name, chromosomes in data.items():
            unpacked = [Chromosome.from_dict(c) for c in chromosomes]
            return cls(name=name, chromosomes=unpacked)

    def __hash__(self) -> int:
        """Make Genome hashable."""
        return hash((self.name, tuple(self.chromosomes)))
