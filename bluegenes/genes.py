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
    def make(cls, n_bases: int, max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {},
             name: str = None) -> Gene:
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
class Allele:
    name: str = field(default_factory=lambda: random_str(3))
    genes: list[Gene] = field(default_factory=list)

    def copy(self) -> Allele:
        """Returns an exact copy of the Allele, copying the underlying
            Genes as well.
        """
        return Allele(name=self.name, genes=[g.copy() for g in self.genes])

    def insert(self, index: int = None, gene: Gene = None, *,
               n_bases: int = None, max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Allele:
        """Inserts the gene at the index. If index is None, the gene is
            inserted at a random index. If gene is None, adds a random
            gene using Gene.make, passing the kwargs. Returns self for
            chaining operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        typert(index, int, "index")
        vert(index < len(self.genes), "index out of range")

        if gene is None:
            n_bases = n_bases or randint(0, max(len(g.bases) for g in self.genes))
            typert(n_bases, int, "n_bases")
            gene = Gene.make(
                n_bases, max_base_size, base_factory, factory_args,
                factory_kwargs
            )
        typert(gene, Gene, "gene")

        self.genes.insert(index, gene)
        return self

    def append(self, gene: Gene = None, /, *, n_bases: int = None,
               max_base_size: int = 10,
               base_factory: Callable[[None], int|float|str] = None,
               factory_args: list = [], factory_kwargs: dict = {}) -> Allele:
        """Adds a gene to the end of the allele. If gene is None, adds a
            random Gene using Gene.make and passing kwargs. Returns self
            for chaining operations.
        """
        return self.insert(
            len(self.genes), gene, n_bases=n_bases, max_base_size=max_base_size,
            base_factory=base_factory, factory_args=factory_args,
            factory_kwargs=factory_kwargs
        )

    def duplicate(self, index: int = None) -> Allele:
        """Duplicates the Gene at the index. If index is None, the Gene
            at a random index is duplicated. Returns self for chaining
            operations.
        """
        index = index if index is not None else randint(0, len(self.genes)-1)
        typert(index, int, "index")
        vert(0 <= index < len(self.genes), "index out of range")
        return self.insert(index, self.genes[index].copy())

    def delete(self, index: int = None) -> Allele:
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
                   factory_args: list = [], factory_kwargs: dict = {}) -> Allele:
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
                n_bases, max_base_size, base_factory, factory_args,
                factory_kwargs
            )
        typert(gene, Gene, "gene")

        self.genes[index] = gene
        return self

    def recombine(self, other: Allele, indices: list[int] = None,
                  recombine_genes: bool = True,
                  match_genes: bool = True) -> Allele:
        """Recombines with the other Allele, swapping at the given
            indices. If indices is None, between 1 and ceil(log(len(self.genes)))
            random indices will be chosen. Recombines individual Genes
            if recombine_genes is True. Recombines only Genes with
            matching names if match_genes is True. Returns the new
            Allele.
        """
        typert(other, Allele, "other")
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

        return Allele(name=name, genes=genes)

    @classmethod
    def make(cls, n_genes: int, n_bases: int, name: str = None, *,
             max_base_size: int = 10,
             base_factory: Callable[[None], int|float|str] = None,
             factory_args: list = [], factory_kwargs: dict = {}) -> Allele:
        """Makes and returns an Allele of randomized Genes."""
        genes = [
            Gene.make(
                n_bases=n_bases, max_base_size=max_base_size,
                base_factory=base_factory, factory_args=factory_args,
                factory_kwargs=factory_kwargs
            )
            for _ in range(n_genes)
        ]
        if name:
            return cls(name=name, genes=genes)
        return cls(genes=genes)

    def to_dict(self) -> dict:
        """Serialize the Allele to a dict."""
        return {
            self.name: [gene.to_dict() for gene in self.genes],
        }

    @classmethod
    def from_dict(cls, data: dict) -> Allele:
        """Deserialize an Allele from a dict."""
        for name, genes in data.items():
            unpacked = [Gene.from_dict(d) for d in genes]
            return cls(name=name, genes=unpacked)

    def __hash__(self) -> int:
        """Make Allele hashable."""
        return hash((self.name, hash(tuple(self.genes))))


@dataclass
class Chromosome:
    name: str = field(default_factory=lambda: random_str(2))
    alleles: list[Allele] = field(default_factory=list)


@dataclass
class Genome:
    name: str = field(default_factory=lambda: random_str(6))
    chromosomes: list[Chromosome] = field(default_factory=list)
