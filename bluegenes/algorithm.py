from .errors import tert
from .genes import Gene, Allele, Chromosome, Genome
from random import choices
from typing import Any, Callable


def optimize_gene(measure_fitness: Callable[[Gene, int|float], int|float],
                  mutate_gene: Callable[[Gene], Gene],
                  initial_population: list[Gene] = None,
                  population_size: int = 100, gene_size: int = 10,
                  fitness_target: int|float = 1.0, max_iterations: int = 1000,
                  base_factory: Callable[[Any], int|float|str] = None,
                  factory_args: list[Any] = None,
                  factory_kwargs: dict[str, Any] = None,
                  gene_name: str = None, parents_per_generation: int = 10,
                  ) -> tuple[int, list[Gene]]:
    """Optimize a Gene given a measure_fitness function, a mutate_gene
        function, a population_size int, a gene_size int, a
        fitness_target float, and a max_iterations int. Supply
        base_factory to produce Gene bases other than random ints
        between 0 and 10, with optional factory_args and factory_kwargs
        which will be passed to each call of base_factory. Supply
        gene_name to assign the name to each Gene in the population.
        Returns the number of iterations and the final population.
    """
    tert(initial_population is None or all(type(p) is Gene for p in initial_population),
         "initial_population must be None or list[Gene]")
    population = [
        Gene.make(
            gene_size, base_factory=base_factory, factory_args=factory_args,
            factory_kwargs=factory_kwargs, name=gene_name
        )
        for _ in range(population_size)
    ] if initial_population is None else initial_population

    count = 0
    fitness_scores: list[tuple[int|float, Gene]] = [
        (measure_fitness(g, fitness_target), g)
        for g in population
    ]
    fitness_scores.sort(key=lambda fs: fs[0])
    fitness_scores.reverse()
    best_fitness = fitness_scores[0][0]

    while count < max_iterations and best_fitness < fitness_target:
        # breed parents at random proportional to their order by score
        parents = [fs[1] for fs in fitness_scores[:parents_per_generation]]
        children = [
            gene_child_from_parents(parents)
            for _ in range(population_size-len(parents))
        ]

        children = [mutate_gene(child) for child in children]
        population = [*children, *parents]
        fitness_scores: list[tuple[int|float, Gene]] = [
            (measure_fitness(g, fitness_target), g)
            for g in population
        ]
        fitness_scores.sort(key=lambda fs: fs[0])
        fitness_scores.reverse()
        best_fitness = fitness_scores[0][0]
        count += 1

    population = [fs[1] for fs in fitness_scores]
    return (count, population)

def gene_child_from_parents(parents: list[Gene]) -> Gene:
    """Select two parents at random semi-proportional to their order in
        the list. Recombine the two chosen parent Genes, and return the
        result.
    """
    weights = [len(parents[i:])/len(parents) for i in range(len(parents))]
    dad, mom = choices(parents, weights, k=2)
    return dad.recombine(mom)


def optimize_allele(measure_fitness: Callable[[Allele, int|float], int|float],
                  mutate_allele: Callable[[Allele], Allele],
                  initial_population: list[Allele] = None,
                  population_size: int = 100, allele_size: int = 2,
                  gene_size: int = 10, fitness_target: int|float = 1.0,
                  max_iterations: int = 1000,
                  base_factory: Callable[[Any], int|float|str] = None,
                  factory_args: list[Any] = None,
                  factory_kwargs: dict[str, Any] = None,
                  allele_name: str = None, parents_per_generation: int = 10,
                  ) -> tuple[int, list[Allele]]:
    """Optimize an Allele given a measure_fitness function, a
        mutate_allele function, a population_size int, a fitness_target
        float, and a max_iterations int. Supply base_factory to produce
        Gene bases other than random ints between 0 and 10,
        with optional factory_args and factory_kwargs which will be
        passed to each call of base_factory. Supply allele_name to
        assign the name to each generated Allele in the population.
        Supply an allele_size int and a gene_size int to customize
        generation of a random initial population, or supply
        initial_population list[Allele] to specify the initial
        population. Returns the number of iterations and the final
        population.
    """
    tert(initial_population is None or all(type(p) is Allele for p in initial_population),
         "initial_population must be None or list[Allele]")
    population = [
        Allele.make(
            allele_size, gene_size, base_factory=base_factory,
            factory_args=factory_args, factory_kwargs=factory_kwargs,
            name=allele_name
        )
        for _ in range(population_size)
    ] if initial_population is None else initial_population

    count = 0
    fitness_scores: list[tuple[int|float, Gene]] = [
        (measure_fitness(g, fitness_target), g)
        for g in population
    ]
    fitness_scores.sort(key=lambda fs: fs[0])
    fitness_scores.reverse()
    best_fitness = fitness_scores[0][0]

    while count < max_iterations and best_fitness < fitness_target:
        # breed parents at random proportional to their order by score
        parents = [fs[1] for fs in fitness_scores[:parents_per_generation]]
        children = [
            allele_child_from_parents(parents)
            for _ in range(population_size-len(parents))
        ]

        children = [mutate_allele(child) for child in children]
        population = [*children, *parents]
        fitness_scores: list[tuple[int|float, Gene]] = [
            (measure_fitness(g, fitness_target), g)
            for g in population
        ]
        fitness_scores.sort(key=lambda fs: fs[0])
        fitness_scores.reverse()
        best_fitness = fitness_scores[0][0]
        count += 1

    population = [fs[1] for fs in fitness_scores]
    return (count, population)

def allele_child_from_parents(parents: list[Allele]) -> Allele:
    """Select two parents at random semi-proportional to their order in
        the list. Recombine the two chosen parent Genes, and return the
        result.
    """
    weights = [len(parents[i:])/len(parents) for i in range(len(parents))]
    dad, mom = choices(parents, weights, k=2)
    return dad.recombine(mom)


def optimize_chromosome(
        measure_fitness: Callable[[Chromosome, int|float], int|float],
        mutate_chromosome: Callable[[Chromosome], Chromosome],
        initial_population: list[Chromosome] = None, population_size: int = 100,
        chromosome_size: int = 2, allele_size: int = 3, gene_size: int = 10,
        fitness_target: int|float = 1.0, max_iterations: int = 1000,
        base_factory: Callable[[Any], int|float|str] = None,
        factory_args: list[Any] = None,
        factory_kwargs: dict[str, Any] = None,
        chromosome_name: str = None, parents_per_generation: int = 10,
    ) -> tuple[int, list[Chromosome]]:
    """Optimize an Chromosome given a measure_fitness function, a
        mutate_chromosome function, a population_size int, a fitness_target
        float, and a max_iterations int. Supply base_factory to produce
        Gene bases other than random ints between 0 and 10,
        with optional factory_args and factory_kwargs which will be
        passed to each call of base_factory. Supply chromosome_name to
        assign the name to each generated Chromosome in the population.
        Supply an allele_size int and a gene_size int to customize
        generation of a random initial population, or supply
        initial_population list[Chromosome] to specify the initial
        population. Returns the number of iterations and the final
        population.
    """
    tert(initial_population is None or all(type(p) is Chromosome for p in initial_population),
         "initial_population must be None or list[Chromosome]")
    population = [
        Chromosome.make(
            chromosome_size, allele_size, gene_size, base_factory=base_factory,
            factory_args=factory_args, factory_kwargs=factory_kwargs,
            name=chromosome_name
        )
        for _ in range(population_size)
    ] if initial_population is None else initial_population

    count = 0
    fitness_scores: list[tuple[int|float, Gene]] = [
        (measure_fitness(g, fitness_target), g)
        for g in population
    ]
    fitness_scores.sort(key=lambda fs: fs[0])
    fitness_scores.reverse()
    best_fitness = fitness_scores[0][0]

    while count < max_iterations and best_fitness < fitness_target:
        # breed parents at random proportional to their order by score
        parents = [fs[1] for fs in fitness_scores[:parents_per_generation]]
        children = [
            chromosome_child_from_parents(parents)
            for _ in range(population_size-len(parents))
        ]

        children = [mutate_chromosome(child) for child in children]
        population = [*children, *parents]
        fitness_scores: list[tuple[int|float, Gene]] = [
            (measure_fitness(g, fitness_target), g)
            for g in population
        ]
        fitness_scores.sort(key=lambda fs: fs[0])
        fitness_scores.reverse()
        best_fitness = fitness_scores[0][0]
        count += 1

    population = [fs[1] for fs in fitness_scores]
    return (count, population)

def chromosome_child_from_parents(parents: list[Chromosome]) -> Chromosome:
    """Select two parents at random semi-proportional to their order in
        the list. Recombine the two chosen parent Genes, and return the
        result.
    """
    weights = [len(parents[i:])/len(parents) for i in range(len(parents))]
    dad, mom = choices(parents, weights, k=2)
    return dad.recombine(mom)
