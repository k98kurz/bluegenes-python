from .genes import Gene, Allele, Chromosome, Genome
from random import choices
from typing import Any, Callable


def optimize_gene(measure_fitness: Callable[[Gene], int|float],
                  mutate_gene: Callable[[Gene], Gene],
                  population_size: int = 100, gene_size: int = 10,
                  fitness_target: int|float = 1.0, max_iterations: int = 1000,
                  base_factory: Callable[[Any], int|float|str] = None,
                  factory_args: list[Any] = None,
                  factory_kwargs: dict[str, Any] = None,
                  gene_name: str = None, parents_per_generation: int = 10,
                  ) -> tuple[int, list[Gene]]:
    """Optimize a gene given a measure_fitness function, a mutat_gene
        function, a population_size int, a gene_size int, a
        fitness_target float, and a max_iterations int. Supply
        base_factory to produce Gene bases other than random ints
        between 0 and 10, with optional factory_args and factory_kwargs
        which will be passed to each call of base_factory. Supply
        gene_name to assign the name to each Gene in the population.
        Returns the number of iterations and the final population.
    """
    population = [
        Gene.make(
            gene_size, base_factory=base_factory, factory_args=factory_args,
            factory_kwargs=factory_kwargs, name=gene_name
        )
        for _ in range(population_size)
    ]

    count = 0
    fitness_scores: tuple[int|float, Gene] = [
        (measure_fitness(g), g)
        for g in population
    ]
    fitness_scores.sort(lambda fs: fs[0])
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
        fitness_scores: tuple[int|float, Gene] = [
            (measure_fitness(g), g)
            for g in population
        ]
        fitness_scores.sort(lambda fs: fs[0])
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
