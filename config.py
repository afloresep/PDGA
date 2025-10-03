"""
Configuration settings for the Peptide Design Genetic Algorithm (PDGA).
Update the values below to change the behavior of the algorithm.
"""

CONFIG = {
    # Query settings
    'query': 'CC[C@H](C)[C@H]1C(=O)N[C@H](C(=O)O[C@@H]([C@@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)NCC(=O)N([C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N1)CC2=CNC3=CC=CC=C32)CCC(=O)N)CCC(=O)O)CCCN=C(N)N)CC(C)C)CC4=CC=CC=C4)C)CO)CCCN=C(N)N)NC(=O)C[C@@H](CC(C)C)O)C)[C@@H](C)O',
    'query_format': 'smiles', 

    # Population settings
    'pop_size': 50,
    'pop_selection': 10,

    # Genetic operators parameters
    'mutation_ratio': 0.8,
    'cutoff': 0.3,
    'fitness_function': 'mbafit',
    'selection_strategy': 'greedy',
    'crossover_method': 'single_point',

    # Algorithm runtime settings
    'n_iterations': 10_000,
    'seed': 42,

    # Output settings
    'run_id': 'run',
    'maximize': False,

    # Fitness params
    'fitness_params': {
        'w_map4c': 0.7,
        'w_mxfp': 0.3,
        'mxfp_scale': 400.0,       
    },
}
