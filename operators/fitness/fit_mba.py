from typing import Any, Dict
from rdkit import Chem
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator
from utils import lev_dist


bbmanager = BuildingBlockManager()


class MBA(FitnessOperator):
 
    name = "mbafit"
    def __init__(self) -> None:
        super().__init__()

    def process_query(self, query: str, query_format: str) -> Any:
        return super().process_query(query, query_format)
    
    def process(self, sequence: str, translation_dict: Dict) -> Any:
        return super().process(sequence, translation_dict)
    
    def fitness(self, individual: Any, query: Any) -> float:
        return super().fitness(individual, query)
        """Here we should check if the sequence is very similar to any MBA known, then
        we will favor those cyclations that according to the literature enhance activity while
        penalizing those who don't. 
        To get a sense of how similar they are to known MBAs we check:
            a) How close they are to the MBA sequence in Levenshtein distance
            b) How similar they are in MXFP.


        When a Ser/Thr is at position 1, side-chain cyclization (sSC) is often the "natural" and more
        active closure, while in other families the b-hydroxi 

        """

