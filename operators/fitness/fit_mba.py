from typing import Any, Dict, Tuple
from rdkit import Chem
from bblocks.bbmanager import BuildingBlockManager
from .fitness import FitnessOperator
from utils import lev_dist
from mapchiral import mapchiral 
from mxfp import mxfp

# Initialize building block manager and MXFP calculator.
bbmanager = BuildingBlockManager()
mxfp_calculator = mxfp.MXFPCalculator()
bbmanager = BuildingBlockManager()

class MBA(FitnessOperator):
    """Fitness operator that calculates fitness relative to known MBAs sequences

    Attributes:
        name(str): Unique identifier for this fitness operator
    """

    name = "mbafit"

    def __init__(self, w_map4c:float=0.5, w_mxfp:float=0.5, mxfp_scale:int=3000) -> None:
        super().__init__()
        self.w_map4c = w_map4c
        self.w_mxfp = w_mxfp 
        self.mxfp_scale = mxfp_scale 
        ## The idea is to at some point have an average similarity to all the mbas sequences
        # self.bb_sequence_mbas = {
        # 'lysocin_e': 'BB0033-BB0049-BB0029-BB0044-BB0122-BB0009-BB0049-BB0040-BB0039-BB0025-BB0013-BB0034',
        # }

    def process_query(self, query: str, query_format: str) -> Tuple:
        """
        Process the query to convert it into a molecular fingerprint.

        Depending on the query format, this method converts the query into a SMILES string,
        then creates an RDKit molecule and computes its fingerprint using the MXFP calculator.

        Args:
            query (str): The query input, which can be a SMILES string, a sequence, or a PDGA sequence.
            query_format (str): The format of the query. Expected values are 'smiles', 'sequence', or 'pdga_sequence'.

        Returns:
            The fingerprint of the query molecule for MAP4C and MXFP.

        Raises:
            ValueError: If an invalid query format is provided.
        """
        if query_format == 'smiles':
            _query = query
        elif query_format == 'sequence':
            _query = Chem.MolToSmiles(Chem.MolFromSequence(query, flavor=1))
        elif query_format == 'pdga_sequence':
            _query = bbmanager.seq_to_smiles(query)
        else:
            raise ValueError(f'Invalid query format: {query_format}')

        try:
            mol = Chem.MolFromSmiles(_query)

            # Compute both map4c and mxfp fingerprints
            mxfp_fp= mxfp_calculator.mxfp_from_mol(mol)
            map4c_fp= mapchiral.encode(mol)

        except Exception as e:
            raise ValueError(f"An excepcion ocurred for query: {_query}\n Raised excepction: {e}\n")

        return (map4c_fp, mxfp_fp)
    
    def process(self, sequence: str) -> Any:
        """Take any sequence object 

        Args:
            sequence (str): _description_
            translation_dict (Dict): _description_

        Returns:
            Any: _description_
        """
        seq = bbmanager.seq_to_smiles(sequence) 
        mol = Chem.MolFromSmiles(seq)

        # Compute both map4c and mxfp fingerprints
        mxfp_fp= mxfp_calculator.mxfp_from_mol(mol)
        map4c_fp= mapchiral.encode(mol)
        
        return  (map4c_fp, mxfp_fp)


    def fitness(self, individual: Any, query: Any) -> float:
        """Here we should check if the sequence is very similar to any MBA known, then
        we will favor those cyclations that according to the literature enhance activity while
        penalizing those who don't. 
        To get a sense of how similar they are to known MBAs we check:
            a) How close they are to the MBA sequence in Levenshtein distance
            b) How similar they are in MXFP.
        For consistency all of the fitness function will be treated as maximization optimization problems

        When a Ser/Thr is at position 1, side-chain cyclization (sSC) is often the "natural" and more
        active closure, while in other families the b-hydroxi 
        """
        ind_map4c, ind_mxfp = individual
        q_map4c, q_mxfp = query

        # distances we MINIMIZE
        jacc = mapchiral.jaccard_similarity(ind_map4c, q_map4c)   # in [0,1]
        d_map4c = 1.0 - float(jacc)                               # in [0,1]

        d_mxfp = float(abs(ind_mxfp - q_mxfp).sum())              # >= 0
        norm_mxfp = d_mxfp / (d_mxfp + self.mxfp_scale)           # in [0,1)

        # weighted sum (lower is better)
        return self.w_map4c * d_map4c + self.w_mxfp * norm_mxfp
