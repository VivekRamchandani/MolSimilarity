from rdkit import Chem
import networkx as nx
from collections.abc import Callable
from typing import Concatenate
from functools import partial

def _createGraph(mol: Chem.rdchem.Mol) -> nx.Graph:
    """Create graph from 2d structure of molecule

    Args:
        mol (Chem.rdchem.Mol): molecule to create graph of

    Returns:
        graph (networkx.Graph): returns graph created using networkx
    """

    graph = nx.Graph()
    graph.add_nodes_from([i for i in range(0, mol.GetNumAtoms())])

    for atom in mol.GetAtoms():
        atom_id = atom.GetIdx()
        neighbors = []
        # Getting neighbor atoms
        for neighbor in atom.GetNeighbors():
            neighbor_id = neighbor.GetIdx()
            # Making sure to not add duplicate edges.
            if neighbor_id > atom_id:
                neighbors.append((atom_id, neighbor_id))
        if neighbors:
            graph.add_edges_from(neighbors)

    return graph

def _createProximityGraph(self, mol: Chem.rdchem.Mol, distance: float) -> nx.Graph:
    pass

class GraphCreator():

    def __init__(self, creator_func: Callable[Concatenate[Chem.rdchem.Mol, ...], nx.Graph], **extra_args: dict[str, any]):
        """Graph Creator Factory

        Args:
            creator_func (Callable): Function that creates the graph using rdkit's `Mol` object.
            **extra_args (dict): extra arguments to be passed when calling creator_func with `Mol` object.
        """
        if not extra_args:
            self.creator = creator_func
        else:
            self.creator: Callable[[Chem.rdchem.Mol], nx.Graph] = partial(creator_func, **creator_func)

        # Options
        self.add_Hs = False

    def add_hydrogen(self):
        self.add_Hs = True
        return self

    def usingMolecule(self, mol: Chem.rdchem.Mol):
        if not isinstance(mol, Chem.rdchem.Mol):
            raise TypeError("`mol` not of type `Chem.rdchem.Mol`")

        if self.add_Hs:
            mol = Chem.AddHs(mol)

        return self.creator(mol)
        
    def fromSmiles(self, smiles: str) -> Chem.rdchem.Mol:
        """Create molecule graph from smiles

        Args:
            smiles (str): SMILES string of molecule

        Returns:
            mol (Chem.rdchem.Mol): returns rdkit representation of Molecule
        """

        mol = Chem.MolFromSmiles(smiles)
        # Validate smiles
        if not mol:
            raise Exception("Error Invalid SMILES")

        if self.add_Hs:
            mol = Chem.AddHs(mol)

        return self.creator(mol)

def create_graph():
    return GraphCreator(_createGraph)

def create_proximity_graph(distance: float):
    return GraphCreator(_createProximityGraph, distance=distance)