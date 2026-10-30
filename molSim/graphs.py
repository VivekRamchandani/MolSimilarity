from rdkit import Chem
from rdkit.Chem import AllChem

import networkx as nx

import pyfglt.fglt as fg

from collections.abc import Callable
from typing import Concatenate
from functools import partial

import numpy as np


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

def _createProximityGraph(mol: Chem.rdchem.Mol, distance: float) -> nx.Graph:
    """Create proximity graph from 3d structure of molecule

    Args:
        mol (Chem.rdchem.Mol): molecule to create graph of
        distance (float): controls the distance to determine whether two atoms are in each other's proximity.

    Returns:
        graph (networkx.Graph): returns graph created using networkx
    """
    
    graph = nx.Graph()
    graph.add_nodes_from([i for i in range(0, mol.GetNumAtoms())])

    conf_id = AllChem.EmbedMolecule(mol)
    if conf_id == -1:
        raise Exception("Failed to generate conformer")

    conf = mol.GetConformer(conf_id)

    positions = conf.GetPositions()
    for atom_idx, atom_pos in enumerate(positions):
        neighbors = []
        for idx, other_pos in enumerate(positions[atom_idx+1:]):
            other_atom_idx = idx + atom_idx + 1
            dis = np.linalg.norm(atom_pos - other_pos)
            if dis <= distance:
                neighbors.append((atom_idx, other_atom_idx))
        if neighbors:
            graph.add_edges_from(neighbors)

    return graph

class GraphCreator():

    def __init__(self, creator_func: Callable[Concatenate[Chem.rdchem.Mol, ...], nx.Graph], **extra_args: dict[str, any]):
        """Graph Creator Factory

        Args:
            creator_func (Callable): Function that creates the graph using rdkit's `Mol` object.
            **extra_args (dict): extra arguments to be passed when calling creator_func with `Mol` object.
        """
        if extra_args:
            self.creator: Callable[[Chem.rdchem.Mol], nx.Graph] = partial(creator_func, **extra_args)
        else:
            self.creator = creator_func
        
    def createuUsingMolecule(self, mol: Chem.rdchem.Mol):
        if not isinstance(mol, Chem.rdchem.Mol):
            raise TypeError("`mol` not of type `Chem.rdchem.Mol`")

        mol = Chem.AddHs(mol)

        graph = self.creator(mol)

        num_rings = mol.GetRingInfo().NumRings()
        graph.graph["rings"] = num_rings

        return graph
        
    def createFromSmiles(self, smiles: str) -> Chem.rdchem.Mol:
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

        mol = Chem.AddHs(mol)

        graph = self.creator(mol)

        num_rings = mol.GetRingInfo().NumRings()
        graph.graph["rings"] = num_rings

        return graph

def graph_creator() -> GraphCreator:
    """Returns a `GraphCreator` object that creates graph of molecule structure.
    """
    return GraphCreator(_createGraph)

def proximity_graph_creator(distance: float) -> GraphCreator:
    """Return a `GraphCreator` object to make proximity graph of molecute structure
    
    Args:
        distance (float): controls the distance to determine whether two atoms are in each other's proximity. 
    """

    return GraphCreator(_createProximityGraph, distance=distance)


def embed_graph(graph: nx.Graph, normalize: bool = False) -> np.ndarray:
    """Embeds a graph

    Args:
        graph (networkx.Graph): graph to embed
        normalize (bool): normalize embeddings (default=False)
    
    Returns:
        embedding (np.ndarray): embeddings for the garph.
    """

    # Computing Graphlet Correlation Matrix from Graphlet Degree Vector
    gdv = fg.compute(graph)
    gcm_df = fg.compute_graphlet_correlation_matrix(gdv).fillna(0.0)
    # Converting dataframe to numpy array
    # and flattening 
    gcm = gcm_df.to_numpy()
    vec = np.r_[ *[gcm[i][i+1:] for i in range(0, 14)] ]

    if normalize:
        vec = vec / np.linalg.norm(vec)

    return vec