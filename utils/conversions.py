from rdkit import Chem
import networkx as nx

class CreateGraph:
    """Networkx Graph creator class"""

    def __init__(self):
        pass

    def _molToGraph(self, mol: Chem.rdchem.Mol) -> nx.classes.graph.Graph:
        """Create graph from 2d structure of molecule

        Args:
            mol (Chem.rdchem.Mol): molecule to create graph of

        Returns:
            graph (networkx.classes.graph.Graph): returns graph created using networkx
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

    def fromSmiles(self, smiles: str) -> nx.classes.graph.Graph:
        """Create molecule graph from smiles

        Args:
            smiles (str): SMILES string of molecule

        Returns:
            graph (networkx.classes.graph.Graph): returns 2D strucutre of molecule as graph
        """

        mol = Chem.MolFromSmiles(smiles)
        # Validate smiles
        if not mol:
            raise Exception("Error Invalid SMILES")

        return self._molToGraph(mol)