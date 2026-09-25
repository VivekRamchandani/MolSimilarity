import pytest
from utils.conversions import create_graph
import networkx as nx

def test_smiles_to_graph():
    smiles = "N1C2C=CC=CC=2CCC(Br)C1=O"
    g = nx.Graph()
    g.add_nodes_from([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    g.add_edges_from([
        (0, 1), (0, 11), (1, 2), (1, 6), 
        (2, 3), (3, 4), (4, 5), (5, 6), 
        (6, 7), (7, 8), (8, 9), (9, 10), 
        (9, 11), (11, 12)
    ])

    graph = create_graph().fromSmiles(smiles)
    assert nx.utils.graphs_equal(graph, g)

def test_smiles_to_graph_execption():
    invalid_smiles = "N1CC=CC=CC=2CCC(Br)C1=O"
    with pytest.raises(Exception):
        create_graph().fromSmiles(invalid_smiles)