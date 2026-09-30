import pyfglt.fglt as fg
import numpy as np
import networkx as nx

def embed_graph(graph: nx.Graph) -> np.array:
    # Computing Graphlet Correlation Matrix from Graphlet Degree Vector
    gdv = fg.compute(graph)
    gcm_df = fg.compute_graphlet_correlation_matrix(gdv).fillna(0.0)
    # Converting dataframe to numpy array
    # and flattening 
    gcm = gcm_df.to_numpy()
    vec = np.r_[*[gcm[i][i+1:] for i in range(0, 14)]]

    return vec