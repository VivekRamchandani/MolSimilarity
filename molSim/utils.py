import pandas as pd
from rdkit import Chem
from . import proximity_graph_creator, graph_creator, embed_graph
import numpy as np
import faiss

def load_smi(file: str) -> pd.DataFrame:
    """Read .smi file format and laod smiles into dataframe

    Args:
        file (str): filepath
    """
    df = pd.read_csv(file, sep=r"\s+", names=["smiles", "id"])
    df = df.dropna()

    if not df.shape[0]:
        raise Exception("Empty file !??")

    df = df.reindex(["id", "smiles"], axis="columns")
    return df

def load_smi(file: str) -> Chem.SDMolSupplier:
    """Returns rdkit's molecule reader from .sdf file

    Args:
        file (str): filepath
    """
    pass

def faiss_index_from_smiles(df: pd.DataFrame, path: str):
    """Creates faiss index from smiles and writes at given path
    
    Args:
        df (pd.DataFrame): dataframe containing smiles
        path (str): folder path to store faiss index
    """

    creator = proximity_graph_creator(distance=2.25)

    rows = len(df)

    df["embeddings"] = None
    total_failed = 0

    for i in range(0, rows, 10000):
        end = min(i + 10000, rows)
        # graphs = np.empty(end - i, dtype=object)
        for idx in range(i, end):
            graph = creator.createFromSmiles(df.loc[idx, "smiles"])

            # Failed creating 
            if graph == None:
                total_failed += 1
                continue

            df.at[i, "embeddings"] = embed_graph(graph)

    df.dropna(inplace=True, ignore_index=True)

    x = np.stack(df["embeddings"].to_numpy())

    # Index config
    dim = x.shape[1]    # vector dimension 105
    m = 7               # 7 subvectors (each subvector dimension 15)
    nbits = 6           # 64 centroid for each subspace
    nlist = 512         # Max centroid we need 448 (7 x 64) 

    quantizer = faiss.IndexFlatL2(dim)
    index = faiss.IndexIVFPQ(quantizer, dim, nlist, m, nbits)

    index.train(x)
    if not index.is_trained():
        raise Exception("What the helllll !!!!")
    index.add(x)

    faiss.write_index(index, path / "index")