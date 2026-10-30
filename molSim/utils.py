import pandas as pd

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