import argparse
from pathlib import Path
from . import utils
import sqlite3

def handle_createdb(args: argparse.Namespace) -> None:
    """Handler for 'createdb' subcommand."""
    db_path: Path = args.dbPath
    data_file: Path = args.file

    if db_path.exists():
        raise FileExistsError("Folder already exists.")

    db_path.mkdir()

    suffix = data_file.suffix

    match suffix:
        case ".smi" | ".smiles":
            df = utils.load_smi(str(data_file))
        case _:
            db_path.rmdir()
            raise Exception(f"File format {suffix} not supported !!!!!")

    conn = sqlite3.connect(db_path / "metadata.db")
    df.to_sql(name="metadata", con=conn)

def handle_querydb(args: argparse.Namespace) -> None:
    """Handler for 'querydb' subcommand."""
    db_path: Path = args.dbPath
    query_file: Path = args.file

    # TODO: Implement database query logic here


def build_parser() -> argparse.ArgumentParser:
    """Build argparser"""
    parser = argparse.ArgumentParser(
        prog="molSim",
        description="CLI tool for the molSim",
    )
    
    subparsers = parser.add_subparsers(
        dest="command",
        title="subcommands",
        description="Valid subcommands (use `molSim <command> --help` for more info)",
        required=True,
    )

    # Command: createdb
    parser_createdb = subparsers.add_parser(
        "createdb",
        help="Create a new database from a data file.",
        description="Creates a new database. Requires a folder path and a data file.",
    )
    parser_createdb.add_argument(
        "--dbPath",
        required=True,
        type=Path,
        help="Path to new database folder (e.g., /path/to/db)",
    )
    parser_createdb.add_argument(
        "--file",
        required=True,
        type=Path,
        help="Path to a molecule file (e.g., data.smi)",
    )
    parser_createdb.set_defaults(func=handle_createdb)

    # Command: querydb
    parser_querydb = subparsers.add_parser(
        "querydb",
        help="Run queries on an existing database.",
        description="Querying the database for similar molecules. Requires a DB folder path and a query file.",
    )
    parser_querydb.add_argument(
        "--dbPath",
        required=True,
        type=Path,
        help="Path to existing database folder (e.g., /path/to/db)",
    )
    parser_querydb.add_argument(
        "--file",
        required=True,
        type=Path,
        help="Path to a molecule file (e.g., query.smi)",
    )
    parser_querydb.set_defaults(func=handle_querydb)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Route execution to handler function
    args.func(args)


if __name__ == "__main__":
    main()