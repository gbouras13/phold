#!/usr/bin/env python3
from pathlib import Path

from phold.utils.external_tools import ExternalTool


def run_foldseek_search(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    temp_db: Path,
    threads: int,
    logdir: Path,
    evalue: float,
    sensitivity: float,
) -> None:
    """
    Run a Foldseek search using given parameters.

    Args:
        query_db (Path): Path to the query database.
        target_db (Path): Path to the target database.
        result_db (Path): Path to store the result database.
        temp_db (Path): Path to store temporary files.
        threads (int): Number of threads to use for the search.
        logdir (Path): Path to the directory where logs will be stored.
        evalue (float): E-value threshold for the search.
        sensitivity (float): Sensitivity threshold for the search.

    Returns:
        None
    """

    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_search)


def create_result_tsv(
    query_db: Path, target_db: Path, result_db: Path, result_tsv: Path, logdir: Path
) -> None:
    """
    Create a TSV file containing the results of a Foldseek search.

    Args:
        query_db (Path): Path to the query database.
        target_db (Path): Path to the target database.
        result_db (Path): Path to the result database generated by the search.
        result_tsv (Path): Path to save the resulting TSV file.
        logdir (Path): Path to the directory where logs will be stored.

    Returns:
        None
    """
    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createtsv {query_db} {target_db} {result_db}  {result_tsv}  ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_search)
