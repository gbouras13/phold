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
    max_seqs: int,
    only_representatives: bool,
    ultra_sensitive: bool,
    extra_foldseek_params: str
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
        max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
        only_representatives (bool): turns off --cluster-search 1 parameter in foldseek
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
        extra_foldseek_params (str): Extra foldseek search params

    Returns:
        None
    """

    if ultra_sensitive:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --exhaustive-search"
    else:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --max-seqs {max_seqs}"

    if only_representatives is False:
        cmd += " --cluster-search 1"

    if extra_foldseek_params:
        cmd += f" {extra_foldseek_params}"

    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
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

    foldseek_createtsv = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createtsv {query_db} {target_db} {result_db}  {result_tsv}  ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createtsv)
