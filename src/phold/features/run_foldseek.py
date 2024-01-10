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
    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_search)


# foldseek createtsv viji_fs_db/viji toy_foldseek_db/toy_prophage_db resultDB/result result.tsv


def create_result_tsv(
    query_db: Path, target_db: Path, result_db: Path, result_tsv: Path, logdir: Path
) -> None:
    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"createtsv {query_db} {target_db} {result_db}  {result_tsv}  ",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_search)
