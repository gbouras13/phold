#!/usr/bin/env python3
from pathlib import Path
from typing import Optional

from pholdlib.prostt5.device import cuda_visible_devices_value, parse_gpus

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
    ultra_sensitive: bool,
    extra_foldseek_params: str,
    foldseek_gpu: bool,
    structures: bool,
    clustered_db: bool,
    gpus: Optional[str] = None,
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
        ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
        extra_foldseek_params (str): Extra foldseek search params
        foldseek_gpu (bool): Run Foldseek-GPU with accelerate ungapped prefilter
        structures (bool): Run Foldseek with structures, not ProstT5 3Dis
        clustered_db (bool): Run Foldseek with clustered DB (for benchmarking)
        gpus (Optional[str]): Comma-separated CUDA indices (e.g. "0,2") to
            restrict foldseek's GPU prefilter to a subset of devices. When
            ``foldseek_gpu`` is True and this resolves to ≥1 CUDA device,
            the foldseek subprocess gets ``CUDA_VISIBLE_DEVICES`` set
            accordingly. None = use all visible CUDA GPUs (foldseek default).
            Ignored when ``foldseek_gpu`` is False.

    Returns:
        None
    """

    if ultra_sensitive:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --exhaustive-search"
    else:
        cmd = f"search {query_db} {target_db} {result_db} {temp_db} --threads {str(threads)} -e {evalue} -s {sensitivity} --max-seqs {max_seqs}"

    # support foldseek gpu only for the regular DB search for now
    if foldseek_gpu:
        cmd = f"search {query_db} {target_db}_gpu {result_db} {temp_db} --threads {str(threads)} -e {evalue}  --gpu 1 --prefilter-mode 1 --max-seqs {max_seqs}"

    if extra_foldseek_params:
        cmd += f" {extra_foldseek_params}"

    # need -a 1 to compute the alignment so tmscore and lddt can be output (if using --structures)
    if structures:
        cmd += f" -a 1"

    if clustered_db:
        cmd += f" --cluster-search 1"

    # Build optional env for multi-GPU foldseek. Only applies when GPU mode is
    # on; foldseek selects devices via CUDA_VISIBLE_DEVICES (per its README).
    env = None
    if foldseek_gpu and gpus is not None:
        # parse_gpus validates indices against device_count and errors fast.
        devices = parse_gpus(cpu=False, gpus=gpus)
        cvd = cuda_visible_devices_value(devices)
        if cvd is not None:
            env = {"CUDA_VISIBLE_DEVICES": cvd}

    foldseek_search = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
        env=env,
    )

    ExternalTool.run_tool(foldseek_search)


# def run_foldseek_align(
#     query_db: Path,
#     target_db: Path,
#     result_db: Path,
#     temp_db: Path,
#     aln_db: Path,
#     threads: int,
#     logdir: Path,
#     foldseek_gpu: bool
# ) -> None:
#     """
#     Run a Foldseek align using given parameters.

#     foldseek align test_str/foldseek_db/phold ../../phold_db_v1_updated_annots_no_phrog_hits/all_phold_structures test_str/result_db/result_db alnNew -a

#     Args:
#         query_db (Path): Path to the query database.
#         target_db (Path): Path to the target database.
#         result_db (Path): Path to store the result database.
#         temp_db (Path): Path to store temporary files.
#         threads (int): Number of threads to use for the search.
#         logdir (Path): Path to the directory where logs will be stored.
#         evalue (float): E-value threshold for the search.
#         sensitivity (float): Sensitivity threshold for the search.
#         max_seqs (int): Maximum results per query sequence allowed to pass the prefilter for foldseek.
#         ultra_sensitive (bool): Whether to skip foldseek prefilter for maximum sensitivity
#         extra_foldseek_params (str): Extra foldseek search params
#         foldseek_gpu (bool): Run Foldseek-GPU with accelerate ungapped prefilter

#     Returns:
#         None
#     """

#     if foldseek_gpu:
#         target_db = f"{target_db}_gpu"

#     cmd = f"align {query_db} {target_db} {result_db} {aln_db} --threads {str(threads)} -a"


#     foldseek_search = ExternalTool(
#         tool="foldseek",
#         input=f"",
#         output=f"",
#         params=f"{cmd}",
#         logdir=logdir,
#     )

#     ExternalTool.run_tool(foldseek_search)


def create_result_tsv(
    query_db: Path,
    target_db: Path,
    result_db: Path,
    result_tsv: Path,
    logdir: Path,
    foldseek_gpu: bool,
    structures: bool,
    threads: int,
) -> None:
    """
    Create a TSV file containing the results of a Foldseek search.

    Args:
        query_db (Path): Path to the query database.
        target_db (Path): Path to the target database.
        result_db (Path): Path to the result database generated by the search.
        result_tsv (Path): Path to save the resulting TSV file.
        logdir (Path): Path to the directory where logs will be stored.
        foldseek_gpu (bool): Run Foldseek-GPU with accelerate ungapped prefilter
        structures (bool): Whether structures were input (not ProstT5)
        threads (int): Number of threads to use.

    Returns:
        None
    """
    if structures:
        format_string = "--format-output query,target,bits,fident,evalue,qstart,qend,qlen,tstart,tend,tlen,alntmscore,lddt"
    else:
        format_string = "--format-output query,target,bits,fident,evalue,qstart,qend,qlen,tstart,tend,tlen"
    if foldseek_gpu:
        target_db = f"{target_db}_gpu"

    cmd = f"convertalis {query_db} {target_db} {result_db} {result_tsv} {format_string} --threads {threads}"

    foldseek_createtsv = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_createtsv)
