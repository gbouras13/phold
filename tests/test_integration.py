"""
Integration tests for phold
Usage: pytest .

"""

# import
import os

# this is needed due to a github actions CI issue https://github.com/pytorch/pytorch/issues/78490
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
import shutil

# import functions
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

# import functions


# test data
test_data = Path("tests/test_data")
database_dir = Path(f"{test_data}/phold_db")
pdb_dir = Path(f"{test_data}/NC_043029_pdbs")
output_dir = Path(f"{test_data}/outputs")
output_dir.mkdir(parents=True, exist_ok=True)
run_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_run_gbk"
predict_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_predict_gbk"
compare_pdb_dir: Path = f"{output_dir}/NC_043029_phold_compare_gbk_pdb"
compare_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_compare_gbk"
predict_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_predict_fasta"
compare_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_compare_fasta"
remote_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_remote_gbk"
remote_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_remote_fasta"
proteins_predict_dir: Path = f"{output_dir}/combined_truncated_phold_proteins_predict"
proteins_compare_dir: Path = f"{output_dir}/combined_truncated_phold_proteins_compare"


logger.add(lambda _: sys.exit(1), level="ERROR")
threads = 1


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


# @pytest.fixture(scope="session")
# def tmp_dir(tmpdir_factory):
#     return tmpdir_factory.mktemp("tmp")


temp_dir = Path(f"{test_data}/fake_out")

# the server can be down
run_remote = False


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_install():
    """test phold install"""
    cmd = f"phold install -d {database_dir} "
    exec_command(cmd)


def test_run_genbank():
    """test phold run with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_dir} -t {threads}  --cpu -d {database_dir} -f"
    exec_command(cmd)


def test_predict_genbank():
    """test phold predict with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold predict -i {input_gbk} -o {predict_gbk_dir} -t {threads}  --cpu -d {database_dir} -f"
    exec_command(cmd)


def test_compare_genbank():
    """test phold compare with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_gbk_dir} --predictions_dir {predict_gbk_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_compare_pdb():
    """test phold compare with pdbs input"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_pdb_dir} -t {threads} -d {database_dir} --pdb --pdb_dir {pdb_dir} -f"
    exec_command(cmd)


def test_predict_fasta():
    """test phold predict with fasta input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold predict -i {input_fasta} -o {predict_fasta_dir} -t {threads} -d {database_dir}  --cpu -f"
    exec_command(cmd)


def test_compare_fasta():
    """test phold compare with fasta input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
    cmd = f"phold compare -i {input_fasta} -o {compare_fasta_dir} --predictions_dir {predict_fasta_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_proteins_predict():
    """test phold proteins-predict"""
    input_fasta: Path = f"{test_data}/phanotate.faa"
    cmd = f"phold proteins-predict -i {input_fasta} -o {proteins_predict_dir} -t {threads} -d {database_dir} --cpu -f"
    exec_command(cmd)


def test_proteins_compare():
    """test phold proteins-compare"""
    input_fasta: Path = f"{test_data}/phanotate.faa"
    cmd = f"phold proteins-compare -i {input_fasta} --predictions_dir {proteins_predict_dir} -o {proteins_compare_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


if run_remote is True:

    def test_remote_genbank():
        """test phold remote with genbank input"""
        input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
        cmd = f"phold remote -i {input_gbk} -o {remote_gbk_dir} -t {threads} -d {database_dir} -f"
        exec_command(cmd)

    def test_remote_fasta():
        """test phold remote with fasta input"""
        input_fasta: Path = (
            f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
        )
        cmd = f"phold remote -i {input_fasta} -o {remote_fasta_dir} -t {threads} -d {database_dir} -f"
        exec_command(cmd)


# class testFails(unittest.TestCase):
#     """Tests for fails"""

#     def test_dupe_header(self):
#         """tests that pharokka exits if a duplicate header is passed"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{standard_data}/dupe_header.fasta"
#             cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m"
#             exec_command(cmd)


remove_directory(output_dir)


