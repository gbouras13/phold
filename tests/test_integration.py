"""
Integration tests for phold
Usage: pytest .

"""

# import
import os
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
database_dir = Path(f"{test_data}/phold_structure_foldseek_db")
model_dir = Path(f"{test_data}/model")
pdb_dir = Path(f"{test_data}/NC_043029_pdbs")
output_dir = Path(f"{test_data}/outputs")
output_dir.mkdir(parents=True, exist_ok=True)
predict_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_predict_gbk"
compare_pdb_dir: Path = f"{output_dir}/NC_043029_phold_compare_gbk_pdb"
compare_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_compare_gbk"

# functions_data = Path(f"{test_data}/functions_files")
# overall_data = Path(f"{test_data}/overall")
# meta_data = Path(f"{overall_data}/Meta_example")
# standard_data = Path(f"{overall_data}/Standard_examples")
# standard_data_output = Path(f"{standard_data}/SAOMS1_Output")
# stop_recoding_data = Path(f"{overall_data}/stop_recoding")
# custom_db = Path(f"{test_data}/custom_db/microvirus.h3m")
# custom_data = Path(f"{overall_data}/custom_examples")
# tmrna_data = Path(f"{overall_data}/tmRNA_example")
# AMR_data = Path(f"{overall_data}/AMR_example")
# CRISPR_data = Path(f"{overall_data}/CRISPR_example")
# VFDB_data = Path(f"{overall_data}/VFDB_example")
# genbank_data = Path(f"{overall_data}/genbank_examples")
logger.add(lambda _: sys.exit(1), level="ERROR")
threads = 1


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


temp_dir = Path(f"{test_data}/fake_out")


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


# def test_download(tmp_dir):
#     """test phold download"""
#     cmd = f"phold download -d {database_dir}"
#     exec_command(cmd)


def test_predict_genbank(tmp_dir):
    """test phold predict with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold predict -i {input_gbk} -o {predict_gbk_dir} -t {threads} -m {model_dir} --cpu -f"
    exec_command(cmd)

def test_compare_genbank(tmp_dir):
    """test phold predict with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_gbk_dir} --predictions_dir {predict_gbk_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)

def test_predict_pdb(tmp_dir):
    """test phold compare with pdbs input"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_pdb_dir} -t {threads} -d {database_dir} --pdb --pdb_dir {pdb_dir} -f"
    exec_command(cmd)



# class testFails(unittest.TestCase):
#     """Tests for fails"""

#     def test_dupe_header(self):
#         """tests that pharokka exits if a duplicate header is passed"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{standard_data}/dupe_header.fasta"
#             cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m"
#             exec_command(cmd)


# remove_directory(predict_gbk_dir)
# remove_directory(database_dir)
