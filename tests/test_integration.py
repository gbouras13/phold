"""
Integration tests for phold

# to run pytest without remote and no gpu
pytest .

# to run with remote and no gpu
pytest --run_remote .  

# to run with remote and with gpu
pytest --run_remote  --gpu_available .

# to run with 8 threads 
pytest --run_remote  --gpu_available --threads 8 .

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
database_dir = Path(f"{test_data}/phold_db")
ext_database_dir = Path(f"{test_data}/phold_db_3M16")
pdb_dir = Path(f"{test_data}/NC_043029_pdbs")
cif_dir = Path(f"{test_data}/NC_043029_af3_cif")
output_dir = Path(f"{test_data}/outputs")
dummy_custom_db = Path(f"{test_data}/dummy_custom_db/dummy_db")
output_dir.mkdir(parents=True, exist_ok=True)
run_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_run_gbk"
run_gbk_hyps: Path = f"{output_dir}/NC_043029_hyps"
run_gbk_ncbi: Path = f"{output_dir}/run_gbk_ncbi_gbk"
run_gbk_bakta: Path = f"{output_dir}/run_gbk_ncbi_bakta"
run_gbk_pharokka_1_4_1_dir: Path = f"{output_dir}/NC_043029_pharokka1.4.1_gbk"
run_gbk_pharokka_1_4_1_dir_extra: Path = f"{output_dir}/NC_043029_pharokka1.4.1_gbk_extra"
run_gbk_pharokka_1_4_1_dir_custom: Path = f"{output_dir}/NC_043029_pharokka1.4.1_gbk_custom"
run_gbk_long_header_dir: Path = f"{output_dir}/long_header_gbk"
run_fasta_long_header_dir: Path = f"{output_dir}/long_header_fasta"
run_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_run_fasta"
run_fasta_efam_dir: Path = f"{output_dir}/KF_efam_phold_run_fasta"
run_fasta_netflax_dir: Path = f"{output_dir}/WP_netflax_phold_run_fasta"
# run_fasta_depolymerase_dir: Path = f"{output_dir}/PP_depolymerase_phold_run_fasta"
predict_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_predict_gbk"
save_embeddings_predict_gbk_dir: Path = (
    f"{output_dir}/combined_truncated_phold_predict_save_embeddings_gbk"
)
compare_pdb_dir: Path = f"{output_dir}/NC_043029_phold_compare_gbk_pdb"
compare_pdb_dir_custom: Path = f"{output_dir}/NC_043029_phold_compare_gbk_pdb_custom"
compare_cif_dir: Path = f"{output_dir}/NC_043029_phold_compare_gbk_pdb"
compare_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_compare_gbk"
predict_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_predict_fasta"
compare_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_compare_fasta"
remote_gbk_dir: Path = f"{output_dir}/combined_truncated_phold_remote_gbk"
remote_fasta_dir: Path = f"{output_dir}/combined_truncated_phold_remote_fasta"
proteins_predict_dir: Path = f"{output_dir}/combined_truncated_phold_proteins_predict"
proteins_compare_dir: Path = f"{output_dir}/combined_truncated_phold_proteins_compare"
proteins_compare_pdb_dir: Path = f"{output_dir}/NC_043029_phold_proteins_compare_pdb"
proteins_compare_cif_dir: Path = f"{output_dir}/NC_043029_phold_proteins_compare_cif"
proteins_compare_cif_dir_custom: Path = f"{output_dir}/NC_043029_phold_proteins_compare_cif_custom_db"
plots_dir: Path = f"{output_dir}/plot_output"


logger.add(lambda _: sys.exit(1), level="ERROR")
# threads = 1

def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

@pytest.fixture(scope="session")
def gpu_available(pytestconfig):
    return pytestconfig.getoption("gpu_available")

@pytest.fixture(scope="session")
def run_remote(pytestconfig):
    return pytestconfig.getoption("run_remote")


@pytest.fixture(scope="session")
def threads(pytestconfig):
    return pytestconfig.getoption("threads")


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


def test_install(threads):
    """test phold install"""
    cmd = f"phold install -d {database_dir} -t {threads}"
    exec_command(cmd)

def test_install_extended(threads):
    """test phold install"""
    cmd = f"phold install -d {ext_database_dir} -t {threads} --extended_db"
    exec_command(cmd)

def test_run_genbank(gpu_available, threads):
    """test phold run with genbank input"""
    input_gbk: Path = f"{test_data}/NC_043029_pharokka1.4.1.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_pharokka_1_4_1_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_hyps(gpu_available, threads):
    """test phold run with --hyps"""
    input_gbk: Path = f"{test_data}/NC_043029_pharokka1.4.1.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_hyps} -t {threads} -d {database_dir} -f --hyps"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_genbank_ncbi(gpu_available, threads):
    """test phold run with genbank input from NCBI"""
    input_gbk: Path = f"{test_data}/NC_043029_ncbi.gbk"
    cmd = (
        f"phold run -i {input_gbk} -o {run_gbk_ncbi} -t {threads} -d {database_dir} -f"
    )
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_genbank_bakta(gpu_available, threads):
    """test phold run with genbank input from bakta"""
    input_gbk: Path = f"{test_data}/NC_043029_bakta.gbk"
    cmd = (
        f"phold run -i {input_gbk} -o {run_gbk_bakta} -t {threads} -d {database_dir} -f"
    )
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_genbank_extra_foldseek_params(gpu_available, threads):
    """test phold run with --extra_foldseek_params"""
    input_gbk: Path = f"{test_data}/NC_043029_pharokka1.4.1.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_pharokka_1_4_1_dir_extra} -t {threads} -d {database_dir} --extra_foldseek_params \"--cov-mode 2\" -f "
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_custom_db(gpu_available, threads):
    """test phold run with --custom_db"""
    input_gbk: Path = f"{test_data}/NC_043029_pharokka1.4.1.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_pharokka_1_4_1_dir_custom} -t {threads} -d {database_dir} --custom_db {dummy_custom_db} -f "
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_run_genbank_old_pharokka(gpu_available, threads):
    """test phold run with genbank input from pharokka prior to v1.5.0 no transl_table field (#34)"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_genbank_long_header(gpu_available, threads):
    """test phold run with pharokka genbank with large ID/locus tag (over 54 chars)"""
    input_gbk: Path = f"{test_data}/long_header.gbk"
    cmd = f"phold run -i {input_gbk} -o {run_gbk_long_header_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_fasta_long_header(gpu_available, threads):
    """test phold run with FASTA with large header"""
    input_fasta: Path = f"{test_data}/long_header.fasta"
    cmd = f"phold run -i {input_fasta} -o {run_fasta_long_header_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_fasta(gpu_available, threads):
    """test phold run with genbank input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
    cmd = f"phold run -i {input_fasta} -o {run_fasta_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_efam(gpu_available, threads):
    """test phold run with a tophit to efam"""
    input_fasta: Path = f"{test_data}/KF623293.1_subset_efam.fasta"
    cmd = f"phold run -i {input_fasta} -o {run_fasta_efam_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_run_netflax(gpu_available, threads):
    """test phold run with a tophit to netflax"""
    input_fasta: Path = f"{test_data}/WP_006719989_subset_test.fasta "
    cmd = f"phold run -i {input_fasta} -o {run_fasta_netflax_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_predict_genbank(gpu_available, threads):
    """test phold predict with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold predict -i {input_gbk} -o {predict_gbk_dir} -t {threads}  -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_predict_save_embeddings(gpu_available, threads):
    """test phold predict with genbank input and save embeddings"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold predict -i {input_gbk} -o {save_embeddings_predict_gbk_dir} -t {threads}  -d {database_dir} -f --save_per_residue_embeddings --save_per_protein_embeddings"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_compare_genbank(threads):
    """test phold compare with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_gbk_dir} --predictions_dir {predict_gbk_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_compare_pdb(threads):
    """test phold compare with pdbs input"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_pdb_dir} -t {threads} -d {database_dir} --structures --structure_dir {pdb_dir} -f"
    exec_command(cmd)

def test_compare_pdb_custom(threads):
    """test phold compare with pdbs input with custom db"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_pdb_dir_custom} -t {threads} -d {database_dir} --structures --structure_dir {pdb_dir} --custom_db {dummy_custom_db}  -f"
    exec_command(cmd)

def test_compare_cif(threads):
    """test phold compare with AF3 cif input"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold compare -i {input_gbk} -o {compare_cif_dir} -t {threads} -d {database_dir} --structures --structure_dir {cif_dir} -f"
    exec_command(cmd)

def test_proteins_compare_pdb(threads):
    """test phold proteins-compare with cifs input"""
    input_faa: Path = f"{test_data}/NC_043029_aa.fasta"
    cmd = f"phold proteins-compare -i {input_faa} -o {proteins_compare_pdb_dir} -t {threads} -d {database_dir} --structures --structure_dir {cif_dir}  -f"
    exec_command(cmd)

def test_proteins_compare_cifs_custom_db(threads):
    """test phold proteins-compare with cifs input"""
    input_faa: Path = f"{test_data}/NC_043029_aa.fasta"
    cmd = f"phold proteins-compare -i {input_faa} -o {proteins_compare_cif_dir_custom} -t {threads} -d {database_dir} --structures --structure_dir {cif_dir} --custom_db {dummy_custom_db}  -f"
    exec_command(cmd)

def test_proteins_compare_cif(threads):
    """test phold proteins-compare with AF3 cif input"""
    input_faa: Path = f"{test_data}/NC_043029_aa.fasta"
    cmd = f"phold proteins-compare -i {input_faa} -o {proteins_compare_cif_dir} -t {threads} -d {database_dir} --structures --structure_dir {cif_dir}  -f"
    exec_command(cmd)


def test_predict_fasta(gpu_available, threads):
    """test phold predict with fasta input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
    cmd = f"phold predict -i {input_fasta} -o {predict_fasta_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)


def test_compare_fasta(threads):
    """test phold compare with fasta input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
    cmd = f"phold compare -i {input_fasta} -o {compare_fasta_dir} --predictions_dir {predict_fasta_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_proteins_predict(gpu_available, threads):
    """test phold proteins-predict"""
    input_fasta: Path = f"{test_data}/phanotate.faa"
    cmd = f"phold proteins-predict -i {input_fasta} -o {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_proteins_predict_gzip(gpu_available, threads):
    """test phold proteins-predict with gzip"""
    input_fasta: Path = f"{test_data}/phanotate.faa.gz"
    cmd = f"phold proteins-predict -i {input_fasta} -o {proteins_predict_dir} -t {threads} -d {database_dir} -f"
    if gpu_available is False:
        cmd = f"{cmd} --cpu"
    exec_command(cmd)

def test_proteins_compare(threads):
    """test phold proteins-compare"""
    input_fasta: Path = f"{test_data}/phanotate.faa"
    cmd = f"phold proteins-compare -i {input_fasta} --predictions_dir {proteins_predict_dir} -o {proteins_compare_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_proteins_compare_gzip(threads):
    """test phold proteins-compare"""
    input_fasta: Path = f"{test_data}/phanotate.faa.gz"
    cmd = f"phold proteins-compare -i {input_fasta} --predictions_dir {proteins_predict_dir} -o {proteins_compare_dir} -t {threads} -d {database_dir} -f"
    exec_command(cmd)


def test_plot():
    """test phold plot"""
    input_gbk: Path = f"{test_data}/NC_043029.gbk"
    cmd = f"phold plot -i {input_gbk}  -o {plots_dir} -f"
    exec_command(cmd)


def test_remote_genbank(run_remote, threads):
    """test phold remote with genbank input"""
    input_gbk: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.gbk"
    if run_remote is True:
        cmd = f"phold remote -i {input_gbk} -o {remote_gbk_dir} -t {threads} -d {database_dir} -f"
        exec_command(cmd)


def test_remote_fasta(run_remote, threads):
    """test phold remote with fasta input"""
    input_fasta: Path = f"{test_data}/combined_truncated_acr_defense_vfdb_card.fasta"
    if run_remote is True:
        cmd = f"phold remote -i {input_fasta} -o {remote_fasta_dir} -t {threads} -d {database_dir} -f"
        exec_command(cmd)


# class testFails(unittest.TestCase):
#     """Tests for fails"""
   
#     def test_run_hyps_ncbi(gpu_available, threads):
#         """test phold run with --hyps but not pharokka input"""
#         input_gbk: Path = f"{test_data}/NC_043029_ncbi.gbk"
#         cmd = f"phold run -i {input_gbk} -o {run_gbk_hyps} -t {threads} -d {database_dir} -f --hpys"
#         if gpu_available is False:
#             cmd = f"{cmd} --cpu"
#         exec_command(cmd)


remove_directory(output_dir)
