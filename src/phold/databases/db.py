"""
to tar DBs

# v0.1.0
GZIP=-9 tar cvzf phold_structure_foldseek_db.tar.gz phold_structure_foldseek_db

# v0.2.0
GZIP=-9 tar cvzf phold_db_v_0_2_0.tar.gz phold_db_v_0_2_0

"""

import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger

from phold.utils.util import remove_directory

# set this if changes
CURRENT_DB_VERSION: str = "0.2.0"

# to hold information about the different DBs
VERSION_DICTIONARY = {
    "0.1.0": {
        "md5": "353a1a6763e1261c5c44e1e2da9d8736",
        "major": 0,
        "minor": 1,
        "minorest": 0,
        "db_url": "https://zenodo.org/records/10675285/files/phold_structure_foldseek_db.tar.gz",
        "dir_name": "phold_structure_foldseek_db",
        "tarball": "phold_structure_foldseek_db.tar.gz",
        "prostt5_backup_url": "https://zenodo.org/records/11234657/files/models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_tarball": "models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_md5": "118c1997e6d2cb5025abda95d36681e0",
    }
}

VERSION_DICTIONARY = {
    "0.2.0": {
        "md5": "99ed8b4bcc41ca6e05e8690ba7e85197",
        "major": 0,
        "minor": 2,
        "minorest": 0,
        "db_url": "https://zenodo.org/records/12735568/files/phold_db_v_0_2_0.tar.gz",
        "dir_name": "phold_db_v_0_2_0",
        "tarball": "phold_db_v_0_2_0.tar.gz",
        "prostt5_backup_url": "https://zenodo.org/records/11234657/files/models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_tarball": "models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_md5": "118c1997e6d2cb5025abda95d36681e0",
    }
}

CURRENT_VERSION = "0.2.0"

PHOLD_DB_NAMES = [
    "acrs_plddt_over_70_metadata.tsv",
    "all_phold_structures",
    "all_phold_structures_clustered_searchDB",
    "all_phold_structures_clustered_searchDB_ca",
    "all_phold_structures_clustered_searchDB_ca.dbtype",
    "all_phold_structures_clustered_searchDB_ca.index",
    "all_phold_structures_clustered_searchDB_clu",
    "all_phold_structures_clustered_searchDB_clu.dbtype",
    "all_phold_structures_clustered_searchDB_clu.index",
    "all_phold_structures_clustered_searchDB.dbtype",
    "all_phold_structures_clustered_searchDB_h",
    "all_phold_structures_clustered_searchDB_h.index",
    "all_phold_structures_clustered_searchDB_h.dbtype",
    "all_phold_structures_clustered_searchDB.index",
    "all_phold_structures_clustered_searchDB.lookup",
    "all_phold_structures_clustered_searchDB_seq.1",
    "all_phold_structures_clustered_searchDB_seq_ca.1",
    "all_phold_structures_clustered_searchDB_seq_ca.dbtype",
    "all_phold_structures_clustered_searchDB_seq_ca.index",
    "all_phold_structures_clustered_searchDB_seq_h.1",
    "all_phold_structures_clustered_searchDB_seq_h.dbtype",
    "all_phold_structures_clustered_searchDB_seq_h.index",
    "all_phold_structures_clustered_searchDB_seq.index",
    "all_phold_structures_clustered_searchDB_seq_ss.1",
    "all_phold_structures_clustered_searchDB_seq_ss.dbtype",
    "all_phold_structures_clustered_searchDB_seq_ss.index",
    "all_phold_structures_clustered_searchDB.source",
    "all_phold_structures_clustered_searchDB_ss",
    "all_phold_structures_clustered_searchDB_ss.dbtype",
    "all_phold_structures_clustered_searchDB_ss.index",
    "all_phold_structures.dbtype",
    "all_phold_structures_h",
    "all_phold_structures_h.index",
    "all_phold_structures_h.dbtype",
    "all_phold_structures.index",
    "all_phold_structures.lookup",
    "all_phold_structures.source",
    "card_plddt_over_70_metadata.tsv",
    "vfdb_description_output.csv",
    "netflax_annotation_table.tsv",
    "phold_annots.tsv",
    "defensefinder_plddt_over_70_metadata.tsv",
]


PROSTT5_MD5_DICTIONARY = {
    "refs": {"main": "962133e8e2bff04ec1768fa58dd788f3"},
    "blobs": {
        "2c19eb6e3b583f52d34b903b5978d3d30b6b7682": "8fd03e945174de0818746ecbde1aad8e",
        "60fe6bb247c90b8545d7b73820cd796ce6dcbd59": "ce32377619690072ced051ec7fc6c5f0",
        "6fc7be92c58e238f20a6cdea5a87b123a4ad35e2": "1deb27443d0d9b90b8e5704e752373e2",
        "74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324": "6ad28d9980aaec37d3935072d204e520",
        "b1a9ffcef73280cc57f090ad6446b4116b574b6c75d83ccc32778282f7f00855": "45f066fc3b0d87fb9f98bb0ddb97a3dc",
        "e9322396e6e75ecf8da41a9527e24dfa4eeea505": "b1cdd31ea50a37bf84cc0d7ef11820c8",
    },
}


def install_database(db_dir: Path) -> None:
    """
    Install the Phold database.

    Args:
        db_dir Path: The directory where the database should be installed.
    """

    # check the database is installed
    logger.info(f"Checking Phold database installation in {db_dir}.")
    downloaded_flag = check_db_installation(db_dir)
    if downloaded_flag == True:
        logger.info("All Phold databases files are present")
    else:
        logger.info("Some Phold databases files are missing")
        logger.info("Downloading the Phold database")

        db_url = VERSION_DICTIONARY[CURRENT_DB_VERSION]["db_url"]
        requiredmd5 = VERSION_DICTIONARY[CURRENT_DB_VERSION]["md5"]

        logger.info(f"Downloading Phold database from {db_url}")

        tarball = VERSION_DICTIONARY[CURRENT_DB_VERSION]["tarball"]
        tarball_path = Path(f"{db_dir}/{tarball}")

        download(db_url, tarball_path)

        md5_sum = calc_md5_sum(tarball_path)

        if md5_sum == requiredmd5:
            logger.info(f"Phold database file download OK: {md5_sum}")
        else:
            logger.error(
                f"Error: corrupt database file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
            )

        logger.info(
            f"Extracting Phold database tarball: file={tarball_path}, output={db_dir}"
        )
        untar(tarball_path, db_dir)
        tarball_path.unlink()


"""
lots of this code from the marvellous bakta https://github.com/oschwengers/bakta, db.py specifically
"""


def download(db_url: str, tarball_path: Path) -> None:
    """
    Download the database from the given URL.

    Args:
        db_url (str): The URL of the database.
        tarball_path (Path): The path where the downloaded tarball should be saved.
    """
    try:
        with tarball_path.open("wb") as fh_out, requests.get(
            db_url, stream=True
        ) as resp:
            total_length = resp.headers.get("content-length")
            if total_length is not None:  # content length header is set
                total_length = int(total_length)
            with alive_bar(total=total_length, scale="SI") as bar:
                for data in resp.iter_content(chunk_size=1024 * 1024):
                    fh_out.write(data)
                    bar(count=len(data))
    except IOError:
        logger.error(
            f"ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}"
        )


def download_zenodo_prostT5(model_dir):
    """
    Download the ProstT5 model from Zenodo

    Args:
        db_url (str): The URL of the database.
        tarball_path (Path): The path where the downloaded tarball should be saved.
    """

    db_url = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_url"]
    requiredmd5 = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_md5"]

    logger.info(f"Downloading ProstT5 model backup from {db_url}")

    tarball = VERSION_DICTIONARY[CURRENT_DB_VERSION]["prostt5_backup_tarball"]
    tarball_path = Path(f"{model_dir}/{tarball}")

    download(db_url, tarball_path)
    md5_sum = calc_md5_sum(tarball_path)

    if md5_sum == requiredmd5:
        logger.info(f"ProstT5 model backup file download OK: {md5_sum}")
    else:
        logger.error(
            f"Error: corrupt file! MD5 should be '{requiredmd5}' but is '{md5_sum}'"
        )

    logger.info(
        f"Extracting ProstT5 model backup tarball: file={tarball_path}, output={model_dir}"
    )

    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(model_dir))

    except OSError:
        logger.warning("Encountered OSError: {}".format(OSError))
        logger.error(f"Could not extract {tarball_path} to {model_dir}")

    tarball_path.unlink()


def check_prostT5_download(model_dir: Path, model_name: str) -> bool:
    """
     Args:
        model_dir (Path): Directory where the model and tokenizer is be stored.
        model_name (str): Name of the pre-trained T5 model.
    Returns:
        bool: bool to tell Phold whether to download ProstT5
    """

    # assumes already has been downloaded
    download = False

    for key in PROSTT5_MD5_DICTIONARY:
        for nested_key in PROSTT5_MD5_DICTIONARY[key]:

            file_path = Path(
                f"{model_dir}/models--Rostlab--ProstT5_fp16/{key}/{nested_key}"
            )

            # check file exists
            if file_path.exists():
                md5_sum = calc_md5_sum(file_path)
                if md5_sum != PROSTT5_MD5_DICTIONARY[key][nested_key]:
                    logger.warning(
                        f"Corrupt model file {file_path}! MD5 should be '{PROSTT5_MD5_DICTIONARY[key][nested_key]}' but is '{md5_sum}'"
                    )
                    download = True
            else:
                logger.warning(f"Model file {file_path} does not exist.")
                download = True

    return download


def calc_md5_sum(tarball_path: Path, buffer_size: int = 1024 * 1024) -> str:
    """
    Calculate the MD5 checksum of the given file.

    Args:
        tarball_path (Path): The path to the file for which the MD5 checksum needs to be calculated.
        buffer_size (int): The buffer size for reading the file.

    Returns:
        str: The MD5 checksum of the file.
    """

    md5 = hashlib.md5()
    with tarball_path.open("rb") as fh:
        data = fh.read(buffer_size)
        while data:
            md5.update(data)
            data = fh.read(buffer_size)
    return md5.hexdigest()


def untar(tarball_path: Path, output_path: Path) -> None:
    """
    Extract the tarball to the output path.

    Args:
        tarball_path (Path): The path to the tarball file.
        output_path (Path): The path where the contents of the tarball should be extracted.
    """
    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(output_path))

        tarpath = Path(output_path) / VERSION_DICTIONARY[CURRENT_DB_VERSION]["dir_name"]

        # Get a list of all files in the directory
        files_to_move = [f for f in tarpath.iterdir() if f.is_file()]

        # Move each file to the destination directory
        for file_name in files_to_move:
            destination_path = output_path / file_name.name
            shutil.move(file_name, destination_path)
        # remove the directory
        remove_directory(tarpath)

    except OSError:
        logger.warning("Encountered OSError: {}".format(OSError))
        logger.error(f"Could not extract {tarball_path} to {output_path}")


def check_db_installation(db_dir: Path) -> bool:
    """
    Check if the Phold database is installed.

    Args:
        db_dir Path: The directory where the database is installed.

    Returns:
        bool: True if all required files are present, False otherwise.
    """
    downloaded_flag = True
    for file_name in PHOLD_DB_NAMES:
        path = Path(db_dir) / file_name
        if not path.is_file():
            logger.warning(f"Phold Database file {path} is missing")
            downloaded_flag = False
            break

    return downloaded_flag


def validate_db(database: str, default_dir: str) -> Path:
    """
    Validates the Phold database is installed.

    Args:
        database str: The directory where the database is installed.
        default_dir str: Default DB location

    Returns:
        bool: True if all required files are present, False otherwise.
    """
    # set default DB if not specified
    if database is not None:
        database: Path = Path(database)
    else:
        database = Path(default_dir)

    # check the database is installed
    logger.info(f"Checking Phold database installation in {database}")
    downloaded_flag = check_db_installation(database)
    if downloaded_flag == True:
        logger.info("All Phold databases files are present")
    else:
        if database == Path(default_dir):  # default
            logger.error(
                f"Phold database not found. Please run phold install to download and install the Phold database"
            )
        else:  # specific
            logger.error(
                f"Phold database not found. Please run phold install -d {database} to download and install the Phold database"
            )

    return database
