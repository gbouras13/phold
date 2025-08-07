import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import requests
from alive_progress import alive_bar
from loguru import logger

from phold.utils.util import remove_directory
from phold.utils.external_tools import ExternalTool

# set this if changes
CURRENT_DB_VERSION: str = "1.0.0"

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
    },

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
    },
    "1.0.0": {
        "md5": "ddbe0d94b1d94a392cfeb4ec113f4362",
        "major": 1,
        "minor": 0,
        "minorest": 0, 
        "db_url": "https://zenodo.org/records/16741548/files/phold_search_db_v_1_0_0.tar.gz",
        "dir_name": "phold_search_db_v_1_0_0",
        "tarball": "phold_search_db_v_1_0_0.tar.gz",
        "prostt5_backup_url": "https://zenodo.org/records/11234657/files/models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_tarball": "models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_md5": "118c1997e6d2cb5025abda95d36681e0",
    }
}

# for the extended DB

VERSION_DICTIONARY_3M16 = {
    "1.0.0": {
        "md5": "10aaafa55b5a28c04c08833f3c787097",
        "major": 1,
        "minor": 0,
        "minorest": 0,
        "db_url": "https://zenodo.org/records/16741548/files/phold_db_3M16_v_1_0_0.tar.gz",
        "dir_name": "phold_db_3M16_v_1_0_0",
        "tarball": "phold_db_3M16_v_1_0_0.tar.gz",
        "prostt5_backup_url": "https://zenodo.org/records/11234657/files/models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_tarball": "models--Rostlab--ProstT5_fp16.tar.gz",
        "prostt5_backup_md5": "118c1997e6d2cb5025abda95d36681e0",
    }
}



PHOLD_DB_NAMES = [
    "acrs_plddt_over_70_metadata.tsv",
    "all_phold_structures",
    "all_phold_structures_ca",
    "all_phold_structures_ca.dbtype",
    "all_phold_structures_ca.index",
    "all_phold_structures.dbtype",
    "all_phold_structures_h",
    "all_phold_structures_h.index",
    "all_phold_structures_h.dbtype",
    "all_phold_structures.index",
    "all_phold_structures.lookup",
    "all_phold_structures.source",
    "all_phold_structures_ss",
    "all_phold_structures_ss.dbtype",
    "all_phold_structures_ss.index",
    "card_plddt_over_70_metadata.tsv",
    "defensefinder_plddt_over_70_metadata.tsv",
    "netflax_annotation_table.tsv",
    "vfdb_description_output.csv",
    "phold_annots.tsv"
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

PROSTT5_FINETUNE_MD5_DICTIONARY = {
    "refs": {"main": "833198a97a993b079841bb51b3096b6a"},
    "blobs": {
        "17ade346a1042cbe0c1436f5bedcbd85c099d582": "9583322ae544dfbdb0001350ab4c3fd9",
        "474673ae9f1fcf1da9b939df0da5888abf0674af6e50546927b4edd4536b8897": "028128ddd5340a8bfb4151f562cb2091",
        "47ce84c8a086256c31ca8aa555e55787f9368b3e": "865fce8e8f6db4ad9bbb0f9596438e1e",
        "74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324": "6ad28d9980aaec37d3935072d204e520",
        "7b4566bb3be907e8bba6c764ecd41d3388b81b79": "67183db6031983c323e7fbb0ae7dac57",
        "a45ead5f5f3e3b842f660aa8e7ab2d66cbb959ad743f8677d0559c13f7f4a0c4": "b7cb2d32552a3d3cabf4db25e54fdb44",
        "c72eadd6ad03b26932928906cd12a29586b9fd7217775cd8207fa38cd52ebaaf" : "5cacb3a7b9789143e6a28279a25a0bad",
        "d91b97e1ea1225fcbc4e5ee14c09bb196251cee9" : "ce3c5c2e2e583314c72b3a5c823118f8"
    },
}


PHOLD_DB_FOLDSEEK_GPU_NAMES = [
    "all_phold_structures_gpu"
]


def install_database(db_dir: Path, foldseek_gpu: bool, extended_db: bool, threads: int) -> None:
    """
    Install the Phold database.

    Args:
        db_dir Path: The directory where the database should be installed.
        foldseek_gpu bool: Whether to install foldseek-gpu compatible phold db
        extended_db bool: Whether to download the extended Phold DB 3.16M with 1.8M unknown function efam and enVhog proteins
        threads int: Number of threads available (makes downloading faster)
    """

    # check the database is installed
    logger.info(f"Checking Phold database installation in {db_dir}.")
    downloaded_flag, gpu_flag = check_db_installation(db_dir, foldseek_gpu)
    if downloaded_flag:
        logger.info("All Phold databases files are present")
    else:
        logger.info("Some Phold databases files are missing")

        if extended_db:
            DICT = VERSION_DICTIONARY_3M16
            db_url = DICT[CURRENT_DB_VERSION]["db_url"]
            logger.info(f"Downloading Phold DB 3.16M from {db_url}")

        else:
            DICT = VERSION_DICTIONARY
            db_url = DICT[CURRENT_DB_VERSION]["db_url"]
            logger.info(f"Downloading Phold Search DB 1.36M from {db_url}")

        requiredmd5 = DICT[CURRENT_DB_VERSION]["md5"]
        tarball = DICT[CURRENT_DB_VERSION]["tarball"]

        
        tarball_path = Path(f"{db_dir}/{tarball}")
        logdir = Path(db_dir) / "logdir"

        download(db_url, tarball_path, logdir, threads)

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
        untar(tarball_path, db_dir, DICT)
        tarball_path.unlink()

    if foldseek_gpu:
        if gpu_flag:
            logger.info("All Phold database files compatible with Foldseek-GPU are present")
        else:
            logger.info("Some Phold database files compatible with Foldseek-GPU are missing")
            logger.info("Creating them")
            foldseek_makepaddedseqdb(db_dir)



"""
lots of this code from the marvellous bakta https://github.com/oschwengers/bakta, db.py specifically
"""


# def download(db_url: str, tarball_path: Path) -> None:
#     """
#     Download the database from the given URL.

#     Args:
#         db_url (str): The URL of the database.
#         tarball_path (Path): The path where the downloaded tarball should be saved.
#     """
#     try:
#         with tarball_path.open("wb") as fh_out, requests.get(
#             db_url, stream=True
#         ) as resp:
#             total_length = resp.headers.get("content-length")
#             if total_length is not None:  # content length header is set
#                 total_length = int(total_length)
#             with alive_bar(total=total_length, scale="SI") as bar:
#                 for data in resp.iter_content(chunk_size=1024 * 1024):
#                     fh_out.write(data)
#                     bar(count=len(data))
#     except IOError:
#         logger.error(
#             f"ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}"
#         )

"""
aria2c bottlenecked by Zenodo but still faster
"""

def download(db_url: str, tarball_path: Path, logdir: Path, threads: int) -> None:
    """
    Download the database from the given URL using aria2c.

    Args:
        db_url (str): The URL of the database.
        tarball_path (Path): The path where the downloaded tarball should be saved.
        logdir (Path): The path to store logs
        threads (int): Number of threads for aria2c
    """

    cmd = f"--dir {str(tarball_path.parent)} --out {tarball_path.name} --max-connection-per-server={str(threads)} --allow-overwrite=true  {db_url}"

    download_db = ExternalTool(
        tool="aria2c",
        input=f"",
        output=f"",
        params=f"{cmd}",
        logdir=logdir,
    )
    try:
        ExternalTool.run_download(download_db)
    except:
        logger.warning("Downloading the database with aria2c failed. Trying now without.")
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






#     try:
#         with tarball_path.open("wb") as fh_out, requests.get(
#             db_url, stream=True
#         ) as resp:
#             total_length = resp.headers.get("content-length")
#             if total_length is not None:  # content length header is set
#                 total_length = int(total_length)
#             with alive_bar(total=total_length, scale="SI") as bar:
#                 for data in resp.iter_content(chunk_size=1024 * 1024):
#                     fh_out.write(data)
#                     bar(count=len(data))
#     except IOError:
#         logger.error(
#             f"ERROR: Could not download file from Zenodo! url={db_url}, path={tarball_path}"
#         )


def download_zenodo_prostT5(model_dir, logdir, threads):
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

    download(db_url, tarball_path, logdir, threads)
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

    if model_name == "Rostlab/ProstT5_fp16":

        model_sub_dir = "models--Rostlab--ProstT5_fp16"
        DICT = PROSTT5_MD5_DICTIONARY

    elif model_name == "gbouras13/ProstT5Phold":

        model_sub_dir = "models--gbouras13--ProstT5Phold"
        DICT = PROSTT5_FINETUNE_MD5_DICTIONARY


    for key in DICT:
        for nested_key in DICT[key]:
            file_path = Path(
                f"{model_dir}/{model_sub_dir}/{key}/{nested_key}"
            )

            # check file exists
            if file_path.exists():
                md5_sum = calc_md5_sum(file_path)
                if md5_sum != DICT[key][nested_key]:
                    logger.warning(
                        f"Corrupt model file {file_path}! MD5 should be '{DICT[key][nested_key]}' but is '{md5_sum}'"
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


def untar(tarball_path: Path, output_path: Path, DICT: dict) -> None:
    """
    Extract the tarball to the output path.

    Args:
        tarball_path (Path): The path to the tarball file.
        output_path (Path): The path where the contents of the tarball should be extracted.
        DICT (dict): version dictionary
    """
    try:
        with tarball_path.open("rb") as fh_in, tarfile.open(
            fileobj=fh_in, mode="r:gz"
        ) as tar_file:
            tar_file.extractall(path=str(output_path))

        tarpath = Path(output_path) / DICT[CURRENT_DB_VERSION]["dir_name"]

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


def check_db_installation(db_dir: Path, foldseek_gpu: bool) -> bool:
    """
    Check if the Phold database is installed.

    Args:
        db_dir Path: The directory where the database is installed.
        foldseek_gpu bool: Whether to install foldseek-gpu compatible phold db

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
    
    gpu_flag = True
    if foldseek_gpu:
        for file_name in PHOLD_DB_FOLDSEEK_GPU_NAMES:
            path = Path(db_dir) / file_name
            if not path.is_file():
                logger.warning(f"Phold Foldseek-GPU Database file {path} is missing")
                gpu_flag = False
                break 

    return downloaded_flag, gpu_flag


def validate_db(database: str, default_dir: str, foldseek_gpu: bool) -> Path:
    """
    Validates the Phold database is installed.

    Args:
        database str: The directory where the database is installed.
        default_dir str: Default DB location
        foldseek_gpu bool: Whether to install foldseek-gpu compatible phold db

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
    downloaded_flag, gpu_flag = check_db_installation(database, foldseek_gpu)
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
    if foldseek_gpu:
        if gpu_flag:
            logger.info("All Phold database files compatible with Foldseek-GPU are present")
        else:
            logger.error(
                f"Phold database files compatible with Foldseek-GPU not found. Please run phold install -d {database} --foldseek_gpu"
            )


    return database

# for now, until we have a better confidence metric, use PyTorch
# def foldseek_gpu_prostt5_download(db_dir: Path) -> None:

#     prostt5_db_path = Path(db_dir) / "prostt5_weights"
#     tmp_dir = Path(db_dir) / "tmp"
#     logdir = Path(db_dir) / "logdir"

#     foldseek_createdb_gpu = ExternalTool(
#         tool="foldseek",
#         input=f"",
#         output=f"",
#         params=f"databases ProstT5 {prostt5_db_path} {tmp_dir}  ",
#         logdir=logdir,
#     )

#     ExternalTool.run_tool(foldseek_createdb_gpu)
#     remove_directory(tmp_dir)

def foldseek_makepaddedseqdb(db_dir: Path) -> None:

    phold_db_search = Path(db_dir) / "all_phold_structures"
    phold_db_search_gpu = Path(db_dir) / "all_phold_structures_gpu"
    logdir = Path(db_dir) / "logdir"

    foldseek_makepaddedseqdb = ExternalTool(
        tool="foldseek",
        input=f"",
        output=f"",
        params=f"makepaddedseqdb {phold_db_search} {phold_db_search_gpu}",
        logdir=logdir,
    )

    ExternalTool.run_tool(foldseek_makepaddedseqdb)



    
