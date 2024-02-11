import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path

import click
from loguru import logger


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def phold_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(phold_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    with open(phold_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
begin and end functions
"""


def begin_phold(params, subcommand: str):
    """
    begins phold
    params: dict params a dictionary of params for hazed
    returns: int start time
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(params["--output"], f"phold_{subcommand}_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")

    print_splash()
    logger.info("phold: annotating phage genomes with protein structures")

    logger.info(f"You are using phold version {get_version()}")
    logger.info("Repository homepage is https://github.com/gbouras13/phold")
    logger.info(f"You are running phold {subcommand}")
    logger.info(f"Listing parameters")
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}")

    return start_time


def end_phold(start_time, subcommand: str):
    """
    finishes phold
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info(f"phold {subcommand} has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")


# need the logo here eventually
def print_splash():
    click.echo(
        """\b

.______    __    __    ______    __       _______  
|   _  \  |  |  |  |  /  __  \  |  |     |       \ 
|  |_)  | |  |__|  | |  |  |  | |  |     |  .--.  |
|   ___/  |   __   | |  |  |  | |  |     |  |  |  |
|  |      |  |  |  | |  `--'  | |  `----.|  '--'  |
| _|      |__|  |__|  \______/  |_______||_______/ 
                                                   

"""
    )


def remove_file(file_path: Path) -> None:
    if file_path.exists():
        file_path.unlink()  # Use unlink to remove the file


def remove_directory(dir_path: Path) -> None:
    if dir_path.exists():
        shutil.rmtree(dir_path)


def touch_file(path: Path) -> None:
    with open(path, "a"):
        os.utime(path, None)


def clean_up_temporary_files(output: Path):
    result_high_tsv: Path = Path(output) / "foldseek_results_high.tsv"
    result_low_tsv: Path = Path(output) / "foldseek_results_low.tsv"
    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    result_db_base: Path = Path(output) / "result_db"
    temp_db: Path = Path(output) / "temp_db"
    remove_directory(result_db_base)
    remove_directory(temp_db)
    remove_file(result_tsv)
    remove_file(result_high_tsv)
    remove_file(result_low_tsv)
