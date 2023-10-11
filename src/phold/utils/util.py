"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import os
import subprocess as sp
import sys
import time

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


def begin_phold(params):
    """
    begins phold
    params: dict params a dictionary of params for hazed
    returns: int start time
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(params["--output"], f"phold_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")

    print_splash()
    logger.info(
        "phold: annotating phage genomes with protein structures."
    )

    logger.info(f"You are using phold version {get_version()}")
    logger.info("Repository homepage is https://github.com/gbouras13/phold")
    logger.info(f"Listing parameters.")
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}.")

    return start_time



def end_phold(start_time):
    """
    finishes phold
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("phold has finished")
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