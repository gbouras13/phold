#!/usr/bin/env python3
"""phold"""

import os
import shutil
from pathlib import Path

import click
from loguru import logger

from phold.io.handle_genbank import get_genbank

# from phold.utils.all import all_process_blast_output_and_reorient
# from phold.utils.bulk import bulk_process_blast_output_and_reorient, run_bulk_blast
# from phold.utils.cds_methods import run_blast_based_method, run_mystery, run_nearest
# from phold.utils.constants import phold_DB
# from phold.utils.external_tools import ExternalTool
from phold.utils.util import (
    begin_phold,
    end_phold,
    get_version,
    print_citation
)

from phold.utils.validation import instantiate_dirs
# from phold.utils.validation import (
#     check_evalue,
#     instantiate_dirs,
#     validate_choice_autocomplete,
#     validate_choice_mode,
#     validate_custom_db_fasta,
#     validate_fasta,
#     validate_fasta_all,
#     validate_fasta_bulk,
#     validate_ignore_file,
# )


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-i",
            "--input",
            help="Path to input file in Genbank format",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "-o",
            "--output",
            default="output_phold",
            show_default=True,
            type=click.Path(),
            help="Output directory ",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads to use with Foldseek",
            default=1,
            show_default=True,
        ),
        click.option(
            "-p",
            "--prefix",
            default="phold",
            help="Prefix for output files",
            show_default=True,
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output directory",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
Chromosome command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-5",
    help="e value threshold for Foldseek",
    show_default=True,
)
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    **kwargs,
):
    """Runs phold"""

    # validates the directory  (need to before I start phold or else no log file is written)
    instantiate_dirs(output, force)

    output: Path = Path(output)
    logdir: Path = Path(output) / "logs"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue
    }


    # initial logging etc
    start_time = begin_phold(params)

    # validates fasta
    gb_dict = get_genbank(input)
    if not gb_dict:
        logger.warning("Error: no sequences found in genbank file")
        logger.error("No sequences found in genbank file. Nothing to annotate")

    for key, value in gb_dict.items():
        logger.info(f"Parameter: {key} {value}.")

    # Create a nested dictionary to store CDS features by contig ID
    cds_dict = {}

    for record_id, record in gb_dict.items():
        
        # set level 1
        cds_dict[record_id] = {}

        #print(record.features)
        for cds_feature in record.features:
            if cds_feature.type == 'CDS':
                print(cds_feature.qualifiers['ID'][0])
                cds_dict[record_id][cds_feature.qualifiers['ID'][0]] = cds_feature

    print(cds_dict)

        #         cds_sequence = cds_feature.extract(record.seq)
        #         cds_sequences.append(str(cds_sequence))

        # Add the CDS sequences to the nested dictionary under the contig ID
        # cds_dict[record_id] = cds_sequences

    # position = [
    #     (int(this_CDS[i].location.start), int(this_CDS[i].location.end))
    #     for i in range(len(this_CDS))
    # ]
    # sense = [
    #     re.split("]", str(this_CDS[i].location))[1][1] for i in range(len(this_CDS))
    # ]
    # protein_id = [
    #     this_CDS[i].qualifiers.get("protein_id") for i in range(len(this_CDS))
    # ]
    # protein_id = [p[0] if p is not None else None for p in protein_id]
    # phrogs = [this_CDS[i].qualifiers.get("phrog") for i in range(len(this_CDS))]
    # phrogs = ["No_PHROG" if i is None else i[0] for i in phrogs]

    # return {
    #     "length": phage_length,
    #     "phrogs": phrogs,
    #     "protein_id": protein_id,
    #     "sense": sense,
    #     "position": position,
    # }

    #print(gb_dict)


    # validates fasta
    #validate_fasta(input)

    # validate e value
    #check_evalue(evalue)



    # end phold
    end_phold(start_time)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


# main_cli.add_command(run)
main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
