from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from loguru import logger


def open_protein_fasta_as_cds_dict(input_fasta):
    """
    Parse a (possibly gzipped) protein FASTA file into a CDS-style dictionary.

    All sequences are converted to uppercase to ensure downstream compatibility
    (e.g. 3Di / ProstT5 models).

    Parameters
    ----------
    input_fasta : str
        Path to protein FASTA file (can be gzipped).
    logger : logging.Logger
        Logger for warnings and errors.

    Returns
    -------
    dict
        Dictionary with structure:
        {
            "proteins": {
                protein_id: SeqFeature(...)
            }
        }
    """
    cds_dict = {"proteins": {}}

    with open_protein_fasta_file(input_fasta) as handle:  # handles gzip too
        records = list(SeqIO.parse(handle, "fasta"))

    if not records:
        logger.warning(f"No proteins were found in your input file {input_fasta}.")
        logger.error(
            f"Your input file {input_fasta} is likely not an amino acid FASTA file. "
            f"Please check this."
        )
        return cds_dict

    for record in records:
        prot_id = record.id
        seq_str = str(record.seq)

        if any(c.islower() for c in seq_str):
            logger.warning("Lower case amino acid input detected")
            logger.warning(
                "This will automatically be converted to uppercase or else "
                "the 3Di prediction will not work"
            )

        feature_location = FeatureLocation(0, len(seq_str))

        seq_feature = SeqFeature(
            feature_location,
            type="CDS",
            qualifiers={
                "ID": record.id,
                "description": record.description,
                "translation": seq_str.upper(),
            },
        )

        cds_dict["proteins"][prot_id] = seq_feature

    if not cds_dict["proteins"]:
        logger.error(f"Error: no AA protein sequences found in {input_fasta}")

    return cds_dict
