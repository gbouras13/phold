from pathlib import Path
from loguru import logger
import re
from tqdm import tqdm
import numpy as np
from pathlib import Path
import os
import shutil
from phold.features.create_foldseek_db import  foldseek_tsv2db
from phold.utils.util import remove_file


def generate_mmseqs_db_from_aa(
    cds_dict: Path,  output: Path, logdir: Path, prefix: str, proteins_flag: bool
) -> None:
    """
    Generate MMSeqs2 database from amino-acid sequences - for use with profiles later

    Args:
        fasta_aa (Path): Path to the amino-acid FASTA file.
        output (Path): Path to output directory.
        logdir (Path): Path to the directory where logs will be stored.
        prefix (str): Prefix for the Foldseek database.
        proteins_flag (bool): True if phold proteins-predict

    Returns:
        None
    """

    lookup = {}

    mmseqs2_db_dir: Path = Path(output) / f"query_profiledb" # this will be subdir where the mmseqs2 db is 
    mmseqs2_db_dir.mkdir(parents=True, exist_ok=True)

    # create MMSeqs2 db names
    short_db_name = f"{prefix}"

    lookup_db_name: Path = Path(mmseqs2_db_dir) / f"{short_db_name}.lookup"

    temp_aa_tsv = Path(output) / "aa.tsv"
    temp_header_tsv = Path(output) / "header.tsv"

    # with open(temp_aa_tsv, "w") as aa_f, open(temp_header_tsv, "w") as h_f, open(lookup_db_name, "w") as l_f:
    with open(temp_aa_tsv, "w") as aa_f, open(temp_header_tsv, "w") as h_f:
        idx = 1

        for contig_id, aa_contig_dict in cds_dict.items():
            for seq_id, feature in aa_contig_dict.items():
                aa_f.write(f"{idx}\t{feature.qualifiers["translation"]}\n")
                h_f.write(f"{idx}\t{seq_id}\n")
                lookup[seq_id] = f"{idx}"
                #l_f.write(f"{idx-1}\t{seq_id}\t0") # I think I dont need this
                idx += 1

    
    aa_db_name: Path = Path(mmseqs2_db_dir) / short_db_name
    header_db_name: Path = Path(mmseqs2_db_dir) / f"{short_db_name}_h"

    # create Foldseek database with foldseek tsv2db

    foldseek_tsv2db(temp_aa_tsv, aa_db_name, 0, logdir)
    foldseek_tsv2db(temp_header_tsv, header_db_name, 12, logdir)

    # clean up
    remove_file(temp_aa_tsv)
    remove_file(temp_header_tsv)

    # return the lookup dict no need to actually write
    return lookup 


def build_lookup(filename):
    """
    builds foldseek profile lookup using mmseqs lookup
    """
    lookup = {}
    with open(filename, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                key, value = parts[1], parts[0]
                lookup[key] = value
    return lookup

#############################
# pssm_to_db functionality
#############################

blosum62_str = """
# BLOSUM62 in 1/2 Bit
# Background (precomputed optional): 0.07422 0.02469 0.05363 0.05431 0.04742 0.07415 0.02621 0.06792 0.05815 0.09891 0.02499 0.04465 0.03854 0.03426 0.05161 0.05723 0.05089 0.07292 0.01303 0.03228 0.00001
# Lambda     (precomputed optional): 0.34657
   A       C       D       E       F       G       H       I       K       L       M       N       P       Q       R       S       T       V       W       Y       X
A  3.9291 -0.4085 -1.7534 -0.8639 -2.2101  0.1596 -1.6251 -1.3218 -0.7340 -1.4646 -0.9353 -1.5307 -0.8143 -0.8040 -1.4135  1.1158 -0.0454 -0.1894 -2.5269 -1.7640 -1.0000
C -0.4085  8.5821 -3.4600 -3.6125 -2.3755 -2.5004 -2.9878 -1.2277 -3.0363 -1.2775 -1.4198 -2.6598 -2.7952 -2.9019 -3.3892 -0.8750 -0.8667 -0.8077 -2.3041 -2.4071 -1.0000
D -1.7534 -3.4600  5.7742  1.5103 -3.4839 -1.3135 -1.1189 -3.1212 -0.7018 -3.6057 -3.0585  1.2717 -1.4801 -0.3134 -1.6058 -0.2610 -1.0507 -3.1426 -4.2143 -3.0650 -1.0000
E -0.8639 -3.6125  1.5103  4.9028 -3.1924 -2.1102 -0.1177 -3.1944  0.7753 -2.8465 -1.9980 -0.2680 -1.1162  1.8546 -0.1154 -0.1469 -0.8633 -2.4423 -2.8354 -2.0205 -1.0000
F -2.2101 -2.3755 -3.4839 -3.1924  6.0461 -3.1074 -1.2342 -0.1609 -3.0787  0.4148  0.0126 -2.9940 -3.5973 -3.1644 -2.7863 -2.3690 -2.1076 -0.8490  0.9176  2.9391 -1.0000
G  0.1596 -2.5004 -1.3135 -2.1102 -3.1074  5.5633 -2.0409 -3.7249 -1.5280 -3.6270 -2.6766 -0.4228 -2.1335 -1.7852 -2.3041 -0.2925 -1.5754 -3.1387 -2.4915 -3.0398 -1.0000
H -1.6251 -2.9878 -1.1189 -0.1177 -1.2342 -2.0409  7.5111 -3.2316 -0.7210 -2.7867 -1.5513  0.5785 -2.1609  0.4480 -0.2499 -0.8816 -1.6859 -3.1175 -2.3422  1.6926 -1.0000
I -1.3218 -1.2277 -3.1212 -3.1944 -0.1609 -3.7249 -3.2316  3.9985 -2.6701  1.5216  1.1268 -3.2170 -2.7567 -2.7696 -2.9902 -2.3482 -0.7176  2.5470 -2.5805 -1.3314 -1.0000
K -0.7340 -3.0363 -0.7018  0.7753 -3.0787 -1.5280 -0.7210 -2.6701  4.5046 -2.4468 -1.3547 -0.1790 -1.0136  1.2726  2.1087 -0.2034 -0.6696 -2.2624 -2.9564 -1.8200 -1.0000
L -1.4646 -1.2775 -3.6057 -2.8465  0.4148 -3.6270 -2.7867  1.5216 -2.4468  3.8494  1.9918 -3.3789 -2.8601 -2.1339 -2.1546 -2.4426 -1.1975  0.7884 -1.6319 -1.0621 -1.0000
M -0.9353 -1.4198 -3.0585 -1.9980  0.0126 -2.6766 -1.5513  1.1268 -1.3547  1.9918  5.3926 -2.1509 -2.4764 -0.4210 -1.3671 -1.4809 -0.6663  0.6872 -1.4248 -0.9949 -1.0000
N -1.5307 -2.6598  1.2717 -0.2680 -2.9940 -0.4228  0.5785 -3.2170 -0.1790 -3.3789 -2.1509  5.6532 -2.0004  0.0017 -0.4398  0.6009 -0.0461 -2.8763 -3.6959 -2.0818 -1.0000
P -0.8143 -2.7952 -1.4801 -1.1162 -3.5973 -2.1335 -2.1609 -2.7567 -1.0136 -2.8601 -2.4764 -2.0004  7.3646 -1.2819 -2.1086 -0.8090 -1.0753 -2.3487 -3.6542 -2.9198 -1.0000
Q -0.8040 -2.9019 -0.3134  1.8546 -3.1644 -1.7852  0.4480 -2.7696  1.2726 -2.1339 -0.4210  0.0017 -1.2819  5.2851  0.9828 -0.1011 -0.6753 -2.1984 -1.9465 -1.4211 -1.0000
R -1.4135 -3.3892 -1.6058 -0.1154 -2.7863 -2.3041 -0.2499 -2.9902  2.1087 -2.1546 -1.3671 -0.4398 -2.1086  0.9828  5.4735 -0.7648 -1.1223 -2.5026 -2.6794 -1.6939 -1.0000
S  1.1158 -0.8750 -0.2610 -0.1469 -2.3690 -0.2925 -0.8816 -2.3482 -0.2034 -2.4426 -1.4809  0.6009 -0.8090 -0.1011 -0.7648  3.8844  1.3811 -1.6462 -2.7519 -1.6858 -1.0000
T -0.0454 -0.8667 -1.0507 -0.8633 -2.1076 -1.5754 -1.6859 -0.7176 -0.6696 -1.1975 -0.6663 -0.0461 -1.0753 -0.6753 -1.1223  1.3811  4.5453 -0.0555 -2.4289 -1.6060 -1.0000
V -0.1894 -0.8077 -3.1426 -2.4423 -0.8490 -3.1387 -3.1175  2.5470 -2.2624  0.7884  0.6872 -2.8763 -2.3487 -2.1984 -2.5026 -1.6462 -0.0555  3.7689 -2.8343 -1.2075 -1.0000
W -2.5269 -2.3041 -4.2143 -2.8354  0.9176 -2.4915 -2.3422 -2.5805 -2.9564 -1.6319 -1.4248 -3.6959 -3.6542 -1.9465 -2.6794 -2.7519 -2.4289 -2.8343 10.5040  2.1542 -1.0000
Y -1.7640 -2.4071 -3.0650 -2.0205  2.9391 -3.0398  1.6926 -1.3314 -1.8200 -1.0621 -0.9949 -2.0818 -2.9198 -1.4211 -1.6939 -1.6858 -1.6060 -1.2075  2.1542  6.5950 -1.0000
X -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000 -1.0000
"""

AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']

backprobs = [0.0489372, 0.0306991, 0.101049, 0.0329671, 0.0276149,
             0.0416262, 0.0452521, 0.030876, 0.0297251, 0.0607036,
             0.0150238, 0.0215826, 0.0783843, 0.0512926, 0.0264886,
             0.0610702, 0.0201311, 0.215998, 0.0310265, 0.0295417,
             0.00001]
backprobs = np.asarray(backprobs, dtype=np.float32)[:20] # done need X
bitFactor = 8

def parse_profile_from_block(block_lines):
    # Skip the first two header lines (header and amino acid header)
    data_lines = block_lines[2:]
    profile_rows = []
    for line in data_lines:
        line = line.strip()
        if line:  # ignore blank lines
            values = list(map(float, line.split()))
            if len(values) != 20:
                print(f"Warning: expected 20 values per line but got {len(values)}. Skipping line: {line}")
                continue
            profile_rows.append(values)
    profile_matrix = np.array(profile_rows)
    return profile_matrix

def parse_profiles(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    profiles = []
    current_block = []
    header_pattern = re.compile(r"Query profile of sequence\s+(.*)")
    
    for line in tqdm(lines, desc="Parsing profile lines"):
        if header_pattern.match(line):
            if current_block:
                profiles.append(current_block)
                current_block = []
        current_block.append(line)
    if current_block:
        profiles.append(current_block)
    
    profile_data = []
    for block in tqdm(profiles, desc="Processing profile blocks"):
        if len(block) < 3:
            continue
        header_line = block[0].strip()
        m = header_pattern.match(header_line)
        key = m.group(1).strip() if m else ""
        profile_matrix = parse_profile_from_block(block)
        profile_data.append((key, profile_matrix))
    
    return profile_data


def computeLogPSSM(profile_matrix):
    consensus = np.argmax(profile_matrix, axis=1).astype(np.int8)

    odds = profile_matrix / backprobs[None, :]

    # log2 only where odds > 0
    log_probs = np.empty_like(odds, dtype=np.float32)
    log_probs.fill(-128.0)          # default for invalid
    mask = odds > 0.0
    log_probs[mask] = np.log2(odds[mask])

    pssm = log_probs * bitFactor

    # round exactly like original
    pssm = np.where(pssm < 0, pssm - 0.5, pssm + 0.5)

    # clip and cast
    pssm = np.clip(pssm, -128, 127).astype(np.int8)

    return pssm.ravel(), consensus


def toBuffer_pssm(profile, consensus):
    L = len(consensus)
    buf = bytearray(L * 25)

    offset = 0
    prof = memoryview(profile)

    for i in range(L):
        # 20 profile bytes
        buf[offset:offset+20] = prof[i*20:(i+1)*20]
        offset += 20

        c = consensus[i]
        buf[offset]   = c
        buf[offset+1] = c
        buf[offset+2] = 0
        buf[offset+3] = 0
        buf[offset+4] = 0
        offset += 5

    return bytes(buf)

#############################
# seq_to_db functionality
#############################

# Note: This version uses a similar approach (collecting packed bytes in a list)
# but with modifications (e.g. scaling scores by 4) as in your provided code.
def parse_substitution_matrix_seq(matrix_str):
    lines = matrix_str.strip().splitlines()
    header = None
    mat = {}
    # for line in tqdm(lines, desc="Parsing substitution matrix"):
    for line in lines: # no tqdm
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if header is None:
            tokens = line.split()
            header = tokens[:20]
        else:
            tokens = line.split()
            row_letter = tokens[0]
            scores = list(map(float, tokens[1:21]))
            mat[row_letter.upper()] = scores
    return header, mat

# Use the same AAs list as before.
# def generate_profile_for_sequence(seq, sub_matrix):
#     profile = []
#     for res in seq:
#         r = res.upper()
#         if r not in sub_matrix:
#             r = 'X'
#         # In your seq_to_db.py version, if r=='X' you skip the residue.
#         if r == 'X':
#             continue
#         row = sub_matrix[r]
#         for score in row:
#             # Scale the score by 4
#             scaled = score * 4
#             if scaled < 0:
#                 val = int(scaled - 0.5)
#             else:
#                 val = int(scaled + 0.5)
#             val = max(-128, min(val, 127))
#             profile.append(np.int8(val))
#         profile.append(np.int8(AAs.index(r)))
#         profile.append(np.int8(AAs.index(r)))
#         profile.append(np.int8(0))
#         profile.append(np.int8(0))
#         profile.append(np.int8(0))
#     return profile

# def pack_profile_seq(profile):
#     parts = []
#     for val in profile:
#         parts.append(struct.pack('b', val))
#     return b"".join(parts)

def pack_profile_seq(profile):
    # vectorised
    profile = np.asarray(profile, dtype=np.int8)
    return profile.tobytes()

def read_sequences(seq_db_file, seq_index_file):
    with open(seq_db_file, "rb") as f:
        db_data = f.read()
    sequences = []
    with open(seq_index_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            key = parts[0]
            offset = int(parts[1])
            length = int(parts[2])
            seq_bytes = db_data[offset: offset+length]
            seq = seq_bytes.decode("ascii").strip()
            # Remove last two characters as in original code.
            sequences.append((key, seq[:-2]))
    # logger.info("First 5 sequences:", sequences[:5])
    return sequences

AA_TO_INDEX = {aa: i for i, aa in enumerate(AAs)}

def precompute_sub_blocks(sub_matrix):
    blocks = {}

    for aa, row in sub_matrix.items():
        if aa == 'X':
            continue

        row = np.asarray(row, dtype=np.float32)

        scaled = row * 4.0
        scaled = np.where(scaled < 0, scaled - 0.5, scaled + 0.5)
        scaled = np.clip(scaled, -128, 127).astype(np.int8)

        idx = AA_TO_INDEX[aa]

        block = bytearray(25)
        block[0:20] = scaled.view(np.uint8).tobytes()
        block[20] = idx
        block[21] = idx
        block[22] = 0
        block[23] = 0
        block[24] = 0

        blocks[aa] = bytes(block)

    return blocks


def generate_profile_for_sequence(seq, sub_blocks):
    out = bytearray()

    for res in seq:
        aa = res.upper()
        block = sub_blocks.get(aa)
        if block is None:
            continue  #  X-skip 
        out.extend(block)

    return out

def build_database_seq(seq_db_file, seq_index_file, output_db="profile", output_index="profile.idx"):
    header, sub_matrix = parse_substitution_matrix_seq(blosum62_str)
    seqs = read_sequences(seq_db_file, seq_index_file)
    parts = []
    index_lines = []
    current_offset = 0
    sub_blocks = precompute_sub_blocks(sub_matrix)
    # for key, seq in tqdm(seqs, desc="Processing sequences"):
    for key, seq in seqs:
        profile = generate_profile_for_sequence(seq, sub_blocks)
        packed = pack_profile_seq(profile)
        parts.append(packed)
        length = len(packed)
        index_lines.append(f"{key}\t{current_offset}\t{length}")
        current_offset += length
    db_buffer = b"".join(parts)
    with open(output_db, "wb") as f:
        f.write(db_buffer)
    with open(output_index, "w") as f:
        f.write("\n".join(index_lines) + "\n")
    # logger.info(f"seq database written to '{output_db}'")
    # logger.info(f"seq index written to '{output_index}'")

#############################
# Additional file copying & index sorting
#############################

def copy_and_create_extras(mmseqs_db):

    # Source file paths (assumed to exist)
    src_h = f"{mmseqs_db}_h"
    src_h_index = f"{mmseqs_db}_h.index"
    src_h_dbtype = f"{mmseqs_db}_h.dbtype"
    src_lookup = f"{mmseqs_db}.lookup"
    
    # Destination file paths:
    dest_profile_h = f"{mmseqs_db}_profile_h"
    dest_profile_h_index =f"{mmseqs_db}_profile_h.index"
    dest_profile_h_dbtype =  f"{mmseqs_db}_profile_h.dbtype"
    dest_lookup = f"{mmseqs_db}_profile.lookup"
    dest_profile_dbtype =  f"{mmseqs_db}_profile.dbtype"
    dest_profile_ss_dbtype =  f"{mmseqs_db}_profile_ss.dbtype"
    
    shutil.copy2(src_h, dest_profile_h)
    shutil.copy2(src_h_index, dest_profile_h_index)
    shutil.copy2(src_h_dbtype, dest_profile_h_dbtype)
    # shutil.copy2(src_lookup, dest_lookup)
    
    # Create the .dbtype files by writing 4 bytes: (2,0,0,0)
    data = bytes([2,0,0,0])
    with open(dest_profile_dbtype, "wb") as f:
        f.write(data)
    with open(dest_profile_ss_dbtype, "wb") as f:
        f.write(data)
    
    # logger.info("Extra files copied and .dbtype files created.")

def sort_index_file(index_path):
    # Read index lines, sort them numerically by the first field.
    with open(index_path, "r") as f:
        lines = f.readlines()
    # Try to sort by converting the first column to int (if possible)
    try:
        sorted_lines = sorted(lines, key=lambda x: int(x.split()[0]))
    except ValueError:
        # If conversion fails, sort lexicographically.
        sorted_lines = sorted(lines)
    with open(index_path, "w") as f:
        f.writelines(sorted_lines)
    # logger.info(f"Sorted index file: {index_path}")