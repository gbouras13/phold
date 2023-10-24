#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from loguru import logger

def get_tophits( result_tsv: Path, top_hits_tsv: Path, evalue: float  ) ->  pd.DataFrame:

    logger.info("Processing Foldseek output.")
    col_list = [
            "query",
            "target",
            "foldseek_alnScore",
            "foldseek_seqIdentity",
            "foldseek_eVal",
            "qStart",
            "qEnd",
            "qLen",
            "tStart",
            "tEnd",
            "tLen",
        ]
    
    foldseek_df = pd.read_csv(
            result_tsv, delimiter="\t", index_col=False, names=col_list
        )


    # Group by 'gene' and find the top hit for each group
    tophits_df = (
            foldseek_df.groupby("query", group_keys=True)
            .apply(lambda group: group.nsmallest(1, "foldseek_eVal"))
            .reset_index(drop=True)
        )
    
    # scientific notation to 3dp
    tophits_df['foldseek_eVal'] = tophits_df['foldseek_eVal'].apply(lambda x: '{:.3e}'.format(float(x)))

    # tolerance = 1e-10
    # filtered_tophits_df = tophits_df[tophits_df['foldseek_eVal'] < evalue + tolerance]
    # print(tophits_df)


    
    # save tophits
    tophits_df.to_csv(
        top_hits_tsv, sep="\t", index=False
    )

    return tophits_df
    
