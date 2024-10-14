from pathlib import Path
from typing import Dict
from urllib.error import HTTPError, URLError

import bioframe as bf
import pandas as pd

from .convert import chromsizes_to_dict
from .schemas import bedgraph_dtypes, bedgraph_schema


def read_bedgraph(bg_filename: str) -> pd.DataFrame:
    """Loads bedgraph file with specific dtypes into pd.DataFrame."""
    return pd.read_csv(bg_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=bedgraph_schema,
                       dtype=bedgraph_dtypes)


def read_annotation(filename, schema='bed6'):
    return bf.read_table(filename, schema=schema)


def _read_chromsizes(chromsizes_filename: str) -> pd.Series:
    return pd.read_csv(chromsizes_filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=('chrom', 'length')).set_index('chrom').squeeze()


def _fetch_chromsizes(genome: str) -> pd.Series:
    return bf.fetch_chromsizes(genome)


def get_chromsizes(chromsizes: str) -> Dict[str, int]:
    if Path(chromsizes).exists():
        chromsizes_dict = chromsizes_to_dict(_read_chromsizes(chromsizes))
    else:
        try:
            chromsizes_dict = chromsizes_to_dict(_fetch_chromsizes(chromsizes))
        except HTTPError:
            raise ValueError(f'{chromsizes} file does not exist and is not a valid UCSC genome name.')
        except URLError:
            raise Exception(f"Couldn't fetch chromsizes for {chromsizes}, check internet connection.")
    return chromsizes_dict
