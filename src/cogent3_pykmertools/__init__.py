from cogent3.app import typing as c3types
from cogent3.app.composable import define_app
import pykmertools as pkt
import numpy as np

__version__ = "0.1.0"


@define_app
def pkt_count_kmers(
    seqs: list[c3types.SeqType],
    k: int = 3,
    count_min_complements: bool = False,
    parallel: bool = False,
) -> np.ndarray:
    """compute counts of overlapping k-mers

    Parameters
    ----------
    seqs
        list of cogent3 sequence objects
    k
        size of the kmer
    count_min_complements
        merge counts from strand complementary k-mers
    parallel
        whether to run in parallel, by default False

    Returns
    -------
        a 2D array of k-mer counts
    """
    ctr = pkt.OligoComputer(k)
    raw_seqs = [str(seq) for seq in seqs]
    if parallel:
        counts = np.array(ctr.vectorise_batch(raw_seqs, False, count_min_complements))
    else:
        counts = np.array(
            [ctr.vectorise_one(s, False, count_min_complements) for s in raw_seqs],
        )
    return counts.astype(int)


@define_app
def pkt_kmer_header(
    k: int = 3,
    count_min_complements: bool = False,
) -> np.ndarray:
    """returns the numpy string array of the kmers

    Parameters
    ----------
    k
        size of the kmer
    count_min_complements
        merge counts from strand complementary k-mers

    Returns
    -------
        a 1D array of kmers
    """
    ctr = pkt.OligoComputer(k)
    return ctr.get_header(count_min_complements)
