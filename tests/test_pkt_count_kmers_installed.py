import pytest
import collections
import cogent3 as c3
import cogent3_pykmertools as c3pkt


def test_pkt_count_kmers_installed():
    app = c3.get_app("pkt_count_kmers")
    assert str(app).startswith("pkt_count_kmers")


@pytest.fixture
def primates():
    return c3.get_dataset("primate-brca1").degap()


@pytest.mark.parametrize("k", [1, 2])
def test_pkt_count_kmers(k):
    data = "AACGTTTCG"
    expect = collections.Counter([data[i : i + k] for i in range(len(data) - k + 1)])
    seq = c3.make_seq(data, moltype="dna")
    app = c3pkt.pkt_count_kmers(k=k)
    counts = app.main([seq])[0]

    # the set of counts should match
    expect_counts = set(expect.values())
    expect_counts = (
        expect_counts if k == 1 else expect_counts | {0}
    )  # k=1 has all kmers
    assert set(counts.tolist()) == expect_counts

    # the kmers mapping to their counts should match
    header = c3pkt.pkt_kmer_header()(k)
    # filter out zero counts
    result = {h: c for h, c in zip(header, counts, strict=True) if c}
    assert result == expect


def test_integration_with_cogent3_seq(primates):
    k = 2
    seq = primates.seqs["Human"]
    c3_counts = seq.count_kmers(k=k, use_hook="cogent3")
    c3_header = seq.moltype.alphabet.get_kmer_alphabet(k=k)
    expect = dict(zip(c3_header, c3_counts, strict=True))

    # expect them to have a different ordering
    pkt_header = c3pkt.pkt_kmer_header()(k)
    assert c3_header != pkt_header
    pkt_counts = seq.count_kmers(k=k, use_hook="cogent3_pykmertools")
    # so their element wise comparisons will differ
    assert (pkt_counts != c3_counts).any()
    # expect their counts to be identical once mapped to kmer strings
    assert dict(zip(pkt_header, pkt_counts, strict=True)) == expect


def test_integration_with_cogent3_seqcoll(primates):
    seqs = primates
    k = 2
    pkt_header = c3pkt.pkt_kmer_header()(k)
    c3_header = seqs.moltype.alphabet.get_kmer_alphabet(k=k)
    pkt_to_c3 = [pkt_header.index(h) for h in c3_header]
    pkt_counts = seqs.count_kmers(k=k, use_hook="cogent3_pykmertools")
    c3_counts = seqs.count_kmers(k=k, use_hook="cogent3")
    # expect the same shape
    assert pkt_counts.shape == c3_counts.shape
    # but different ordering
    assert (pkt_counts != c3_counts).any()
    # but then re-order the pkt counts to match c3 counts
    reordered_pkt_counts = pkt_counts[:, pkt_to_c3]
    assert (reordered_pkt_counts == c3_counts).all()


def test_cogent3_seqcoll_parallel(primates):
    seqs = primates
    serial = seqs.count_kmers(k=2, use_hook="cogent3_pykmertools")
    parallel = seqs.count_kmers(k=2, use_hook="cogent3_pykmertools", parallel=True)
    assert (serial == parallel).all()
