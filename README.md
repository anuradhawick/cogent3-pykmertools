# `cogent3_pykmertools`: kmer counting plugins for `cogent3`

## Installation

We recommend installing `cogent3_pykmertools` from PyPI as follows

```
pip install cogent3_pykmertools
```

Or using `uv`

```
uv pip install cogent3_pykmertools
```

The installation process registers the apps and plugins with `cogent3`.

## The provided apps

- `pkt_count_kmers()`: does what the name says! When this package is installed, this app will be used by the `cogent3` `Sequence.count_kmers()` and `SequenceCollection.count_kmers()` methods. This app can also be run in parallel, which can be useful if you have a lot of sequences in the `SequenceCollection`.
- `pkt_kmer_header`: returns the k-mers as strings, so you know what their order is.

## Getting help on the apps

Use the `cogent3.app_help()` function. For example, `cogent3.app_help("pkt_count_kmers")`.
