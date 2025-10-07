

import pytest
from cogent3 import get_app
from count quick import Pkt_Count_Kmers

@pytest.mark.xfail(reason="Expected to fail because it is using cookiecutter values")
def test_pkt_count_kmers_installed():
    app = get_app("Pkt_Count_Kmers")
    got = app("test")
    expected = None  # replace with expected result
    assert got == expected, f"got {got}, expected {expected}"
