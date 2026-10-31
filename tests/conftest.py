import random

import pytest


@pytest.fixture
def rng():
    return random.Random(12345)


@pytest.fixture
def reference(rng):
    """A random 5 kb 'genome' (random DNA has almost no repeats)."""
    return "".join(rng.choice("ACGT") for _ in range(5000))
