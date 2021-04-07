import pytest

import vvhgvs.easy


@pytest.fixture(scope="session")
def parser():
    return vvhgvs.easy.parser

@pytest.fixture(scope="session")
def am38():
    return vvhgvs.easy.am38


