import pytest
import os

def pytest_runtest_teardown():
    from vmd import molecule
    for _ in molecule.listall():
        molecule.delete(_)
