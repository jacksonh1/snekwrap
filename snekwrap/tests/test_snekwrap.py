"""
Unit and regression test for the snekwrap package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import snekwrap


def test_snekwrap_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "snekwrap" in sys.modules
