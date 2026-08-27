import subprocess
import sys

import pytest

import lqc


def test_all_names_resolve_lazily():
    for name in lqc.__all__:
        assert getattr(lqc, name) is not None


def test_unknown_attribute_raises():
    with pytest.raises(AttributeError):
        lqc.__getattr__('definitely_not_real')


def test_submodules_are_not_eagerly_imported():
    code = (
        'import sys, lqc;'
        'print("figure", "lqc.report_figure" in sys.modules);'
        'print("html", "lqc.report_html" in sys.modules)'
    )
    out = subprocess.check_output(
        [sys.executable, '-c', code], text=True
    )
    assert 'figure False' in out
    assert 'html False' in out