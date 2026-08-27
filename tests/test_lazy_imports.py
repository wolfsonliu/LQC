import subprocess
import sys

import pytest

import lqc


def test_all_names_resolve_lazily():
    for name in lqc.__all__:
        assert getattr(lqc, name) is not None


def test_unknown_attribute_raises():
    name = 'definitely_not_real'
    with pytest.raises(AttributeError):
        getattr(lqc, name)


def test_public_surface_matches_lazy_registry():
    assert set(lqc._LAZY) == set(lqc.__all__)


def test_dir_exposes_lazy_names():
    assert set(lqc.__all__) <= set(dir(lqc))


def test_submodules_are_not_eagerly_imported():
    code = (
        'import sys, lqc;'
        'print("figure", "lqc.report_figure" in sys.modules);'
        'print("html", "lqc.report_html" in sys.modules);'
        'print("stat", "lqc.stat" in sys.modules)'
    )
    out = subprocess.check_output(
        [sys.executable, '-c', code], text=True
    )
    assert 'figure False' in out
    assert 'html False' in out
    assert 'stat False' in out
