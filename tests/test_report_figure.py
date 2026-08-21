import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

matplotlib.use('Agg')

from lqc.readstat import ReadStat
from lqc.report_figure import (
    determine_figure_size,
    get_facet_row_col,
    plot_readstat_bar,
)


def test_get_facet_row_col():
    assert get_facet_row_col(1) == (1, 1)
    assert get_facet_row_col(4) == (2, 2)
    assert get_facet_row_col(5) == (2, 3)


def test_determine_figure_size():
    assert determine_figure_size(2, 3) == (12, 6)


def test_plot_readstat_bar_returns_figure():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
    fig = plot_readstat_bar([r], 'Read count')
    assert isinstance(fig, Figure)
    plt.close('all')