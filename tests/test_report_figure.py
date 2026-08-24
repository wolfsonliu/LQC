import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

matplotlib.use('Agg')

from lqc.indel import Indel
from lqc.readstat import ReadStat
from lqc.report_figure import (
    determine_figure_size,
    get_facet_row_col,
    plot_indel_hist_length,
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


def test_plot_readstat_bar_multi_contig_cycles_palette():
    rs = []
    for i in range(5):
        r = ReadStat(f'chr{i}')
        r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
        rs.append(r)
    fig = plot_readstat_bar(rs, 'Read count')
    bars = fig.axes[0].patches
    colors = {p.get_facecolor() for p in bars}
    assert len(colors) > 2
    plt.close('all')


def test_plot_indel_hist_length_readable():
    indels = [Indel('chr1')]
    # one short and one long indel: the long tail must not hide the short bin
    indels[0].add_indel('a', 0.1)
    indels[0].add_indel('t' * 1000, 0.5)
    fig = plot_indel_hist_length(indels, width=5, height=4)
    ax = fig.axes[0]
    assert ax.get_yscale() == 'log'
    plt.close('all')