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
    plot_mapping_aligned_fraction_hist,
    plot_mapping_aligned_vs_query,
    plot_mapping_mapq_hist,
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
    for i in range(9):
        r = ReadStat(f'chr{i}')
        r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
        rs.append(r)
    fig = plot_readstat_bar(rs, 'Read count')
    bars = fig.axes[0].patches
    colors = {p.get_facecolor() for p in bars}
    # 9 bars over an 8-color palette wrap back: exactly 8 distinct colors
    assert len(colors) == 8
    plt.close('all')


def test_plot_indel_hist_length_readable():
    indels = [Indel('chr1')]
    # one short and one long indel: the long tail must not hide the short bin
    indels[0].add_indel('a', 0.1)
    indels[0].add_indel('t' * 1000, 0.5)
    fig = plot_indel_hist_length(indels, width=5, height=4)
    ax = fig.axes[0]
    assert ax.get_yscale() == 'log'
    heights = [p.get_height() for p in ax.patches]
    # one integer-width bin per length (1..1000), no catch-all tail bin
    assert len(heights) == 1000
    assert heights[0] == 1
    assert heights[-1] == 1
    assert sum(heights) == 2
    plt.close('all')


def test_plot_indel_hist_length_one_bin_per_length():
    indels = [Indel('chr1')]
    for length in range(1, 21):
        indels[0].add_indel('a' * length, 0.5)
    fig = plot_indel_hist_length(indels, width=5, height=4)
    ax = fig.axes[0]
    heights = [p.get_height() for p in ax.patches]
    # 20 distinct lengths -> 20 bins, each with exactly one count (no plateau)
    assert len(heights) == 20
    assert heights == [1] * 20
    plt.close('all')


def test_plot_indel_hist_length_empty_indel_renders_blank():
    indels = [Indel('chr1')]  # no indels -> empty lengths
    fig = plot_indel_hist_length(indels, width=5, height=4)
    ax = fig.axes[0]
    assert len(ax.patches) == 0
    plt.close('all')


def _one_mapped_readstat():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
               mapping_quality=60, aligned_length=180)
    return r


def test_plot_mapping_mapq_hist_returns_figure():
    fig = plot_mapping_mapq_hist([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')


def test_plot_mapping_aligned_fraction_hist_returns_figure():
    fig = plot_mapping_aligned_fraction_hist([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')


def test_plot_mapping_aligned_vs_query_returns_figure():
    fig = plot_mapping_aligned_vs_query([_one_mapped_readstat()])
    assert isinstance(fig, Figure)
    plt.close('all')


def test_plot_mapping_mapq_hist_caps_high_mapq():
    r = ReadStat('chr1')
    r.add_read(200, insertion=0, deletion=0, mismatch=0, intron=0,
               mapping_quality=255, aligned_length=200)
    fig = plot_mapping_mapq_hist([r])
    assert len(fig.axes[0].patches) < 100
    plt.close('all')


def test_plot_mapping_hist_empty_readstat_renders_blank():
    empty = ReadStat('chr1')
    fig = plot_mapping_mapq_hist([empty])
    assert len(fig.axes[0].patches) == 0
    fig2 = plot_mapping_aligned_fraction_hist([empty])
    assert len(fig2.axes[0].patches) == 0
    plt.close('all')