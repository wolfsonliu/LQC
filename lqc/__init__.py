__version__ = "0.0.5"

from lqc.cs import CS
from lqc.indel import Indel
from lqc.mismatch import Mismatch
from lqc.report_figure import (
    plot_element_total_count,
    plot_indel_hist_length,
    plot_indel_hist_location,
    plot_mismatch_hist_location,
    plot_mismatch_type_count,
    plot_readstat_bar,
    plot_readstat_bar_mean_element_per_read,
    plot_readstat_bar_mean_element_per_read_per_kb,
    plot_readstat_bar_ratio_with_element,
    plot_readstat_cumulative_length,
    plot_readstat_length_hist,
    plot_splice_type_count,
)
from lqc.report_html import (
    html_add_bootstrap,
    html_add_data,
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mismatch_table,
    html_add_readstat_table,
    html_add_splice_table,
    inline_figures,
)
from lqc.report_table import (
    create_indel_summary_table,
    create_mismatch_normalized_read_location_table,
    create_readstat_table,
    create_splice_table,
)
from lqc.splice import Splice
from lqc.stat import stat_element_from_bam_by_contig
from lqc.template import copy_logo, get_html_template
from lqc.utils import (
    bam_or_sam,
    check_bam_with_cs_or_md,
    list_bam_contigs,
    write_readcs,
)

__all__ = [
    'CS',
    'Indel',
    'Mismatch',
    'Splice',
    'bam_or_sam',
    'check_bam_with_cs_or_md',
    'copy_logo',
    'create_indel_summary_table',
    'create_mismatch_normalized_read_location_table',
    'create_readstat_table',
    'create_splice_table',
    'get_html_template',
    'html_add_bootstrap',
    'html_add_data',
    'html_add_deletion_table',
    'html_add_insertion_table',
    'html_add_mismatch_table',
    'html_add_readstat_table',
    'html_add_splice_table',
    'inline_figures',
    'list_bam_contigs',
    'plot_element_total_count',
    'plot_indel_hist_length',
    'plot_indel_hist_location',
    'plot_mismatch_hist_location',
    'plot_mismatch_type_count',
    'plot_readstat_bar',
    'plot_readstat_bar_mean_element_per_read',
    'plot_readstat_bar_mean_element_per_read_per_kb',
    'plot_readstat_bar_ratio_with_element',
    'plot_readstat_cumulative_length',
    'plot_readstat_length_hist',
    'plot_splice_type_count',
    'stat_element_from_bam_by_contig',
    'write_readcs'
]
