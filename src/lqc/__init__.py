import importlib

__version__ = "0.0.9"

__all__ = [
    'CS', 'Indel', 'Mismatch', 'Splice',
    'bam_or_sam', 'check_bam_with_cs_or_md', 'copy_logo',
    'create_indel_summary_table', 'create_mapping_table',
    'create_mismatch_normalized_read_location_table',
    'create_readstat_table', 'create_splice_all_table',
    'create_splice_table',
    'get_html_template', 'html_add_bootstrap', 'html_add_data',
    'html_add_deletion_table', 'html_add_insertion_table',
    'html_add_mapping', 'html_add_mismatch_table',
    'html_add_readstat_table', 'html_add_splice_table',
    'inline_figures', 'list_bam_contigs', 'open_alignment_file',
    'plot_element_total_count', 'plot_indel_hist_length',
    'plot_indel_hist_location', 'plot_mapping_aligned_fraction_hist',
    'plot_mapping_aligned_vs_query', 'plot_mapping_mapq_hist',
    'plot_mismatch_hist_location', 'plot_mismatch_type_count',
    'plot_readstat_bar', 'plot_readstat_bar_mean_element_per_read',
    'plot_readstat_bar_mean_element_per_read_per_kb',
    'plot_readstat_bar_ratio_with_element',
    'plot_readstat_cumulative_length', 'plot_readstat_length_hist',
    'plot_splice_type_count', 'stat_element_from_bam_by_contig',
    'write_readcs',
]

_LAZY = {
    'CS': 'lqc.cs',
    'Indel': 'lqc.indel',
    'Mismatch': 'lqc.mismatch',
    'Splice': 'lqc.splice',
    'bam_or_sam': 'lqc.utils',
    'check_bam_with_cs_or_md': 'lqc.utils',
    'list_bam_contigs': 'lqc.utils',
    'open_alignment_file': 'lqc.utils',
    'write_readcs': 'lqc.utils',
    'stat_element_from_bam_by_contig': 'lqc.stat',
    'copy_logo': 'lqc.template',
    'get_html_template': 'lqc.template',
    'create_indel_summary_table': 'lqc.report_table',
    'create_mapping_table': 'lqc.report_table',
    'create_mismatch_normalized_read_location_table': 'lqc.report_table',
    'create_readstat_table': 'lqc.report_table',
    'create_splice_all_table': 'lqc.report_table',
    'create_splice_table': 'lqc.report_table',
    'html_add_bootstrap': 'lqc.report_html',
    'html_add_data': 'lqc.report_html',
    'html_add_deletion_table': 'lqc.report_html',
    'html_add_insertion_table': 'lqc.report_html',
    'html_add_mapping': 'lqc.report_html',
    'html_add_mismatch_table': 'lqc.report_html',
    'html_add_readstat_table': 'lqc.report_html',
    'html_add_splice_table': 'lqc.report_html',
    'inline_figures': 'lqc.report_html',
    'plot_element_total_count': 'lqc.report_figure',
    'plot_indel_hist_length': 'lqc.report_figure',
    'plot_indel_hist_location': 'lqc.report_figure',
    'plot_mapping_aligned_fraction_hist': 'lqc.report_figure',
    'plot_mapping_aligned_vs_query': 'lqc.report_figure',
    'plot_mapping_mapq_hist': 'lqc.report_figure',
    'plot_mismatch_hist_location': 'lqc.report_figure',
    'plot_mismatch_type_count': 'lqc.report_figure',
    'plot_readstat_bar': 'lqc.report_figure',
    'plot_readstat_bar_mean_element_per_read': 'lqc.report_figure',
    'plot_readstat_bar_mean_element_per_read_per_kb': 'lqc.report_figure',
    'plot_readstat_bar_ratio_with_element': 'lqc.report_figure',
    'plot_readstat_cumulative_length': 'lqc.report_figure',
    'plot_readstat_length_hist': 'lqc.report_figure',
    'plot_splice_type_count': 'lqc.report_figure',
}


def __getattr__(name):
    module_name = _LAZY.get(name)
    if module_name is not None:
        value = getattr(importlib.import_module(module_name), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'lqc' has no attribute {name!r}")


def __dir__():
    return list(globals()) + list(_LAZY)
