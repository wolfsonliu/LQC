"""Command-line interface for LQC.

All orchestration lives in main(). See README.md for usage.
"""

import argparse
import logging
import multiprocessing as mp
import os
import pickle
import sys
from functools import partial

import matplotlib
import matplotlib.pyplot as plt

from lqc import __version__ as VERSION
from lqc import (
    check_bam_with_cs_or_md,
    copy_logo,
    create_indel_summary_table,
    create_mapping_table,
    create_mismatch_normalized_read_location_table,
    create_readstat_table,
    create_splice_all_table,
    create_splice_table,
    get_html_template,
    html_add_bootstrap,
    html_add_data,
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mapping,
    html_add_mismatch_table,
    html_add_readstat_table,
    html_add_splice_table,
    inline_figures,
    list_bam_contigs,
    plot_element_total_count,
    plot_indel_hist_length,
    plot_indel_hist_location,
    plot_mapping_aligned_fraction_hist,
    plot_mapping_aligned_vs_query,
    plot_mapping_mapq_hist,
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
from lqc._base import concat_stats
from lqc.constants import TOTAL_LABEL
from lqc.formatting import FLOAT_FORMAT
from lqc.report_table import total_row
from lqc.stat import plan_tasks, reduce_blocks_to_contigs, stat_region

matplotlib.use('Agg')

logger = logging.getLogger(__name__)


def build_directories(dir_dict):
    for b in dir_dict.values():
        os.makedirs(b, exist_ok = True)


def savefig(fig, prefix):
    fig.savefig(prefix + '.png')
    fig.savefig(prefix + '.pdf')


READSTAT_BAR_SPECS = [
    ('Read count', 'readstat_bar_Read_count'),
    ('Median read length', 'readstat_bar_Median_read_length'),
    ('Mean read length', 'readstat_bar_Mean_read_length'),
    ('Insertions per read', 'readstat_bar_insertions_per_read'),
    ('Insertions per read per kb', 'readstat_bar_insertions_per_read_per_kb'),
    ('Deletions per read', 'readstat_bar_deletions_per_read'),
    ('Deletions per read per kb', 'readstat_bar_deletions_per_read_per_kb'),
    ('Mismatches per read', 'readstat_bar_mismatches_per_read'),
    ('Mismatches per read per kb', 'readstat_bar_mismatches_per_read_per_kb'),
    ('Mean intron number', 'readstat_bar_Mean_intron_number'),
    ('N50', 'readstat_bar_N50'),
    ('L50', 'readstat_bar_L50'),
]

MULTI_FIG_SPECS = [
    ('readstat_bar_mean_element_per_read', plot_readstat_bar_mean_element_per_read, 'readstat'),
    ('readstat_bar_mean_element_per_read_per_kb', plot_readstat_bar_mean_element_per_read_per_kb, 'readstat'),
    ('readstat_line_cumulative_length', plot_readstat_cumulative_length, 'readstat'),
    ('readstat_bar_ratio_with_element', plot_readstat_bar_ratio_with_element, 'readstat'),
    ('readstat_hist_length', plot_readstat_length_hist, 'readstat'),
    ('insertion_hist_length', plot_indel_hist_length, 'insertion'),
    ('deletion_hist_length', plot_indel_hist_length, 'deletion'),
    ('insertion_hist_location', plot_indel_hist_location, 'insertion'),
    ('deletion_hist_location', plot_indel_hist_location, 'deletion'),
    ('mismatch_type', plot_mismatch_type_count, 'mismatch'),
    ('mismatch_hist_location', plot_mismatch_hist_location, 'mismatch'),
    ('splice_type', plot_splice_type_count, 'splice'),
    ('mapping_hist_mapq', plot_mapping_mapq_hist, 'readstat'),
    ('mapping_hist_aligned_fraction', plot_mapping_aligned_fraction_hist, 'readstat'),
    ('mapping_scatter_aligned_vs_query', plot_mapping_aligned_vs_query, 'readstat'),
]

ELEMENT_BAR_SPECS = [
    ('insertion_bar_count', 'insertion', 'Insertion'),
    ('deletion_bar_count', 'deletion', 'Deletion'),
    ('mismatch_bar_count', 'mismatch', 'Mismatch'),
    ('intron_bar_count', 'splice', 'Intron'),
]


def generate_multiple_figs(plot_func,
                           data_list, data_sum,
                           filelabel,
                           width = 5, height = 4):
    fig = plot_func(data_list)
    savefig(fig, filelabel)
    fig = plot_func(
        [data_sum],
        width = width, height = height
    )
    savefig(fig, filelabel + '.' + 'Total')
    for i in range(len(data_list)):
        seqname = data_list[i].label
        fig = plot_func(
            [data_list[i]],
            width = width, height = height
        )
        savefig(
            fig, filelabel + '.' + seqname
        )
        plt.close('all')


def main(argv = None) -> int:

    chromosomes = [
        'chr1', 'chr2', 'chr3', 'chr4',
        'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16',
        'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX', 'chrY'
    ]
    parser = argparse.ArgumentParser(
        prog = 'lqc',
        description='The Long-read RNA-seq quality control software.'
    )
    parser.add_argument(
        '-b', '--bam-file',
        help = 'input bam file, with cs tags, sorted and indexed',
        type = str,
        default = None,
        required = True
    )
    parser.add_argument(
        '--genome-fasta',
        help = 'path of genome fasta file',
        type = str,
        default = None,
        required = False
    )
    parser.add_argument(
        '-o', '--output_dir',
        help = 'directory to store output files',
        type = str,
        default = 'out'
    )
    parser.add_argument(
        '--output-cs',
        help = 'output processed cs tags',
        action = 'store_true'
    )
    parser.add_argument(
        '--output-pickle',
        help = 'output pickle file of results',
        action = 'store_true'
    )
    parser.add_argument(
        '-c', '--contig',
        help = 'contigs to be analyzed',
        nargs = '*',
        type = str,
        default = chromosomes
    )
    parser.add_argument(
        '-t', '--thread',
        help = 'threads to be used in calculation',
        type = int,
        default = min(os.cpu_count() or 1, 4)
    )
    parser.add_argument(
        '--log-level',
        help = 'logging level (default INFO): [DEBUG, INFO]',
        type = str,
        default = 'INFO'
    )
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {VERSION}'
    )
    args = vars(parser.parse_args(argv))

    loglevel = logging.INFO
    if args['log_level'] == 'DEBUG':
        loglevel = logging.DEBUG
    logging.basicConfig(
        encoding='utf-8',
        format='%(asctime)s %(levelname)s %(message)s',
        level=loglevel
    )
    ##########

    # check input bam file
    message = 'Determine whether the SAM/BAM file has CS or MD tags.'
    logger.debug(message)
    bam_type = check_bam_with_cs_or_md(args['bam_file'])

    if bam_type is None:
        message = 'The SAM/BAM file input should have cs tags or MD tags.'
        logger.error(message)
        raise ValueError(message)
    elif bam_type == 'MD':
        message = 'The SAM/BAM has MD tag.'
        logger.debug(message)
        if args['genome_fasta'] is None:
            message = 'The SAM/BAM file has MD tags but no cs tags, ' +\
                'so the genome fasta file should be provided.'
            logger.error(message)
            raise ValueError(message)
        else:
            message = 'Genome fasta provided.'
            logger.debug(message)
    elif bam_type == 'both':
        message = 'The SAM/BAM has both cs tag and MD tag. Use cs tag.'
        logger.debug(message)
        bam_type = 'cs'
    else:
        message = 'The SAM/BAM has cs tag.'
        logger.debug(message)


    # output path
    o_dirs = {
        'base': args['output_dir'],
        'table': os.path.join(
            args['output_dir'], 'table'
        ),
        'fig': os.path.join(
            args['output_dir'], 'fig'
        ),
    }

    o_files = {
        'pickle': os.path.join(
            o_dirs['base'], 'result.pickle'
        ),
        'cs': os.path.join(
            o_dirs['base'], 'read.cs'
        ),
        't_readstat': os.path.join(
            o_dirs['table'],
            'read_stat.txt'
        ),
        't_insertion': os.path.join(
            o_dirs['table'],
            'insertion.txt'
        ),
        't_deletion': os.path.join(
            o_dirs['table'],
            'deletion.txt'
        ),
        't_mismatch': os.path.join(
            o_dirs['table'],
            'mismatch.txt'
        ),
        't_splice': os.path.join(
            o_dirs['table'],
            'splice.txt'
        ),
        't_mapping': os.path.join(
            o_dirs['table'], 'mapping.txt'
        ),
        't_splice_all': os.path.join(
            o_dirs['table'], 'splice_all.txt'
        ),
        'html': os.path.join(
            o_dirs['base'], 'LQC_report.html'
        )
    }

    # build output directory structure
    build_directories(o_dirs)

    # keep only requested contigs that actually exist in the input; fetching a
    # contig absent from the BAM header raises ValueError in pysam.
    bam_contigs = set(list_bam_contigs(args['bam_file']))
    contigs = []
    seen = set()
    for contig in args['contig']:
        if contig in seen:
            continue
        seen.add(contig)
        if contig in bam_contigs:
            contigs.append(contig)
        else:
            logger.warning('Contig %s is not in the BAM; skipped.', contig)
    if not contigs:
        message = 'None of the requested contigs are present in the BAM.'
        logger.error(message)
        raise ValueError(message)

    if args['output_cs']:
        for stale in os.listdir(o_dirs['base']):
            if stale.startswith('.readcs-') and stale.endswith('.tmp'):
                os.remove(os.path.join(o_dirs['base'], stale))
        omitted = bam_contigs - set(contigs)
        if omitted:
            logger.warning(
                'read.cs contains only the analyzed contigs; %d other header '
                'reference(s) are omitted: %s',
                len(omitted), ', '.join(sorted(omitted))
            )

    # run jobs by balanced per-contig coordinate windows (worker-side fetch)
    message = 'Element statistic process starts.'
    logger.info(message)
    tasks = plan_tasks(args['bam_file'], contigs)

    with mp.Pool(args['thread']) as p:
        block_results = p.map(
            partial(
                stat_region,
                bam_file = args['bam_file'],
                genome_file = args['genome_fasta'],
                method = bam_type,
                cs_dir = o_dirs['base'] if args['output_cs'] else None
            ),
            tasks
        )

    result = reduce_blocks_to_contigs(block_results, contigs)
    message = 'Element statistic process finished.'
    logger.info(message)

    # drop contigs that are in the header but have no mapped reads; an empty
    # ReadStat has undefined mean/median/N50 and would break table/plots.
    result = [
        r for r in result
        if r.readstat.get_read_count() > 0
    ]
    if not result:
        message = 'No reads were found for the requested contigs.'
        logger.error(message)
        raise ValueError(message)

    l_readstat = [cs.readstat for cs in result]
    l_insertion = [cs.insertion for cs in result]
    l_deletion = [cs.deletion for cs in result]
    l_mismatch = [cs.mismatch for cs in result]
    l_splice = [cs.splice for cs in result]

    message = 'Sum of statistics from each contig.'
    logger.debug(message)
    sreadstat = concat_stats(l_readstat, TOTAL_LABEL)
    sinsertion = concat_stats(l_insertion, TOTAL_LABEL)
    sdeletion = concat_stats(l_deletion, TOTAL_LABEL)
    smismatch = concat_stats(l_mismatch, TOTAL_LABEL)
    ssplice = concat_stats(l_splice, TOTAL_LABEL)

    data_sets = {
        'readstat': (l_readstat, sreadstat),
        'insertion': (l_insertion, sinsertion),
        'deletion': (l_deletion, sdeletion),
        'mismatch': (l_mismatch, smismatch),
        'splice': (l_splice, ssplice),
    }

    message = 'Generate summary tables.'
    logger.info(message)
    t_readstat = create_readstat_table(l_readstat, sreadstat)
    t_insertion = create_indel_summary_table(l_insertion, sinsertion)
    t_deletion = create_indel_summary_table(l_deletion, sdeletion)
    t_mismatch = create_mismatch_normalized_read_location_table(l_mismatch, smismatch)
    t_splice = create_splice_table(l_splice, ssplice)
    t_mapping = create_mapping_table(l_readstat, sreadstat)
    t_splice_all = create_splice_all_table(l_splice, ssplice)

    ####################
    # output

    # write a pickle for results
    if args['output_pickle']:
        message = 'Output pickle file.'
        logger.info(message)
        outdict = {
            'readstat_contig': l_readstat,
            'readstat_sum': sreadstat,
            'insertion_contig': l_insertion,
            'insertion_sum': sinsertion,
            'deletion_contig': l_deletion,
            'deletion_sum': sdeletion,
            'mismatch_contig': l_mismatch,
            'mismatch_sum': smismatch,
            'splice_contig': l_splice,
            'splice_sum': ssplice
        }
        with open(o_files['pickle'], 'wb') as f:
            pickle.dump(outdict, f)

        message = 'Output pickle file finished.'
        logger.debug(message)
    else:
        pass

    if args['output_cs']:
        message = 'Output processed cs tags.'
        logger.info(message)
        tmp_cs = o_files['cs'] + '.tmp'
        with open(tmp_cs, 'w') as out:
            out.write(
                'read_name\tcontig\tlow\thigh\tcs_mark\tcs_value\n'
            )
            for block in block_results:
                try:
                    with open(block.cs_path, 'r') as fh:
                        out.write(fh.read())
                finally:
                    os.remove(block.cs_path)
        os.replace(tmp_cs, o_files['cs'])
        message = 'Output processed cs tags finished.'
        logger.debug(message)
    else:
        pass

    ####################
    # write output tables
    message = 'Output summary tables.'
    logger.info(message)
    for key, table in [
            ('t_readstat', t_readstat),
            ('t_insertion', t_insertion),
            ('t_deletion', t_deletion),
            ('t_mismatch', t_mismatch),
            ('t_splice', t_splice),
            ('t_mapping', t_mapping),
            ('t_splice_all', t_splice_all),
    ]:
        table.to_csv(
            o_files[key], sep = '\t', index = False,
            float_format = FLOAT_FORMAT
        )
    message = 'Output summary tables finished.'
    logger.debug(message)

    ####################
    # plot figures
    message = 'Output figures.'
    logger.info(message)
    figdir = o_dirs['fig']
    for feature, stem in READSTAT_BAR_SPECS:
        fig = plot_readstat_bar(l_readstat, feature)
        savefig(fig, os.path.join(figdir, stem))
        plt.close('all')

    for stem, plot_func, data_key in MULTI_FIG_SPECS:
        data_list, data_sum = data_sets[data_key]
        generate_multiple_figs(
            plot_func, data_list, data_sum,
            os.path.join(figdir, stem), width = 5, height = 4
        )

    for stem, data_key, kind in ELEMENT_BAR_SPECS:
        fig = plot_element_total_count(data_sets[data_key][0], kind)
        savefig(fig, os.path.join(figdir, stem))
        plt.close('all')

    message = 'Output figures finished.'
    logger.debug(message)

    ####################
    # generate report
    message = 'Output html report page.'
    logger.info(message)

    copy_logo(o_dirs['fig'])

    html_string = get_html_template()

    new_html_string = html_add_readstat_table(
        html_string, t_readstat
    )

    new_html_string = html_add_mapping(
        new_html_string, t_mapping
    )

    mismatch_type_counter = smismatch.get_type_count()

    new_html_string = html_add_mismatch_table(
        new_html_string, t_mismatch,
        smismatch.get_total_count(),
        total_row(t_readstat, 'mean_mismatch_per_read_per_kb'),
        mismatch_type_counter
    )

    new_html_string = html_add_insertion_table(
        new_html_string, t_insertion,
        total_row(t_readstat, 'mean_insertion_per_read_per_kb')
    )

    new_html_string = html_add_deletion_table(
        new_html_string, t_deletion,
        total_row(t_readstat, 'mean_deletion_per_read_per_kb')
    )

    new_html_string = html_add_splice_table(
        new_html_string, t_splice,
        total_row(t_readstat, 'mean_intron_per_read')
    )

    new_html_string = html_add_bootstrap(new_html_string, VERSION)

    new_html_string = html_add_data(
        new_html_string,
        {
            'readstat': t_readstat,
            'mapping': t_mapping,
            'insertion': t_insertion,
            'deletion': t_deletion,
            'mismatch': t_mismatch,
            'splice': t_splice,
        }
    )

    new_html_string = inline_figures(new_html_string, o_dirs['fig'])

    with open(o_files['html'], 'w', encoding='utf-8') as f:
        f.write(new_html_string)

    message = 'All done!'
    logger.info(message)
    return 0


if __name__ == '__main__':
    sys.exit(main())


########################################
