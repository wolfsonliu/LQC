"""Command-line interface for LQC.

All orchestration lives in main(). See README.md for usage.
"""

import argparse
import logging
import multiprocessing as mp
import os
import pickle
import sys
from copy import copy
from functools import partial

import matplotlib
import matplotlib.pyplot as plt

from lqc import __version__ as VERSION
from lqc import (
    check_bam_with_cs_or_md,
    copy_logo,
    create_indel_summary_table,
    create_mismatch_normalized_read_location_table,
    create_readstat_table,
    create_splice_table,
    get_html_template,
    html_add_bootstrap,
    html_add_data,
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mismatch_table,
    html_add_readstat_table,
    html_add_splice_table,
    inline_figures,
    list_bam_contigs,
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
from lqc.stat import prefetch_records, reduce_blocks_to_contigs, stat_records

matplotlib.use('Agg')

logger = logging.getLogger(__name__)


def get_stat_list(result, stat_type):
    stat_type_ids = {
        'readstat': 0,
        'insertion': 1,
        'deletion': 2,
        'mismatch': 3,
        'splice': 4
    }
    assert stat_type in stat_type_ids,\
        'stat_type should be in [{}].'.format(
            ','.join(
                stat_type_ids.keys()
            )
        )
    idx = stat_type_ids[stat_type]

    return [result[i][idx] for i in range(len(result))]


def _chunk(records, n):
    n = max(1, int(n))
    if not records:
        return []
    size = (len(records) + n - 1) // n
    return [records[i:i + size] for i in range(0, len(records), size)]


def build_directories(dir_dict):
    for b in dir_dict.values():
        os.makedirs(b, exist_ok = True)


def savefig(fig, prefix):
    fig.savefig(prefix + '.png')
    fig.savefig(prefix + '.pdf')


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
        'f_readstat_bar_Read count': os.path.join(
            o_dirs['fig'], 'readstat_bar_Read_count'
        ),
        'f_readstat_bar_Median read length': os.path.join(
            o_dirs['fig'], 'readstat_bar_Median_read_length'
        ),
        'f_readstat_bar_Mean read length': os.path.join(
            o_dirs['fig'], 'readstat_bar_Mean_read_length'
        ),
        'f_readstat_bar_Insertions per read': os.path.join(
            o_dirs['fig'], 'readstat_bar_insertions_per_read'
        ),
        'f_readstat_bar_Insertions per read per kb': os.path.join(
            o_dirs['fig'], 'readstat_bar_insertions_per_read_per_kb'
        ),
        'f_readstat_bar_Deletions per read': os.path.join(
            o_dirs['fig'], 'readstat_bar_deletions_per_read'
        ),
        'f_readstat_bar_Deletions per read per kb': os.path.join(
            o_dirs['fig'], 'readstat_bar_deletions_per_read_per_kb'
        ),
        'f_readstat_bar_Mismatches per read': os.path.join(
            o_dirs['fig'], 'readstat_bar_mismatches_per_read'
        ),
        'f_readstat_bar_Mismatches per read per kb': os.path.join(
            o_dirs['fig'], 'readstat_bar_mismatches_per_read_per_kb'
        ),
        'f_readstat_bar_Mean intron number': os.path.join(
            o_dirs['fig'], 'readstat_bar_Mean_intron_number'
        ),
        'f_readstat_bar_N50': os.path.join(
            o_dirs['fig'], 'readstat_bar_N50'
        ),
        'f_readstat_bar_L50': os.path.join(
            o_dirs['fig'], 'readstat_bar_L50'
        ),
        'f_readstat_bar_mean_element_per_read': os.path.join(
            o_dirs['fig'], 'readstat_bar_mean_element_per_read'
        ),
        'f_readstat_bar_mean_element_per_read_per_kb': os.path.join(
            o_dirs['fig'], 'readstat_bar_mean_element_per_read_per_kb'
        ),
        'f_readstat_line_cumulative_length': os.path.join(
            o_dirs['fig'], 'readstat_line_cumulative_length'
        ),
        'f_readstat_bar_ratio_with_element': os.path.join(
            o_dirs['fig'], 'readstat_bar_ratio_with_element'
        ),
        'f_readstat_hist_length': os.path.join(
            o_dirs['fig'], 'readstat_hist_length'
        ),
        'f_element_Insertion_bar_count': os.path.join(
            o_dirs['fig'], 'insertion_bar_count'
        ),
        'f_element_Deletion_bar_count': os.path.join(
            o_dirs['fig'], 'deletion_bar_count'
        ),
        'f_element_Mismatch_bar_count': os.path.join(
            o_dirs['fig'], 'mismatch_bar_count'
        ),
        'f_element_Intron_bar_count': os.path.join(
            o_dirs['fig'], 'intron_bar_count'
        ),
        'f_insertion_hist_length': os.path.join(
            o_dirs['fig'], 'insertion_hist_length'
        ),
        'f_deletion_hist_length': os.path.join(
            o_dirs['fig'], 'deletion_hist_length'
        ),
        'f_insertion_hist_location': os.path.join(
            o_dirs['fig'], 'insertion_hist_location'
        ),
        'f_deletion_hist_location': os.path.join(
            o_dirs['fig'], 'deletion_hist_location'
        ),
        'f_mismatch_type': os.path.join(
            o_dirs['fig'], 'mismatch_type'
        ),
        'f_mismatch_hist_location': os.path.join(
            o_dirs['fig'], 'mismatch_hist_location'
        ),
        'f_splice_type': os.path.join(
            o_dirs['fig'], 'splice_type'
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

    # run jobs by chunked reads (data-parallel, independent of contig count)
    message = 'Element statistic process starts.'
    logger.info(message)
    per_contig = prefetch_records(
        args['bam_file'], contigs, bam_type
    )
    tasks = []
    for contig, records in per_contig:
        if not records:
            continue
        for chunk in _chunk(records, args['thread']):
            tasks.append((len(tasks), contig, chunk))

    with mp.Pool(args['thread']) as p:
        block_results = p.map(
            partial(
                stat_records,
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
        if r[0].get_read_count() > 0
    ]
    if not result:
        message = 'No reads were found for the requested contigs.'
        logger.error(message)
        raise ValueError(message)

    l_readstat = get_stat_list(result, 'readstat')
    l_insertion = get_stat_list(result, 'insertion')
    l_deletion = get_stat_list(result, 'deletion')
    l_mismatch = get_stat_list(result, 'mismatch')
    l_splice = get_stat_list(result, 'splice')

    message = 'Sum of statistics from each contig.'
    logger.debug(message)
    # The summed object is a shallow relabeled copy: ``sum([x])`` returns ``x``
    # itself, so relabeling in place would rename the per-contig row too. A
    # shallow copy shares the (read-only after this point) data while carrying
    # its own ``label``, which also avoids duplicating multi-GB event lists.
    def relabel(obj):
        new_obj = copy(obj)
        new_obj.label = 'Total'
        return new_obj
    sreadstat = relabel(sum(l_readstat))
    sinsertion = relabel(sum(l_insertion))
    sdeletion = relabel(sum(l_deletion))
    smismatch = relabel(sum(l_mismatch))
    ssplice = relabel(sum(l_splice))

    message = 'Generate summary tables.'
    logger.info(message)
    t_readstat = create_readstat_table(l_readstat, sreadstat)
    t_insertion = create_indel_summary_table(l_insertion, sinsertion)
    t_deletion = create_indel_summary_table(l_deletion, sdeletion)
    t_mismatch = create_mismatch_normalized_read_location_table(l_mismatch, smismatch)
    t_splice = create_splice_table(l_splice, ssplice)

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
    t_readstat.to_csv(
        o_files['t_readstat'],
        sep = '\t', index = False
    )
    t_insertion.to_csv(
        o_files['t_insertion'],
        sep = '\t', index = False
    )
    t_deletion.to_csv(
        o_files['t_deletion'],
        sep = '\t', index = False
    )
    t_mismatch.to_csv(
        o_files['t_mismatch'],
        sep = '\t', index = False
    )
    t_splice.to_csv(
        o_files['t_splice'],
        sep = '\t', index = False
    )
    message = 'Output summary tables finished.'
    logger.debug(message)

    ####################
    # plot figures
    message = 'Output figures.'
    logger.info(message)
    # readstat: feature
    message = 'Output figures: readstat summary features.'
    logger.debug(message)
    for feature in [
            'Read count', 'Median read length',
            'Mean read length',
            'Insertions per read', 'Insertions per read per kb',
            'Deletions per read', 'Deletions per read per kb',
            'Mismatches per read', 'Mismatches per read per kb',
            'Mean intron number', 'N50', 'L50'
    ]:
        filelabel = 'f_readstat_bar_' + feature
        fig = plot_readstat_bar(l_readstat, feature)
        savefig(fig, o_files[filelabel])
        plt.close('all')

    # readstat: mean element per read
    message = 'Output figures: readstat_bar_mean_element_per_read.'
    logger.debug(message)
    filelabel = 'f_readstat_bar_mean_element_per_read'
    generate_multiple_figs(
        plot_readstat_bar_mean_element_per_read,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
    message = 'Output figures: readstat_bar_mean_element_per_read_per_kb.'
    logger.debug(message)
    filelabel = 'f_readstat_bar_mean_element_per_read_per_kb'
    generate_multiple_figs(
        plot_readstat_bar_mean_element_per_read_per_kb,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # readstat: cumulative length
    message = 'Output figures: readstat_line_cumulative_length.'
    logger.debug(message)
    filelabel = 'f_readstat_line_cumulative_length'
    generate_multiple_figs(
        plot_readstat_cumulative_length,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # readstat: ratio with error element
    message = 'Output figures: readstat_bar_ratio_with_element.'
    logger.debug(message)
    filelabel = 'f_readstat_bar_ratio_with_element'
    generate_multiple_figs(
        plot_readstat_bar_ratio_with_element,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # readstat: length hist
    message = 'Output figures: readstat_hist_length.'
    logger.debug(message)
    filelabel = 'f_readstat_hist_length'
    generate_multiple_figs(
        plot_readstat_length_hist,
        data_list = l_readstat,
        data_sum = sreadstat,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # error element
    message = 'Output figures: error element barplots.'
    logger.debug(message)
    filelabel = 'f_element_Insertion_bar_count'
    fig = plot_element_total_count(
        l_insertion, 'Insertion'
    )
    savefig(fig, o_files[filelabel])
    filelabel = 'f_element_Deletion_bar_count'
    fig = plot_element_total_count(
        l_deletion, 'Deletion'
    )
    savefig(fig, o_files[filelabel])
    filelabel = 'f_element_Mismatch_bar_count'
    fig = plot_element_total_count(
        l_mismatch, 'Mismatch'
    )
    savefig(fig, o_files[filelabel])
    plt.close('all')
    filelabel = 'f_element_Intron_bar_count'
    fig = plot_element_total_count(
        l_splice, 'Intron'
    )
    savefig(fig, o_files[filelabel])
    plt.close('all')

    # indel: length hist
    # insertion
    message = 'Output figures: insertion_hist_length.'
    logger.debug(message)
    filelabel = 'f_insertion_hist_length'
    generate_multiple_figs(
        plot_indel_hist_length,
        data_list = l_insertion,
        data_sum = sinsertion,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
    # deletion
    message = 'Output figures: deletion_hist_length.'
    logger.debug(message)
    filelabel = 'f_deletion_hist_length'
    generate_multiple_figs(
        plot_indel_hist_length,
        data_list = l_deletion,
        data_sum = sdeletion,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # indel: read location hist
    # insertion
    message = 'Output figures: insertion_hist_location.'
    logger.debug(message)
    filelabel = 'f_insertion_hist_location'
    generate_multiple_figs(
        plot_indel_hist_location,
        data_list = l_insertion,
        data_sum = sinsertion,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )
    # deletion
    message = 'Output figures: deletion_hist_location.'
    logger.debug(message)
    filelabel = 'f_deletion_hist_location'
    generate_multiple_figs(
        plot_indel_hist_location,
        data_list = l_deletion,
        data_sum = sdeletion,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # mismatch type
    message = 'Output figures: mismatch_type.'
    logger.debug(message)
    filelabel = 'f_mismatch_type'
    generate_multiple_figs(
        plot_mismatch_type_count,
        data_list = l_mismatch,
        data_sum = smismatch,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # mismatch: read location hist
    message = 'Output figures: mismatch_hist_location.'
    logger.debug(message)
    filelabel = 'f_mismatch_hist_location'
    generate_multiple_figs(
        plot_mismatch_hist_location,
        data_list = l_mismatch,
        data_sum = smismatch,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

    # splice type
    message = 'Output figures: splice_type.'
    logger.debug(message)
    filelabel = 'f_splice_type'
    generate_multiple_figs(
        plot_splice_type_count,
        data_list = l_splice,
        data_sum = ssplice,
        filelabel = o_files[filelabel],
        width = 5, height = 4
    )

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

    mismatch_type_counter = smismatch.get_type_count()

    new_html_string = html_add_mismatch_table(
        new_html_string, t_mismatch,
        smismatch.get_total_count(),
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_mismatch_per_read_per_kb'
        ].values[0],
        mismatch_type_counter
    )

    new_html_string = html_add_insertion_table(
        new_html_string, t_insertion,
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_insertion_per_read_per_kb'
        ].values[0]
    )

    new_html_string = html_add_deletion_table(
        new_html_string, t_deletion,
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_deletion_per_read_per_kb'
        ].values[0]
    )

    new_html_string = html_add_splice_table(
        new_html_string, t_splice,
        t_readstat.loc[
            t_readstat['label'] == 'Total',
            'mean_intron_per_read'
        ].values[0]
    )

    new_html_string = html_add_bootstrap(new_html_string, VERSION)

    new_html_string = html_add_data(
        new_html_string,
        {
            'readstat': t_readstat,
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
