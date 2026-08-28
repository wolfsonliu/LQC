import base64
import json
import os
import re

from lqc.constants import MISMATCH_TYPES, TOTAL_LABEL


def _replace_tokens(html_string, tokens):
    """Apply each ``{token: replacement}`` via ``re.sub`` in insertion order."""
    result = html_string
    for token, replacement in tokens.items():
        result = re.sub(token, replacement, result)
    return result


def html_add_readstat_table(html_string, readstat_table):
    rowstring_list = []
    for _, row in readstat_table.iterrows():
        tmprow_list = [
            '<th scope="row">{}</th>'.format(row['label']),
            '<td>{}</td>'.format(row['read_count']),
            '<td>{}</td>'.format(row['total_base']),
            '<td>{:.4g}</td>'.format(row['read_length_mean']),
            '<td>{}</td>'.format(row['read_length_median']),
            '<td>{}</td>'.format(row['read_length_N50']),
            '<td>{}</td>'.format(row['read_length_L50'])
        ]
        if row['label'] == "Total":
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + "</tr>"

            )
            total_read_count = row['read_count']
            total_base = row['total_base']
            total_median_read_length = row['read_length_median']
            total_mean_read_length = row['read_length_mean']
            total_N50 = row['read_length_N50']
            total_L50 = row['read_length_L50']
            total_insertion_per_read_per_kb = row['mean_insertion_per_read_per_kb']
            total_deletion_per_read_per_kb = row['mean_deletion_per_read_per_kb']
            total_mismatch_per_read_per_kb = row['mean_mismatch_per_read_per_kb']
        else:
            rowstring_list.append(
                "<tr>" +
                '\n'.join(tmprow_list) + "</tr>"
            )

    new_html_string = _replace_tokens(html_string, {
        r"\{%readstat_total_read_count%\}": f'{total_read_count}',
        r"\{%readstat_total_base%\}": f'{total_base}',
        r"\{%readstat_median_read_length%\}": f'{total_median_read_length}',
        r"\{%readstat_mean_read_length%\}": f'{total_mean_read_length:.4g}',
        r"\{%readstat_N50%\}": f'{total_N50}',
        r"\{%readstat_L50%\}": f'{total_L50}',
        r"\{%readstat_insertion_per_read_per_kb%\}" :
            f'{total_insertion_per_read_per_kb:.4g}',
        r"\{%readstat_deletion_per_read_per_kb%\}" :
            f'{total_deletion_per_read_per_kb:.4g}',
        r"\{%readstat_mismatch_per_read_per_kb%\}" :
            f'{total_mismatch_per_read_per_kb:.4g}',
        r"\{%readstat_table%\}": '\n'.join(rowstring_list),
    })
    return new_html_string


def html_add_mismatch_table(html_string, mismatch_table,
                            total_mismatch_count,
                            mean_mismatch_per_read_per_kb,
                            mismatch_type_counter):
    mis_types = list(MISMATCH_TYPES)
    rowstring_list = []
    for _, row in mismatch_table.iterrows():
        for mis in mis_types:
            tmprow_list = [
                    '<th scope="row">{}</th>'.format(row['label'])
            ]
            tmprow_list.append(
                '<td>{}</td>'.format(row['bin'])
            )
            tmprow_list.append(
                f'<td>{mis}</td>'
            )
            if mis in row:
                tmprow_list.append(
                    f'<td>{row[mis]}</td>'
                )
            else:
                tmprow_list.append('<td>0</td>')
            if row['label'] == "Total":
                rowstring_list.append(
                    '<tr class="table-secondary">' +
                    '\n'.join(tmprow_list) + "</tr>"

                )
            else:
                rowstring_list.append(
                    "<tr>" +
                    '\n'.join(tmprow_list) + "</tr>"
                )

    mistype_list1 = [
        f"<li>{mis}: {mismatch_type_counter[mis]}</li>"
        for mis in mis_types[0:6]
    ]
    mistype_list2 = [
        f"<li>{mis}: {mismatch_type_counter[mis]}</li>"
        for mis in mis_types[6:]
    ]
    new_html_string = _replace_tokens(html_string, {
        r"\{%mismatch_total_mismatch_number%\}": f'{total_mismatch_count}',
        r"\{%mismatch_mean_mismatch_per_read_per_kb%\}" :
            f'{mean_mismatch_per_read_per_kb:.4g}',
        r"\{%mismatch_type_list1%\}" :
            "<ul>" + '\n'.join(mistype_list1) + "</ul>",
        r"\{%mismatch_type_list2%\}" :
            "<ul>" + '\n'.join(mistype_list2) + "</ul>",
        r"\{%mismatch_table%\}": '\n'.join(rowstring_list),
    })
    return new_html_string


def _html_add_indel_table(html_string, indel_table, per_kb, kind):
    token = kind  # 'insertion' or 'deletion'
    rowstring_list = []
    total_count = total_length = 0
    mean_length = median_length = 0.0
    for _, row in indel_table.iterrows():
        tmprow_list = [
            '<th scope="row">{}</th>'.format(row['label']),
            '<td>{}</td>'.format(row['total_count']),
            '<td>{}</td>'.format(row['total_length']),
            '<td>{:.4g}</td>'.format(float(row['mean_length'])),
            '<td>{:.4g}</td>'.format(float(row['median_length']))
        ]
        if row['label'] == TOTAL_LABEL:
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + "</tr>"
            )
            total_count = row['total_count']
            total_length = row['total_length']
            mean_length = float(row['mean_length'])
            median_length = float(row['median_length'])
        else:
            rowstring_list.append(
                "<tr>" + '\n'.join(tmprow_list) + "</tr>"
            )

    return _replace_tokens(html_string, {
        r"\{%" + token + r"_total_" + token + r"_number%\}": f'{total_count}',
        r"\{%" + token + r"_total_" + token + r"_length%\}": f'{total_length}',
        r"\{%" + token + r"_mean_length%\}": f'{mean_length:.4g}',
        r"\{%" + token + r"_median_length%\}": f'{median_length:.4g}',
        r"\{%" + token + r"_mean_" + token + r"_per_read_per_kb%\}":
            f'{per_kb:.4g}',
        r"\{%" + token + r"_table%\}": '\n'.join(rowstring_list),
    })


def html_add_insertion_table(html_string, insertion_table,
                             mean_insertion_per_read_per_kb):
    return _html_add_indel_table(
        html_string, insertion_table, mean_insertion_per_read_per_kb,
        'insertion'
    )


def html_add_deletion_table(html_string, deletion_table,
                            mean_deletion_per_read_per_kb):
    return _html_add_indel_table(
        html_string, deletion_table, mean_deletion_per_read_per_kb,
        'deletion'
    )


def html_add_splice_table(html_string, splice_table,
                          mean_intron_per_read):

    rowstring_list = []
    total_gtag = 0
    total_gtagp = 0.0
    total_gcag = 0
    total_gcagp = 0.0
    total_atac = 0
    total_atacp = 0.0
    total_other = 0
    total_otherp = 0.0
    for _, row in splice_table.iterrows():
        gtag = row['gt-ag']
        gtagp = row['gt-ag_pct']
        gcag = row['gc-ag']
        gcagp = row['gc-ag_pct']
        atac = row['at-ac']
        atacp = row['at-ac_pct']
        other = row['other']
        otherp = row['other_pct']
        tmprow_list = [
            f'<th scope="row">{row["label"]}</th>',
            f'<td>{gtag}</td>',
            f'<td>{gtagp:.4g}</td>',
            f'<td>{gcag}</td>',
            f'<td>{gcagp:.4g}</td>',
            f'<td>{atac}</td>',
            f'<td>{atacp:.4g}</td>',
            f'<td>{other}</td>',
            f'<td>{otherp:.4g}</td>'
        ]
        if row['label'] == "Total":
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + "</tr>"
            )
            total_gtag = gtag
            total_gtagp = gtagp
            total_gcag = gcag
            total_gcagp = gcagp
            total_atac = atac
            total_atacp = atacp
            total_other = other
            total_otherp = otherp
        else:
            rowstring_list.append(
                "<tr>" +
                '\n'.join(tmprow_list) + "</tr>"
            )

    splice_list = [
        f"<li>gt-ag: {total_gtag} ({total_gtagp:.4g}%)</li>",
        f"<li>gc-ag: {total_gcag} ({total_gcagp:.4g}%)</li>",
        f"<li>at-ac: {total_atac} ({total_atacp:.4g}%)</li>",
        f"<li>other: {total_other} ({total_otherp:.4g}%)</li>"
    ]
    new_html_string = _replace_tokens(html_string, {
        r"\{%splice_type_list%\}" :
            "<ul>" + '\n'.join(splice_list) + "</ul>",
        r"\{%intron_total_intron_number%\}" :
            f'{total_gtag + total_gcag + total_atac + total_other}',
        r"\{%intron_mean_intron_per_read%\}" : f'{mean_intron_per_read:.4g}',
        r"\{%splice_table%\}": '\n'.join(rowstring_list),
    })
    return new_html_string


def html_add_mapping(html_string, mapping_table):
    rowstring_list = []
    total_aligned_fraction_mean = 0.0
    total_aligned_fraction_median = 0.0
    total_mapq_mean = 0.0
    total_mapq_median = 0.0
    for _, row in mapping_table.iterrows():
        tmprow_list = [
            '<th scope="row">{}</th>'.format(row['label']),
            '<td>{}</td>'.format(row['read_count']),
            '<td>{}</td>'.format(row['query_base']),
            '<td>{}</td>'.format(row['aligned_base']),
            '<td>{:.4g}</td>'.format(float(row['aligned_fraction_mean'])),
            '<td>{:.4g}</td>'.format(float(row['aligned_fraction_median'])),
            '<td>{:.4g}</td>'.format(float(row['mapq_mean'])),
            '<td>{:.4g}</td>'.format(float(row['mapq_median'])),
            '<td>{}</td>'.format(row['reads_aligned_fraction_lt_0.9']),
            '<td>{}</td>'.format(row['reads_fully_aligned'])
        ]
        if row['label'] == 'Total':
            rowstring_list.append(
                '<tr class="table-secondary">' +
                '\n'.join(tmprow_list) + '</tr>'
            )
            total_aligned_fraction_mean = float(row['aligned_fraction_mean'])
            total_aligned_fraction_median = float(row['aligned_fraction_median'])
            total_mapq_mean = float(row['mapq_mean'])
            total_mapq_median = float(row['mapq_median'])
        else:
            rowstring_list.append('<tr>' + '\n'.join(tmprow_list) + '</tr>')

    new_html_string = _replace_tokens(html_string, {
        r"\{%mapping_aligned_fraction_mean%\}" :
            f'{total_aligned_fraction_mean:.4g}',
        r"\{%mapping_aligned_fraction_median%\}" :
            f'{total_aligned_fraction_median:.4g}',
        r"\{%mapping_mapq_mean%\}" : f'{total_mapq_mean:.4g}',
        r"\{%mapping_mapq_median%\}" : f'{total_mapq_median:.4g}',
        r"\{%mapping_table%\}": '\n'.join(rowstring_list),
    })
    return new_html_string


########################################


_BOOTSTRAP_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'template'
)


def _read_text(path):
    with open(path, 'r', encoding='utf-8') as fh:
        return fh.read()


def html_add_bootstrap(html_string, version):
    css = _read_text(os.path.join(_BOOTSTRAP_DIR, 'bootstrap.min.css'))
    # Escape the sequence that would terminate the inline <script> early.
    js = _read_text(
        os.path.join(_BOOTSTRAP_DIR, 'bootstrap.min.js')
    ).replace('</script>', '<\\/script>')

    new_html = re.sub(
        r'\{%bootstrap_css%\}',
        lambda _match: '<style>' + css + '</style>',
        html_string
    )
    new_html = re.sub(
        r'\{%bootstrap_js%\}',
        lambda _match: '<script>' + js + '</script>',
        new_html
    )
    new_html = re.sub(
        r'\{%version%\}',
        lambda _match: version,
        new_html
    )
    return new_html


def html_add_data(html_string, tables):
    # to_json -> json.loads round-trip converts numpy scalars to native JSON
    # types so json.dumps never hits a non-serializable int64/float64.
    payload = {
        name: json.loads(table.to_json(orient='records'))
        for name, table in tables.items()
    }
    raw = json.dumps(payload).replace('</', '<\\/')
    script = (
        '<script type="application/json" id="lqc-data">'
        + raw
        + '</script>'
    )
    # Callable replacement: ``script`` can contain backslashes (the ``<\/``
    # escaping above), which a string re.sub replacement would mis-parse.
    return re.sub(r'\{%lqc_data%\}', lambda _match: script, html_string)


_MIME_BY_EXT = {
    '.png': 'image/png',
    '.svg': 'image/svg+xml',
    '.jpg': 'image/jpeg',
    '.jpeg': 'image/jpeg',
    '.gif': 'image/gif',
}


def inline_figures(html_string, fig_dir):
    def _to_data_uri(match):
        rel = match.group(1)
        ext = os.path.splitext(rel)[1].lower()
        if ext not in _MIME_BY_EXT:
            raise ValueError(f'Unsupported figure type: {rel}')
        path = os.path.join(fig_dir, rel)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f'Figure referenced by HTML is missing: {path}'
            )
        with open(path, 'rb') as fh:
            b64 = base64.b64encode(fh.read()).decode('ascii')
        return f'src="data:{_MIME_BY_EXT[ext]};base64,{b64}"'

    return re.sub(r'src="fig/([^"]+)"', _to_data_uri, html_string)
