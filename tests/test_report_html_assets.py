from lqc import get_html_template


def test_template_has_asset_tokens_and_no_remote_resources():
    html = get_html_template()
    assert '{%bootstrap_css%}' in html
    assert '{%bootstrap_js%}' in html
    assert '{%lqc_data%}' in html
    assert '{%version%}' in html
    assert 'cdn.jsdelivr.net' not in html
    assert 'img.shields.io' not in html


def test_html_add_bootstrap_inlines_vendored_assets():
    from lqc.report_html import html_add_bootstrap

    html = html_add_bootstrap(get_html_template(), '0.0.5')
    assert '{%bootstrap_css%}' not in html
    assert '{%bootstrap_js%}' not in html
    assert '{%version%}' not in html
    assert '<style>' in html
    assert '</style>' in html
    assert 'Bootstrap v5.1.3' in html
    assert 'LQC v0.0.5' in html


def test_html_add_data_injects_parseable_json():
    import json

    import pandas as pd

    from lqc.report_html import html_add_data

    readstat = pd.DataFrame([{'label': 'Total', 'read_count': 2}])
    html = html_add_data(
        '<html>{%lqc_data%}</html>',
        {'readstat': readstat},
    )
    assert '{%lqc_data%}' not in html
    assert '<script type="application/json" id="lqc-data">' in html
    start = html.index('id="lqc-data">') + len('id="lqc-data">')
    end = html.index('</script>', start)
    payload = json.loads(html[start:end])
    assert payload == {'readstat': [{'label': 'Total', 'read_count': 2}]}


def test_html_add_data_escapes_script_terminator():
    import pandas as pd

    from lqc.report_html import html_add_data

    df = pd.DataFrame([{'label': '</script>'}])
    html = html_add_data('<html>{%lqc_data%}</html>', {'block': df})
    inner = html.split('id="lqc-data">')[1].split('</script>')[0]
    assert '</script>' not in inner
    assert '<\\/script>' in inner