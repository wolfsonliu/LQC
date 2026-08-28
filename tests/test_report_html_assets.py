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
    assert 'Bootstrap v5.3.8' in html
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


def test_inline_figures_replaces_png_and_svg(tmp_path):
    import base64

    from lqc.report_html import inline_figures

    png = b'\x89PNG\r\n\x1a\nfake'
    svg = b'<svg></svg>'
    (tmp_path / 'a.png').write_bytes(png)
    (tmp_path / 'b.svg').write_bytes(svg)
    html = '<img src="fig/a.png" /><img src="fig/b.svg" />'
    out = inline_figures(html, str(tmp_path))
    assert 'src="fig/' not in out
    assert (
        f'src="data:image/png;base64,{base64.b64encode(png).decode("ascii")}"'
        in out
    )
    assert (
        f'src="data:image/svg+xml;base64,{base64.b64encode(svg).decode("ascii")}"'
        in out
    )


def test_inline_figures_raises_on_missing_file(tmp_path):
    import pytest

    from lqc.report_html import inline_figures

    with pytest.raises(FileNotFoundError):
        inline_figures('<img src="fig/missing.png" />', str(tmp_path))


def test_inline_figures_raises_on_unknown_type(tmp_path):
    import pytest

    from lqc.report_html import inline_figures

    (tmp_path / 'x.foo').write_bytes(b'x')
    with pytest.raises(ValueError):
        inline_figures('<img src="fig/x.foo" />', str(tmp_path))