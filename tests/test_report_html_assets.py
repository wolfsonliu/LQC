from lqc import get_html_template


def test_template_has_asset_tokens_and_no_remote_resources():
    html = get_html_template()
    assert '{%bootstrap_css%}' in html
    assert '{%bootstrap_js%}' in html
    assert '{%lqc_data%}' in html
    assert '{%version%}' in html
    assert 'cdn.jsdelivr.net' not in html
    assert 'img.shields.io' not in html