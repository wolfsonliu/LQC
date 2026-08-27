from lqc import constants
from lqc.report_table import MISMATCH_TYPES_IN_ORDER


def test_total_label_sentinel():
    assert constants.TOTAL_LABEL == 'Total'


def test_mismatch_types_order():
    assert constants.MISMATCH_TYPES == (
        'ac', 'ag', 'at', 'ca', 'cg', 'ct',
        'ga', 'gc', 'gt', 'ta', 'tc', 'tg',
    )


def test_report_table_uses_shared_mismatch_types():
    assert tuple(MISMATCH_TYPES_IN_ORDER) == constants.MISMATCH_TYPES