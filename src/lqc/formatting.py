"""Shared output-formatting constants."""


# Significant digits for float columns in the TSV summary tables. Uses the same
# 4-significant-digit precision as the ``{:.4}`` formatting in report_html.py,
# though the numeric spelling is not byte-identical to that format.
FLOAT_FORMAT = '%.4g'