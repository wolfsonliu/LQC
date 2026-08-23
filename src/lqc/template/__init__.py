import os
import shutil

file_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(file_dir)


def get_html_template():
    with open(
        os.path.join(file_dir, "template.html"),
        "r"
    ) as html:
        html_string = "\n".join(
            [read.strip() for read in html]
        )
    return html_string


def copy_logo(aim_dir):
    logo_path = os.path.join(file_dir, "logo_white.svg")
    shutil.copy(logo_path, aim_dir)
    card_path = os.path.join(file_dir, "card.svg")
    shutil.copy(card_path, aim_dir)

########################################
