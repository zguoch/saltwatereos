# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sys, os
sys.path.append(os.path.abspath('_extensions'))

# -- Project information -----------------------------------------------------

project = 'SaltWater EOS'
latex_name='SaltWaterEOS'
copyright = "Zhikui Guo, Lars Ruepke"
author = 'Zhikui Guo'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax', 
            'jinja',
            'sphinx.ext.ifconfig',
            'sphinx_inline_tabs',
            'sphinxcontrib.bibtex']
source_encoding = 'utf-8-sig'
source_suffix = '.rst'
master_doc = 'index'
templates_path = ['_templates']

# internationalization
language = 'en'
locale_dirs = ['locale/']
gettext_compact = True
gettext_auto_build=True
# Set smartquotes_action to 'qe' to disable Smart Quotes transform of -- and ---
smartquotes_action = 'qe'
# customize OpenFOAM syntax highlight
from sphinx.highlighting import lexers
from pygments_OpenFOAM.foam import OpenFOAMLexer
lexers['foam'] = OpenFOAMLexer(startinline=True)
# default language to highlight source code
# highlight_language = 'foam'
# default language to highlight source code
highlight_language = 'cpp'
pygments_style = 'monokai'
bibtex_bibfiles = ['manual.bib'] 
# -- Project configuration ------------------------------------------------

# The version shown at the top of the sidebar
version = '1.0'
# The full version shown in the page title
release = '1.0'
# Make the "Edit on GitHub" button link to the correct branch
# Default to master branch if Azure Pipelines environmental variable BUILD_SOURCEBRANCHNAME is not defined
# github_version = os.getenv("BUILD_SOURCEBRANCHNAME", 'master')


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'rtd'
html_theme_path = ["_themes"]
html_theme_options = {
    'sticky_navigation': False,
    'includehidden': False,
    'logo_only' : True,
    'sticky_navigation': True,
    'titles_only': True,
    'display_version': False,
    'prev_next_buttons_location': 'both',
    'style_nav_header_background': 'skyblue',
    # 'gitlab_url': 'https://gitlab.com/gmdpapers/hydrothermalfoam'
}

html_context = {
    "menu_links": [
        (
            '<i class="fa fa-book fa-fw"></i> License',
            "xxx",
        ),
        (
            '<i class="fa fa-comment fa-fw"></i> Contact',
            "xxx",
        ),
        (
            '<i class="fa fa-github fa-fw"></i> Source Code',
            "https://gitlab.com/hydrothermal-openfoam/saltwatereos",
        ),
    ],
    'project':project,
    'downloads_url':'hydrothermal-openfoam.gitlab.io/saltwatereos/manual/downloads',
    'latex_main':  latex_name, 
    'pdf_versions': [
        (
            'latest',
            'hydrothermal-openfoam.gitlab.io/saltwatereos/manual/downloads'
        ),
        (
            '2.0',
            '#'
        ),
        ]
}

# favicon of the docs
html_favicon = "_static/logo.png"
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'
# If true, links to the reST sources are added to the pages.
html_logo = "_static/logo.png" 
html_show_sourcelink = False 
# List of custom CSS files (needs sphinx>=1.8)
html_css_files = ["style.css"]

# Redefine supported_image_types for the HTML builder
from sphinx.builders.html import StandaloneHTMLBuilder
StandaloneHTMLBuilder.supported_image_types = [
  'image/gif', 'image/jpeg', 'image/png', 'image/svg+xml'
]


# -- Options for LaTeX output ---------------------------------------------
latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'utf8extra': '',
    'inputenc': '',
    'babel': r'''\usepackage[english]{babel}''',
    'preamble': r'''\usepackage{ctex}
\definecolor{EOS_green}{rgb}{0.27,0.49,0.36}
    ''',
}
# latex_logo='_static/latex_logo.pdf'
# 设置公式和图片编号依赖于章节
numfig = True
math_numfig = True
math_eqref_format = '({number})'
# 只对make latex有效
# numfig_format = 'Figure. %s'
numfig_secnum_depth = 1
imgmath_latex = 'dvilualatex'
imgmath_image_format = 'svg'
imgmath_dvipng_args = ['-gamma', '1.5', '-bg', 'Transparent']
# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_name='SaltWaterEOS'
latex_documents = [
    (master_doc, latex_name+'.tex', '\\textcolor{EOS_green}{S}alt \\textcolor{EOS_green}{W}ater \\textcolor{EOS_green}{E}quation \\textcolor{EOS_green}{o}f \\textcolor{EOS_green}{S}tate Manual',
     'Zhikui Guo, Lars Rüpke', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'SaltWaterEOS', 'SaltWaterEOS Manual',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Salt Water EOS', 'Salt Water EOS Documentation',
     author, 'Salt Water EOS', 'One line description of project.',
     'Miscellaneous'),
]

def setup(app):
    app.add_css_file("style.css")
    # app.add_javascript("js/custom.js")
    app.add_js_file(
        "https://cdn.jsdelivr.net/npm/clipboard@1/dist/clipboard.min.js")
    
# new defined cite style
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.plugin import register_plugin
from collections import Counter
import re
import unicodedata

from pybtex.style.labels import BaseLabelStyle

_nonalnum_pattern = re.compile('[^A-Za-z0-9 \-]+', re.UNICODE)

def _strip_accents(s):
    return "".join(
        (c for c in unicodedata.normalize('NFD', s)
            if not unicodedata.combining(c)))

def _strip_nonalnum(parts):
    """Strip all non-alphanumerical characters from a list of strings.

    >>> print(_strip_nonalnum([u"ÅA. B. Testing 12+}[.@~_", u" 3%"]))
    AABTesting123
    """
    s = "".join(parts)
    return _nonalnum_pattern.sub("", _strip_accents(s))

class APALabelStyle(BaseLabelStyle):
    def format_labels(self, sorted_entries):
        labels = [self.format_label(entry) for entry in sorted_entries]
        count = Counter(labels)
        counted = Counter()
        for label in labels:
            if count[label] == 1:
                yield label
            else:
                yield label + chr(ord('a') + counted[label])
                counted.update([label])

    def format_label(self, entry):
        label = "Anonymous"
        if 'author' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['author'])
        elif 'editor' in entry.persons:
            label = self.format_author_or_editor_names(entry.persons['editor'])
        elif 'organization' in entry.fields:
            label = entry.fields['organization']
            if label.startswith("The "):
                label = label[4:]

        if 'year' in entry.fields:
            return "{}, {}".format(label, entry.fields['year'])
        else:
            return "{}, n.d.".format(label)

    def format_author_or_editor_names(self, persons):
        if len(persons) is 1:
            return _strip_nonalnum(persons[0].last_names)
        elif len(persons) is 2:
            return "{} & {}".format(
                _strip_nonalnum(persons[0].last_names),
                _strip_nonalnum(persons[1].last_names))
        else:
            return "{} et al.".format(
                _strip_nonalnum(persons[0].last_names))

class APAStyle(UnsrtStyle):

    default_label_style = APALabelStyle

register_plugin('pybtex.style.formatting', 'apa', APAStyle)