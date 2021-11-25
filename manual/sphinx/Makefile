# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = swEOS
SOURCEDIR     = source
BUILDDIR      = build
BUILDDIR_ZH      = build_zh

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).

# build zh: first run `make gettext`
html:
	sphinx-intl build
	sphinx-build -b $@ source build/$@/en
	# sphinx-build -b $@ source build/$@/zh -Dlanguage=zh
	cp index.html build/$@
latex:
	sphinx-intl build
	sphinx-build -b $@ source build/$@/en
	sphinx-build -b $@ source build/$@/zh -Dlanguage=zh

gettext:
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	sphinx-intl update --language=zh 

clean:
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)