# Config
CSL=american-physics-society.csl
BIBLIOGRAPHY=bibliography.bibtex

# Input files
MARKDOWN=$(shell ls *.md)
PREAMBLE=$(MARKDOWN:%.md=preamble/%.tex)
TEMPLATE=$(MARKDOWN:%.md=template/%.tex)

GRAPHS=$(ls graphs/**/*)

# Output files
PDF=$(MARKDOWN:%.md=%.pdf)

all: $(PDF) $(CSL) $(TEMPLATE) $(PREAMBLE)

%.pdf: %.md preamble/%.tex template/%.tex $(BIBLIOGRAPHY) $(CSL)
	pandoc $< \
		-o $@ \
		--citeproc \
		--bibliography=$(BIBLIOGRAPHY) \
		--csl=$(CSL) \
		--pdf-engine=xelatex \
		--template template/$(*F).tex \
		-H preamble/$(*F).tex

%.csl:
	wget https://raw.githubusercontent.com/citation-style-language/styles/master/$@

template/%.tex:
	pandoc --print-default-template=latex > $@

preamble/%.tex:
	touch $@
