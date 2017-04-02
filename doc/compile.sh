# /bin/bash

# Options:
# --no-generator          Do not include a generator credit.
# --date, -d              Include the date at the end of the document (UTC).
# --toc-top-backlinks     Enable backlinks from section headers to the top of
#                         the table of contents.
# --footnote-backlinks    Enable backlinks from footnotes and citations to their
#                         references.  This is the default.
# --section-numbering     Enable Docutils section numbering (default: enabled).
# --tab-width=<width>     Set number of spaces for tab expansion (default 8).
# --documentclass=DOCUMENTCLASS
#                         Specify documentclass.  Default is "article".
# --documentoptions=DOCUMENTOPTIONS
#                         Specify document options.  Multiple options can be
#                         given, separated by commas.  Default is
#                         "10pt,a4paper".
# --use-latex-footnotes   Use LaTeX footnotes. LaTeX supports only numbered
#                         footnotes (does it?). Default: no, uses figures.
# --footnote-references=<format>
#                         Format for footnote references: one of "superscript"
#                         or "brackets".  Default is "superscript".
# --stylesheet=<file>     Specify a stylesheet file. The file will be "input" by
#                         latex in the document header.  Default is no
#                         stylesheet ("").  Overrides --stylesheet-path.
# --use-verbatim-when-possible
#                         When possibile, use verbatim for literal-blocks.
#                         Default is to always use the mbox environment.
# --table-style=<format>  Table style. "standard" with horizontal and vertical
#                         lines, "booktabs" (LaTeX booktabs style) only
#                         horizontal lines above and below the table and below
#                         the header or "nolines".  Default: "standard"
# --font-encoding=FONT_ENCODING
#                         LaTeX font encoding. Possible values are "T1", "OT1",
#                         "" or some other fontenc option. The font encoding
#                         influences available symbols, e.g. "<<" as one
#                         character. Default is "" which leads to package "ae"
#                         (a T1 emulation using CM fonts).

DOCSBASELIST="user-guide"
STYLESHEET="--stylesheet=RST3"
DOCOPTIONS="12pt,a4paper"

for DOCSBASENAME in ${DOCSBASELIST} ; do
    if [ -f "${DOCSBASENAME}.txt" ] ; then
	rst2latex --language=en --report=4 \
	    --section-numbering ${STYLESHEET} --use-verbatim-when-possible \
	    --graphicx-option=pdftex --font-encoding=T1 \
	    --documentoptions=${DOCOPTIONS} \
	    ${DOCSBASENAME}.txt ${DOCSBASENAME}.tex

	./tolistings.awk ${DOCSBASENAME}.tex > pp.tex
	mv pp.tex ${DOCSBASENAME}.tex

	pdflatex -src-specials -file-line-error-style ${DOCSBASENAME}.tex
	while grep 'Rerun to get' ${DOCSBASENAME}.log > /dev/null
	do
	    pdflatex -src-specials -file-line-error-style ${DOCSBASENAME}.tex
	done
	pdflatex -src-specials -file-line-error-style ${DOCSBASENAME}.tex
    fi
done
