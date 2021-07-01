all: committee-doc article

committee-doc: committee-doc.pdf
docs/committee/committee-doc.pdf: committee-doc.md fig/*
	pandoc -s --filter pandoc-crossref --citeproc committee-doc.md -o committee-doc.pdf
