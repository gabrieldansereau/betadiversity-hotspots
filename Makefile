doc/committee-doc.pdf: doc/committee-doc.md
	(cd doc; pandoc -s --filter pandoc-citeproc committee-doc.md -o committee-doc.pdf)
