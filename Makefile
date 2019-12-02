doc/committee-doc.pdf: doc/committee-doc.md doc/fig/*
	(cd doc; pandoc -s --filter pandoc-crossref --filter pandoc-citeproc committee-doc.md -o committee-doc.pdf)
