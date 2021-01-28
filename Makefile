all: committee-doc article

committee-doc: docs/committee/committee-doc.pdf
docs/committee/committee-doc.pdf: docs/committee/committee-doc.md docs/committee/fig/*
	(cd docs/committee/; pandoc -s --filter pandoc-crossref --citeproc committee-doc.md -o committee-doc.pdf)

article: docs/article/article.pdf
docs/article/article.pdf: docs/article/article.md docs/article/fig/*
	(cd docs/article/; pandoc -s --filter pandoc-crossref --citeproc article.md -o article.pdf)
