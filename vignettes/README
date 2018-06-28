# For producing a JSS article with R/Latex, please start from the following template: jss-article-rnw.zip (including all necessary style files). The source files for the JSS contribution are the article.Rnw manuscript and the ref.bib bibliography in BibTeX format. To produce the PDF manuscript the article must first be weaved:

  Sweave("article.Rnw")
  # yielding article.tex along with the graphics file (here article-visualization.pdf). This LaTeX file can be compiled in the "usual" way to the PDF using pdfLaTeX, i.e., using the shell, some LaTeX editor (see also above), or simply with R:
  
  library("tools")
  texi2pdf("article.tex")
# An alternative route for RStudio users is to use the "Compile PDF" button directly for the article.Rnw. However, make sure not to use concordance mode: Tools > Global Options > Sweave > Uncheck: Always enable Rnw concordance.

# The replication code can be obtained by tangling:

Stangle("article.Rnw")
# yielding article.R. However, editing the comments of this .R file typically makes it more accessible to the readers/reviewers.
# 
# If knitr is used to prepare the article.Rnw file, the render_sweave() hook should be used to ensure JSS style code formatting. The weaving and tangling can then by done by knit() and purl(), respectively.