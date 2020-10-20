# Make vignette previews to view on GitHub.

library(rmarkdown)

render("Setup.Rmd", 
       output_format = github_document(toc = FALSE, html_preview = FALSE))

render("Pedigree_verification.Rmd", 
       output_format = github_document(toc = FALSE, html_preview = FALSE))
