library(rmarkdown)
library(quarto)
fs <- list.files(pattern="\\.qmd$")
for (f in fs) {
    print(f)
    quarto::quarto_render(f)
}

