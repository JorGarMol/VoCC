# VoCC: The Velocity of Climate Change and related climatic metrics

Jorge Garcia Molinos et al. 18 July 2019

This package is now in release version (v 1.0.0). Please contact [me](jorgegmolinos@arc.hokudai.ac.jp) for questions or feed-back.


## Installation

To install the package, open R and type:

    install.packages("devtools")

Then, you can install vocc:

    devtools::install_github("JorGarMol/VoCC", dependencies = TRUE, build_vignettes = TRUE)
    # Alternatively call without building the vignette to save installation time
    library(VoCC)
    citation("VoCC")

## License

GNU Affero General Public License v3.0

## Citation

To cite the package itself please use:

García Molinos, J., Schoeman, D. S., Brown, C. J. and Burrows, M. T. (2019). VoCC: The Velocity of Climate Change and
  related climatic metrics. R package version 1.0.0. https://doi.org/10.5281/zenodo.3382092

The following paper explains the package functionality and provides the examples covered by the code in the package vignette

García Molinos, J., Schoeman, D. S., Brown, C. J. and Burrows, M. T. (2019), VoCC: An R package for calculating the velocity of climate change and related climatic metrics. Methods Ecol Evol. doi:10.1111/2041-210X.13295
