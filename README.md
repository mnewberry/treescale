# treescale

This repository reproduces the analysis of a study of fractal branching of
trees in art.

![three trees](./samples/grove5.svg)

## Data and Reproducibility

To reproduce the analysis from the command-line, with `R` and `ggplot2`
installed, run

```
Rscript analysis.R
```

in the repository root directory. This will generate the file `Rplots.pdf`
which contains each data figure.

To generate a random grove of three trees with different values of alpha,
simply load `grove.svg` in a web browser such as Firefox.

For transparency of the data collection process, the original hand-annotated
SVG files for each artwork by each author (authors `j`, `m` and anonymous
participant `a`) are given in the `samples` directory. `count-svg.js` contains
a javascript snippet which can be used to automatically measure the length of
each annotation line. To use `count-svg.js`, load one of the SVG files such as
`gtm.svg` into a browser and paste the contents of `count-svg.js` into the
javascript console (sometimes called Web Developer Console). This generates the
contents of the corresponding `samples/gtm.nlsv` file as a javascript variable.

## References

Lin, Q. & Newerry, M.G., (2023) Seeing through noise in power laws. J. 
Roy. Soc. Interface, 20(205):20230310.
https://doi.org/10.1098/rsif.2023.0310

Newberry, M. G., & Savage, V. M. (2019). Self-similar processes follow a power
law in discrete logarithmic space. Physical review letters, 122(15), 158303.
https://doi.org/10.1103/PhysRevLett.122.158303

Brummer A.B., Lymperopoulos P., Shen J., Tekin E., Bentley L.P., Buzzard V.,
Gray A., Oliveras I., Enquist B.J. & Savage V.M. (2021) Branching principles of
animal and plant networks identified by combining extensive data, machine
learning and modelling. J. R. Soc. Interface. 18:20200624.
http://doi.org/10.1098/rsif.2020.0624

## Licensing and usage authorization

All code and data collected in this study are copyright (C) Jingyi Gao and
Mitchell Newberry 2023 and hereby released under the under GPL license version
3.0, see LICENSE.

The files `samples/balsa.nlsv`, `samples/pinon.nlsv`, and
`samples/ponderosa.nlsv` are reproduced from Brummer et al. (2021)
supplementary data files, published by the Royal Society under the terms of the
Creative Commons Attribution License.

The images `samples/*.jpg` were obtained from Wikimedia Commons and are
faithful reproductions of works in the public domain and hence used here
without further attribution.
