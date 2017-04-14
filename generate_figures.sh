#!/bin/bash

tar -xvf data.tar.gz

Rscript src/figure1.R
Rscript src/figure2.R
Rscript src/figure3.R
Rscript src/figure4.R
Rscript src/figure5.R

rm Rplots.pdf