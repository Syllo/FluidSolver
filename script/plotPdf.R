#!/bin/Rscript

# Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

if (!require("ggplot2", quiet = TRUE)) {
  install.packages("ggplot2", dependencies=TRUE)
  library("ggplot2")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Input data file name not provided\n", call. = FALSE)
}

if (length(args) == 2) {
  out_file <- args[2]
} else {
  out_file <- paste(args[1], ".pdf", sep = "")
}

heatTable <- read.table(args[1], col.names=c("xpos","ypos","fluidDensity","fluidVelX","fluidVelY"))
pInitial <- ggplot(heatTable, aes(x=xpos, y=ypos)) + scale_x_continuous("Physical domain (x dimension)") + scale_y_continuous("Physical domain (y dimension)") + ggtitle("Fluid Density Heatmap")
p <- pInitial + geom_raster(aes(fill=fluidDensity),show.legend=TRUE) + scale_fill_viridis_c("Fluid Density", limits=range(heatTable$fluidDensity),option="C")

pdf(out_file)
p
invisible(dev.off())
#embedFonts(out_file)
