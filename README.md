# ACRR
ACRR (aggregated causal relation reliability) program

We implement ACRR based on pc from pcalg package in R, by modifying some functions in the pc. For the pcalg package, please reference to https://CRAN.R-project.org/package=pcalg

### Mods

In mods, there are three R scripts, "pc_mod.R," "skeleton_mod.R," and "skeleton_mod.R." Each of them is a function modified from the original function from pcalg, to implement ACRR.

### Experiments

We also include some examples of experiment with program generated data, codes to run these experiments are in "Experiments."

### How to Run

pcalg package (https://CRAN.R-project.org/package=pcalg) and graph package (https://www.bioconductor.org/packages/release/bioc/html/graph.html) are needed for running ACRR.

Run three R scripts from "Mods" first, following steps are similar to codes in "Experiments."
