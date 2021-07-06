# Figure 1 and supplementary figure 1

Figure 1B and 1C as well as the supplementary figure 1B and 1C were created based on the script in this folder.

## Figure 1B and supplementary figure 1
See [figure_1C.html](./figure_1C.html). The dependencies are in [renv.lock](./renv.lock).
Supplementary figure 1D is taken from [https://danio-code.zfin.org/visualization/#sizes_chart](https://danio-code.zfin.org/visualization/#sizes_chart).
## Figure 1C
Bigwig files for the tracks were first shifted with [bin/multitrack.R](bin/multitrack.R) and converted into a data_hub file for the [WashU Epigenome Browser](http://epigenomegateway.wustl.edu/) using [bin/generate_datahub.py]. You can see an interactive session with all tracks by uploading the session file [WashU-track-session.json](./WashU-track-session.json) to the WashU Epigenome browser.


