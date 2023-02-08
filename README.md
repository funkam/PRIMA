# NMR Metabolomics Quality Control Panel
A tool for investigating the quality and stability of metabolites in peripheral blood.

Through an experimental setup which consists of leaving samples with different pre- and post-centrifugation times for up to 8 hours, linear mixed models were created for each metabolite and its stability with pre- and post-centrifugation times were established.

These models were used to create an interface using Shiny that allows to highlight the different percentual changes over time for a given combination of pre- and post-centrifugation times. The user can select the color and also the color thresholds (default 15 and 30%) to highlight specific changes. The table can be downloaded in form of a interactive html report (R Markdown), but also the table directly in a variety of formats.

This Shiny app is also permantently hosted on Shinyapps.io [here](https://funkam.shinyapps.io/QC-Tool/).

