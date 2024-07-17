# PRIMA-Panel
Welcome to the pages of the <b>Pr</b>e-Analytical <b>Investigator</b> for NMR-based <b>M</b>et<b>a</b>bolomics

The PRIMA-Panel is tool to investigate the effect processing delays on metabolic parameters in samples of peripheral blood (plasma / serum).
Linear mixed models were used to estimate the change for each metabolic parameter. The data is split into pre- and postcentrinfugation delays
The data is presented in so called stability timepoints. Such a timepoint is defined as the time it takes for a parameter to change by a specific percentage. For example: A value of 0.2 for lactic acid in serum would mean it would take 0.2 hours when the % threshold is set to 20 % change.
Further, the effect of these processing delays can be explored entering real pre-analytical data and observing the direct effect of the dealys on the parameterts. Here, interactive HTML reports can be created.

# Modules

## Data
The Data tab shows different ways of highlighting the different stability time-points in minutes. The time-points are sorted according to their SPREC classification. In addition, there is an alternative way of presenting the date in form of lollipop plots.
The data is split according the two different delays (pre- and post-centrifugation

## Performance Reports
The 'Single' tab allows the user to set a pre-centrifugation time and a post-centrifugation using a Slider. A table is then generated that higlights a minor and a major change for each metabolite in that given timeframe. The colors, as well as the % threshholds, can be adjusted. The table can be downloaded as a .HTML report or the table directly as .csv (or other formats).

## Tools
An additional tab for Tools is available. Currently it consists of a tool for calculating the pre- and post-centrifugation times from the differences of date+time stamps. 

---

# Online Version
This Shiny app is also permantently hosted on Shinyapps.io [here](https://funkam.shinyapps.io/QC-Tool/).

