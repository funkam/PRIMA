# PRIMA-Panel
Welcome to the pages of the <b>Pr</b>e-Analytical <b>I</b>nvestigator for NMR-based <b>M</b>et<b>a</b>bolomics (PRIMA-Panel)
<img align="right" width="200" height="200" src="https://github.com/funkam/QC-Tool/blob/main/www/logo.png">

The PRIMA-Panel is a tool to investigate the effect of processing delays on metabolic parameters in samples of peripheral blood (plasma / serum). <br>
The panel functions as a data expoloration tool. It allows to investigate the these effects interactively. Additionally, the PRIMA-Panel allows the creation of so called performance reports for sample cohorts. Here a data table with pre-analytical information can be uploaded


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

