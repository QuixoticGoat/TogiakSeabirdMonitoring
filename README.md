---
title: "README"
author: McCrea Cobb <mccrea_cobb@fws.gov>
date: 6/7/2018
output: html_notebook
---
The United States Fish and Wildlife Service (FWS) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. FWS has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by FWS. The FWS seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by FWS or the United States Government.

The abundance and reproductive success of black-legged kittiwakes, common murres, and pelagic cormorants was monitored annually at Cape Peirce from 1990-2014 and 2016-2017. The goal of this work is to quantify temporal trends in abundance and identify environmental factors associated with changes in abundance. 

------------------------

**Contents**

*./Scripts/1_ImportFormat.R*  
Imports a .csv file containg counts of seabirds and next by plot annually. Reformats data for plotting and analysis. Outputs the formatted data as a .RData file.  

*./Scripts/2_Plots.R*
Creates summary plots of the data, including a lineplots and boxplots of annual seabird and nest abundance. 

*./Scripts/3_Models.R*  
Generalized mixed effects models for quantifying the response of seabird and nest abundance to time (Year) and environmental variables.
