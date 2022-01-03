# covidvariants
[![DOI](https://zenodo.org/badge/377944722.svg)](https://zenodo.org/badge/latestdoi/377944722)

This repository contains data and code that support
Miller, K., K. Elenberg, and A. Dubrawski. "Forecasting emergence of COVID-19 variants of concern." PLOS ONE [under review].

This library contains a slightly modified version of MutAntiGen https://github.com/davidrasm/MutAntiGen
This version reports infections by antigenic type
out.infectionsByPhenotype <-- number of infections and recoveries by day and type

launch.py samples parameter values and runs many MutAntiGen simulations.

ProcessSims.R summarizes the simulation results.

