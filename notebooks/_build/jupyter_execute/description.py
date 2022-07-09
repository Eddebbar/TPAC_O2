#!/usr/bin/env python
# coding: utf-8

# # Description 
# 
# This book describes the workflow for project titled: "the tropical Pacific Oxygen at the Mesoscale"   
# 
# This book is based off the [Jupyter Book](https://jupyterbook.org/en/stable/intro.html) template and uses  [papermill](https://papermill.readthedocs.io/en/latest/index.html) which enables parameterizing and executing notebooks.
# 
# To run all notebooks, excecute [_run.ipynb](https://github.com/Eddebbar/TPAC_O2/blob/main/notebooks/_run.ipynb).

# # Model outputs are found in the following local directories

# - CESM JRA55 0.1º daily runs: '/glade/campaign/cesm/development/bgcwg/projects/hi-res_JRA/cases/g.e22.G1850ECO_JRA_HR.TL319_t13.004/output/ocn/proc/tseries/day_1/' 
# - CESM JRA55 0.1º monthly runs: '/glade/campaign/cesm/development/bgcwg/projects/hi-res_JRA/cases/g.e22.G1850ECO_JRA_HR.TL319_t13.004/output/ocn/proc/tseries/monthly_1/'     
# - CESM CORE 0.1º 5yr Runs:  '/glade/campaign/cesm/development/bgcwg/projects/hi-res_CESM1_CORE/g.e11.G.T62_t12.eco.006/ocn/hist/'
# - B-TPOSE: 
#     - Mitgcm + BLING 2010-2018: http://www.ecco.ucsd.edu/DATA/TROPAC/bgc/
#     - TPOSE + BLING : http://www.ecco.ucsd.edu/DATA/TROPAC/tpose6_bgc/  
#     - TPOSE Physics only: http://www.ecco.ucsd.edu/DATA/TROPAC/tpose6/
#     - wget into '/glade/work/yeddebba/TPOSE/
# 
