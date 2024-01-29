#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:10:12 2022

@author: kevin
"""

from astropy.io.votable import parse,from_table,writeto
import numpy as np
import matplotlib.pyplot as plt

name='Ruprecht_147'
votable = parse(name+"_pdf_iso_diff_old.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)
number=2
# data=data[:10]
data=data[number]
plt.figure()
temperature_star=data['teff_raw']
logg_star=data['logg_raw']
teff_error=data['e_teff_raw']
logg_error=data['e_logg_raw']
# plt.scatter(temperature_star,logg_star)
plt.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)
xx=data['teff_array']
yy=data['logg_array']
zz=data['probability_grid']
plt.contour(xx,yy,zz)
plt.xlabel('teff')
plt.ylabel('logg')    
