# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 22:28:12 2023

@author: M Elghorab
"""

from engine.fr_entities_tools import *

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns



sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(dpi=750)

Oel = R_Observation('C:/Users/User/Downloads/Oelsnitz.csv').fr.plot(label = 'Oelsnitz')

Zit = R_Observation('C:/Users/User/Downloads/Zittau.csv').fr.plot(label ='Zittau')

Doh = R_Observation('C:/Users/User/Downloads/Dohna.csv').fr.plot(label = 'Dohna')

ax.legend()
ax.set_ylabel('Discharge (m\u00b3/s)')