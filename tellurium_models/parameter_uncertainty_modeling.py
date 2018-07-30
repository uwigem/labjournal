#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 19:10:41 2018

@author: zackmcnulty
"""
# -*- coding: utf-8 -*-
# Dallas Warren
# Washington iGEM - Simulations
# Model of chemically induced dimerization.
# Nanobody kinetic rates based on antibody kinetic rates.
# Nanobody concentrations based on a protein's average concentration in cell.
# Molecule concentration as an arbitrary number.
# Gene kinetic rates and concentrations fudged.
"""
Tellurium oscillation
"""
import numpy as np
import tellurium as te
import roadrunner
import antimony
import matplotlib.pyplot as plt
import time

r = te.loada ('''
model feedback()
  // Reactions:
  J0: Nan1 + Mol -> Nan1Mol; (K1*Nan1*Mol);
  J1: Nan1Mol -> Nan1 + Mol; (K_1*Nan1Mol); 
  J2: Nan1Mol + Nan2 -> Nan1MolNan2; (K2*Nan1Mol*Nan2)
  J3: Nan1MolNan2 + GeneOff -> GeneOn; (K3*Nan1MolNan2*GeneOff);
  //J4: GeneOn -> Nan1MolNan2 + GeneOff; (K_3*GeneOn);

  // Species initializations:
  Nan1 = 0.0001692; 
  Mol = 0.0001692/2; 
  Nan2 = 0.0001692; 
  Nan1Mol = 0;
  Nan1MolNan2 = 0; 
  GeneOff = 5*10^-5; 
  GeneOn = 0;

  // Variable initialization:
  const K1, K_1, K2, K3, K_3;
  K1 = 6.1*10^5; 
  K_1 = 8*10^-5; 
  K2 = 3.3*10^5; 
  K3 = 1*10^5; 
  K_3 = 0;
end''')

def plot_param_uncertainty(model, startVal, name, num_sims):
    stdDev = 0.3*startVal
    #vals = np.linspace((1-stdDev)*startVal, (1+stdDev)*startVal, 100)
    vals = np.random.normal(loc = startVal, scale=stdDev, size = (num_sims, ))
    for val in vals:
        exec("r.%s = %d" % (name, val))
        result = r.simulate(0, .5, 1000)
        r.reset();
        plt.plot(result[:,0],result[:,7])
        plt.title("Response to uncertainty in " + name)
    plt.legend(["GeneOn"])
    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration")
    plt.ylim([0, 4*10**-5])
    plt.xlim([0, 0.5])




startVals = r.getGlobalParameterValues();
names = r.getGlobalParameterIds();
n = len(names) + 1;
dim = math.ceil(math.sqrt(n))
for i in range(1,n):
    plt.subplot(dim,dim,i)
    plot_param_uncertainty(r, startVals[i-1], names[i-1], 100)



plt.tight_layout()
plt.gcf().set_size_inches(10,10)
plt.show()