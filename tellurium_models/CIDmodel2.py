# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:39:36 2018

@author: Joshua Ip - Work
"""

import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import roadrunner
import antimony
import time

r = te.loada("""    
    J0: $AncBinDNA -> AncBinRNA ; a_rna
    J1: $DimBinDNA -> DimBinRNA ; a_rna
    # transcription
    # units of RNA/sec
    
    J2: AncBinRNA -> ; d_rna * AncBinRNA
    J3: DimBinRNA -> ; d_rna * DimBinRNA
    # mRNA decay
    # units of RNA/sec

    J4: AncBinRNA -> AncBinRNA + AncBinder ; a_nb * AncBinRNA
    J5: DimBinRNA -> DimBinRNA + DimBinder ; a_nb * DimBinRNA
    # translation
    # units of RNA/sec
        
    J6: AncBinder -> ; d_nb * AncBinder
    J7: DimBinder -> ; d_nb * DimBinder
    J15: Dimer -> DecDimer; d_nb * Dimer
    # protein decay
    

    
    J8: Mol + AncBinder -> Complex ; a_NBC * (Mol / CytoplasmVol)  * (AncBinder / CytoplasmVol) - d_NBC * (Complex / CytoplasmVol)
    # the anchor binder binds to molecule of interest to form a complex.
    # nanobody complexes may dissociate over time

    J9: Complex + DimBinder -> Dimer ; a_d * (Complex / CytoplasmVol) * (DimBinder / CytoplasmVol) - d_d * (Dimer / CytoplasmVol)
    # dimerization binder binds to complex to form dimers       
    # dimers may dissociate, but much less often than complexes
    
    J10: Dimer + GeneOff -> GeneOn; (GeneOff / NucleusVol) * (Dimer / NucleusVol) * a_g - (GeneOn / NucleusVol) * d_g
    # dimer acts as transcription factor for a gene
    
    J11: GeneOn -> GeneOn + RepRNA; a_rna * GeneOn
    J12: RepRNA -> ; d_rna * RepRNA 
    J13: RepRNA -> RepRNA + Rep; a_nb * RepRNA
    J14: Rep -> ; d_nb * Rep
    # the activated gene transcribes a reporter
    
    # *****************************************************************************************************************************
    # Parameters
    
    TotalCellVol = (30.3 * 10^-6);
    NucleusVol = (4.3 * 10^-6);
    CytoplasmVol = TotalCellVol - NucleusVol;
    # all volumes given in units of L, 
    # volumes from http://bionumbers.hms.harvard.edu/bionumber.aspx?id=106557&ver=1&trm=yeast%20cytoplasm%20volume&org=
    
    scalingFactor = 60 * 60;
    # since all our rates/rate constants are in seconds, we can scale time by multiplying each time-dependent parameter by a scaling factor
    # this particular value scales the parameters for time units of hours
    
    a_rna = (0.002) * scalingFactor;
    # median transcription rate = 0.12 mRNA molecules/min = 0.002 mRNA molecules/sec
    # median transcription rate from http://bionumbers.hms.harvard.edu/bionumber.aspx?id=106766&ver=3&trm=transcription%20rate%20yeast&org=
    # KEY ASSUMPTION: the rate of transcription of our nanobody gene is constant. 
    # in reality, it may not be safe to assume that our molecule is transcribed by the median transcription rate
    
    d_rna = (5.6 * 10^-4) * scalingFactor;        
    # 5.6 * 10 ^ -4 = mRNA decay rate constant in units of sec^-1
    # mRNA decay constant found from http://bionumbers.hms.harvard.edu/bionumber.aspx?id=105510&ver=5&trm=mrna%20s.%20cerevisiae&org=
    
    a_nb = (0.0185) * scalingFactor;
    # yeast has no rough ER, so translation occurs in the cytoplasm
    # median time for translation initiation = 4.0 * 10^2 s
    # median elongation rate = 9.5 aa/s
    # nanobody average amino acids = 130 aa
    # time for elongation = (130 aa)/(9.5 aa/s) = 14 sec
    # total time for 1 mRNA transcript = 14 sec + 40 sec = 54 sec
    # rate at which mRNA is transcribed = 1 protein/(54 sec * 1 mRNA) / ~ 0.0185 protein/(sec mRNA)
    # it is notable that translation initiation rate can vary between mRNA by orders of magnitude
    # all data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694300/

    d_nb = (2.6 * 10^-4) * 5 * scalingFactor;
    # d_nb represents degradation of the nanobodies. Hannah found http://www.pnas.org/content/103/35/13004 (doi: https://doi.org/10.1073/pnas.0605420103) https://www.nature.com/articles/nature10098,
    # which shows that the median half-life of a protein in a budding yeast cell is 43 minutes
    # half-life is cop6e-4 = median rate constant of degradation of proteins in a yeast cell
    
    a_NBC = 0.05 * 10^6 * scalingFactor; 
    d_NBC = 20 * 10^-6 * scalingFactor; 
    # this is one of the binding affinities that we will do a parameter sweep to learn more about
    # current stand-in value is the rate constant of association/dissociation for antibody binding to cytochrome C

    # Binding affinity of the dimerization nanobody
    # will be changed using a parameter sweep, units of M/s
    a_d = 0.01 * 10^5 * scalingFactor; 
    d_d = 100 * 10^-5 * scalingFactor;
    
    # binding affinity of the completed transcription factor. This depends mostly on the DNA-binding domain chosen.
    # 7.0e9 = , association constant of the lac repressor to the lac operon, units of M/s
    a_g = 7.0 * 10 ^ 9 * scalingFactor;
    
    # 12e-12 = k_d of Egr1 DNA binding domain, units of M/s
    # found at http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=5&id=104597
    d_g = 12 * 10 ^ -12 * scalingFactor;
    
    # *****************************************************************************************************************************************
    # Initial values
    # These are all in copies
    AncBinDNA = 1; 
    DimBinDNA = 1;
    Mol = 0;
    GeneOff = 1;
    
""")
r.reset()
#r.draw(width=800,height=300,overlap = "false", splines = "true")

<<<<<<< HEAD
prepertubation = r.simulate(0, 1, 10000)
r.Mol = 50
pertubation = r.simulate(1, 2, 10000)
#r.Mol = 0
postpertubation = r.simulate(2, 10, 10000)
=======
prepertubation = r.simulate(0, 12, 1000)

r.x[0] = 50
pertubation = r.simulate(12, 36, 1000)
r.Mol = 0
postpertubation = r.simulate(36, 48, 1000)
>>>>>>> d807695917cb2def7dc63a60ef963b9667a40ad6
result = numpy.vstack((prepertubation, pertubation, postpertubation))

plt.figure(1)
plt.plot(result[:,0], result[:,1:8])
plt.xlabel("Time (hours)")
plt.ylabel("Copies")
plt.legend(["Anchor Binder RNA", "Dimerization Binder RNA", "Anchor Binder", "Dimerization Binder", "Molecule", "Complex", "Dimer"])

plt.figure(2)
plt.plot(result[:,0],result[:,8:10])
plt.xlabel("Time (hours)")
plt.ylabel("Copies")
plt.legend(["GeneOff", "GeneOn"])

plt.figure(3)
plt.plot(result[:,0],result[:,11])
plt.xlabel("Time (hours)")
plt.ylabel("Reporter, copies")
