# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:39:36 2018

@author: Joshua Ip - Work
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te
from tellurium import ParameterScan as ps
import roadrunner
import antimony
import time


antimonyString = ("""    
    J0: $AncDNA -> AncRNANuc ; a_rna * AncDNA
    J1: $DimDNA -> DimRNANuc ; a_rna * DimDNA
    # transcription
    # units of (mRNA copies)/(sec)

    J3: AncRNANuc -> AncRNACyt ; diffusion_rna * AncRNANuc - diffusion_rna * AncRNACyt   
    J2: DimRNANuc -> DimRNACyt ; diffusion_rna * DimRNANuc - diffusion_rna * DimRNACyt
    # mRNA transport out of the nucleus into the cytoplasm
    # units of (mRNA copies)/(sec)
    
    J4: AncRNACyt -> ; d_rna * AncRNACyt
    J5: DimRNACyt -> ; d_rna * DimRNACyt
    J6: AncRNANuc -> ; d_rna * AncRNANuc
    J7: DimRNANuc -> ; d_rna * DimRNANuc
    # mRNA decay
    # units of 1/(sec) * (mRNA copies) = (mRNA copies)/(sec)

    J8: -> AncBinder ; a_nb * AncRNACyt
    J9: -> DimBinder ; a_nb * DimRNACyt
    # translation
    # units of (protein copies)/(sec * mRNA copies) * (mRNA copies) = (protein copies / sec)

        
    J10: AncBinder -> ; d_nb * AncBinder
    J11: DimBinder -> ; d_nb * DimBinder
    J12: DimerCyt -> ; d_nb * DimerCyt
    J13: DimerNuc -> ; d_nb * DimerNuc
    # protein decay
    # units of (1 / sec) * (protein copies) = (protein copies / sec)
    

    
    J14: Mol + AncBinder -> Complex ; AncBinder * ( (Mol ) / (K_d_NBC + Mol / AvoNum / CytoplasmVol) ) 
    # the anchor binder binds to molecule of interest to form a complex.
    # nanobody complexes may dissociate over time

    J15: Complex + DimBinder -> DimerCyt ; DimBinder * ( (Complex) / (K_d_Dimerization  + Complex / AvoNum / CytoplasmVol) ) 
    # dimerization binder binds to complex to form dimers       
    # dimers may dissociate, but much less often than complexes
    # This is governed by a variant of the hill equation
    
    J16: DimerCyt -> DimerNuc; diffusion_nb * DimerCyt
    J17: DimerNuc -> DimerCyt; diffusion_nb * DimerNuc
    # dimer must be transported into the cell to act as a transcription factor
    #
    
    J18: DimerNuc + GeneOff -> GeneOn; GeneOff * ( (DimerNuc) / (K_d_Gene   + DimerNuc / AvoNum / NucleusVol) )
    # dimer acts as transcription factor for a gene
    # units: (copies) / (copies)
    
    J19: -> RepRNANuc; a_rna * GeneOn
    J20: RepRNANuc -> RepRNACyt; diffusion_rna * RepRNANuc - diffusion_rna * RepRNACyt
    J21: RepRNANuc -> ; d_rna * RepRNANuc
    J22: RepRNACyt -> ; d_rna * RepRNACyt
    J23: -> Rep; a_nb * RepRNACyt
    J24: Rep -> ; d_nb * Rep
    # the activated gene transcribes a reporter
    
    # *****************************************************************************************************************************
    # Parameters
    
    AvoNum = 6.02 * 10^23;
    
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
    # median time for translation initiation = 4.0 * 10^2 s * mRNA / protein
    # median elongation rate = 9.5 aa/s
    # nanobody average amino acids = 130 aa
    # time for elongation = (130 aa / protein)/(9.5 aa/s) = 14 sec / protein
    # total time for 1 mRNA transcript = 14 sec / protein + 40 sec = 54 sec
    # rate at which mRNA is transcribed = 1 protein/(54 sec * 1 mRNA) / ~ 0.0185 protein/(sec mRNA)
    # it is notable that translation initiation rate can vary between mRNA by orders of magnitude
    # all data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694300/

    d_nb = (2.6 * 10^-4) * 5 * scalingFactor;
    # which shows that the median half-life of a protein in a budding yeast cell is 43 minutes
    # median rate constant of degradation of proteins in a yeast cell = 2.6e-4 1/sec
    # data from http://www.pnas.org/content/103/35/13004 (doi: https://doi.org/10.1073/pnas.0605420103) https://www.nature.com/articles/nature10098,
    
    bar = 20;
    K_d_NBC = bar * 10^-6 * scalingFactor;
    
    #kOn_NBC = 1.0 * 10^3 * scalingFactor; 
    #kOff_NBC = 20 * 10^-3 * scalingFactor; 
    # this is one of the binding affinities that we will do a parameter sweep to learn more about
    # current stand-in value is the rate constant of association/dissociation for antibody binding to cytochrome C
    
    foo = 100
    K_d_Dimerization = foo * 10 ^ -9
    
    # Binding affinity of the dimerization nanobody
    # will be changed using a parameter sweep, units of 1/(sec * M)
    # kOn_Dimer = 1.0 * 10^2 * scalingFactor; 
    # kOff_Dimer = 100 * 10^-7 * scalingFactor;
        
    K_d_Gene = 12 * 10 ^ -12 * scalingFactor;
    # 12e-12 = k_d of Egr1 DNA binding domain, units of 1/(sec * M)
    # data from http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=5&id=104597

    
    diffusion_rna = 1;
    diffusion_nb = 1;
    #
    # data from
    
    # *****************************************************************************************************************************************
    # Initial values
    # These are all in copies
    AncDNA = 1; 
    DimDNA = 1;
    Mol = 0;
    GeneOff = 1;
    Setting = 50;
    
    
    at time>=4: Mol=Setting;
    

    
""");

r = te.loadAntimonyModel(antimonyString)
#r.draw(width = '1800')
r.simulate(0, 12, 1000)
p = te.ParameterScan(r,
                        # Settings
                        startTime = 0,
                        endTime = 10,
                        numberOfPoints = 10,
                        value = 'bar',
                        startValue = 100,
                        endValue = 2000,
#
#                        value = 'Setting',
#                        startValue = 0,
#                        endValue = 50, 
                        independent = ['K_d_NBC','K_d_Dimerization'],
                        dependent = 'Rep',
                        xlabel = 'Time',
                        ylabel = 'K_d_NBC (uM)',
                        zlabel = 'Reporter',
                        title = "Model",
                        selection = 'Rep')
#                        value = 'Setting',

p.plotSurface()

