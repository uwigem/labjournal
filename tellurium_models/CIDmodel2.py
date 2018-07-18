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
    scalingFactor = 60 * 60;

    AR1: $DNA1 -> RNA1 ; a_rna1 * DNA1
    AR2: $DNA2 -> RNA2 ; a_rna2 * DNA2
    # DNA to mRNA, DNA doesn't change in concentration


    DR1: RNA1 -> ; d_rna1 * RNA1
    DR2: RNA2 -> ; d_rna2 * RNA2
    # mRNA degrades significantly quicker than mRNA

    AN1: -> NB1 ; a_nb1 * RNA1 
    AN2: -> NB2 ; a_nb2 * RNA2 
    # RNA to Protein, the protein is our nanobody
    # 
        
    DN1: NB1 -> ; d_nb1 * NB1
    DN2: NB2 -> ; d_nb2 * NB2
    # nanobodies can degrade or be denatured
    
    AC1: -> COMPOUND ; a_c 
    DC1: COMPOUND -> ; d_c * COMPOUND
    #Our compounds can be created or destroyed, but currently a_c
    
    ANBC1: COMPOUND + NB1 -> NBC1 ; a_nbc1 * COMPOUND  * NB1
    # ANBC2: C2 + NB2 -> NBC2; a_nbc2 * C2 * NB2
    # nanobodies bind to compound of interest
    
    DNBC1: NBC1 -> NB1 + COMPOUND ; d_nbc1 * NBC1
    # nanobody complexes may dissociate over time
    
    AD: NBC1 + NB2 -> DIMER ; a_dimer * NBC1 * NB2
    # nanobodies and nanobody-complex combine to form dimers 
    
    # DD: DIMER -> NBC1 + NB2 ; d_dimer * DIMER
    # nanobodies and nanobody-complexes may degrade
    
    GC: DIMER + GENEOFF -> GENEON ; 
    # Once formed, our dimer acts as a transcription factor for our gene of interest
    
    GD: GENEON -> GENEOFF + DIMER ;
    # The transcription factor can unbind from the gene
    
    #These represent the repressilator that is activated once the gene is turned on.
    
    AMR1: $DNA_REP1 -> RNA_REP1 ; DNA_REP1 * a_rna_rep1
    AMR2: $DNA_REP2 -> RNA_REP2 ; DNA_REP2 * a_rna_rep2
    AMR3: $DNA_REP3 -> RNA_REP3 ; DNA_REP3 * a_rna_rep3

    DMR1: RNA_REP1 -> ; DNA_REP1 * d_rna_rep1
    DMR2: RNA_REP2 -> ; DNA_REP2 * d_rna_rep2
    DMR3: RNA_REP3 -> ; DNA_REP3 * d_rna_rep3

    AR1: -> REP1 ; RNA_REP1 * a_rep1
    AR2: -> REP2 ; RNA_REP2 * a_rep2 
    AR3: -> REP3 ; RNA_REP3 * a_rep3

    DR1: REP1 -> ; REP1 * d_rep1
    DR2: REP2 -> ; REP2 * d_rep2
    DR3: REP3 -> ; REP3 * d_rep3

    

    # Parameters
    
    # a_rna represents the time for mRNA to be created for the nanobody
    # (4000 genes in E.coli K12/ 20 genes/min initiation rate) = 20 minutes for initiation by e.coli 
    # (length of 351 nucleotides per nanobody / 45 nts/sec elongation) = 7.8 seconds for mRNA formation
    # overall it should take 1207.8 seconds for 1 mRNA, so 1/1207.8 = 0.000828
    a_rna1 = 0.000828 * scalingFactor; 
    a_rna2 = 0.000828 * scalingFactor; 

    d_rna1 = 0.00465 * scalingFactor;
    d_rna2 = 0.00465 * scalingFactor; # half-life of 2.5 minutes, from bionumbers
    
    a_nb1 = 0.0346 * scalingFactor;
    a_nb2 = 0.0346 * scalingFactor; # average initiation time of 15 seconds + mean elongation time for 1 codon of 0.119 sec * 117 codons
    d_nb1 = 0.000138 * scalingFactor;
    d_nb2 = 0.000138 * scalingFactor; # half-life of protein at 2 hours, from bionumbers
    
    a_nbc1 = 6 * 10 ^ 5 * scalingFactor; # rate constant of association of antibody binding to cytochrome C
    d_nbc1 = 8 * 10 ^ -5 * scalingFactor; # rate constant of dissociation for anitbody binding to cytochrome C   
    
    # a_dimer was derived from the rate constant of assocation for hemoglobin dimers, 0.4 * 10^5 
    a_dimer = 4.0 * 10^5 * scalingFactor; 
    #d_dimer = ?;

    # we could simulate our compound being created or destroyed, but for now we assume it is constant
    a_c = 0 * scalingFactor; 
    d_c = 0 * scalingFactor; 

    
    #Though I used the same values here as for the nanobody itself, the promoter sequence may be stronger or weaker in a real system. 
    #That is why I have set up a seperate variable to track the strength of the repressilator promoter.
    a_rna_rep1 = 0.000828 * scalingFactor; 
    a_rna_rep2 = 0.000828 * scalingFactor; 
    a_rna_rep3 = 0.000828 * scalingFactor; 

    d_rna_rep1 = 0.00465 * scalingFactor; 
    d_rna_rep2 = 0.00465 * scalingFactor; 
    d_rna_rep3 = 0.00465 * scalingFactor; 
    
    a_rep1 = 0.0346 * scalingFactor;
    a_rep2 = 0.0346 * scalingFactor;
    a_rep3 = 0.0346 * scalingFactor;
    
    d_rep1 = 0.000138 * scalingFactor;
    d_rep2 = 0.000138 * scalingFactor;
    d_rep3 = 0.000138 * scalingFactor;
    
    
    # Initial values
    DNA1 = 1; # A colony of cells may have anywhere from 10^4 to 10^8 copies of the genome
    DNA2 = 1;
    COMPOUND = 100; 
    DNA_REP1 = 1;
    DNA_REP2 = 1;
    DNA_REP3 = 1;
    
    RNA1 = 0;
    RNA2 = 0;
    NB1 = 0;
    NB2 = 0;
    NBC1 = 0;
    NBC2 = 0;
    DIMER = 0;
    
    RNA_REP1 = 0;
    RNA_REP2 = 0;
    RNA_REP3 = 0;
    
    REP1 = 0;
    REP2 = 0;
    REP3 = 0;
    
""")
r.reset()
# r.draw(width=800,height=300,overlap = "false", splines = "true")
result = r.simulate(0,240, 1000)
r.plot()

plt.figure(1)
plt.plot(result[:,0], result[:,6:8])
plt.xlabel("Time (hours)")
plt.ylabel("Copies")
plt.legend(["Complexes", "Dimers"])


plt.figure(2)
plt.plot(result[:,0],result[:,7])
plt.xlabel("Time (hours)")
plt.ylabel("Dimers, copies")