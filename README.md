# RapidEvolutionNTP

This code simulates evolutionary dynamics at a two-loci, plasticity modifier and a structural locus, in a subdivided population (deme 1 and deme 2) over a patchy (spatially heterogenous) habitat, forward in time as described in Promy et al. 2023. <br> <br>
**The program requires an input file “pars1” which contains a single entry (separated by enter) for each of:**<br>
__________________________________
N1 - size of deme1 = 10000 <br>
N2 - size of deme2 = 10000<br>
TIME - simulation time (100N where N = N1 + N2) = 2000000 <br>
REP - number of replicates (100N)  = 2000000 <br>
interval - interval between migration occurs in generation = 1 <br>
NMU - mutation rate = 0.1 <br>
MIGS - migrants between the demes = 1.0 <br>
sel - selective pressure = 0.01 <br>
REC - recombination rate = 0.5 <br>
pla - relative strength of plasticity, 0-1 = 1.0 <br>
sigma - randomness of selective pressure = 0.0 <br>
 <br>
**Output for the runs is given in “out2” and reports:** <br>
__________________________________
size of deme 1 & 2;<br> migrants;<br> selective pressure in the environment of deme 1 & 2;<br> relative strength of plastisity in ancestral and derived allele of deme 1 & 2;<br> heterozygosity of target & plastic loci;<br> average fixation time of target & plastic loci;<br> average loss time of target & plastic loci; <br> probability of fixation of target & plastic loci.<br>
