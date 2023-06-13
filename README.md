# RapidEvolutionNTP

This code simulates evolutionary dynamics at a two-loci, plasticity modifier and a structural locus, in a subdivided population (deme 1 and deme 2) over a patchy (spatially heterogenous) habitat, forward in time as described in Promy et al 2023. The program requires an input file “pars1” which contains a single entry (separated by enter) for each of:<br>
N1 – size of deme1 <br>
N2 - size of deme2 <br>
TIME - simulation time (100N where N = N1 + N2) <br>
REP - number of replicates (100N) <br>
interval - interval between migration occurs in generation <br>
NMU - mutation rate <br>
MIGS - migrants between the demes <br>
sel - selective pressure <br>
REC - recombination rate <br>
pla - relative strength of plasticity, 0-1 <br>
sigma - randomness of selective pressure <br>
A sample pars1 is provided.<br>
Output for the runs is given in “out2” and reports: <br>
size of deme 1 & 2,<br> migrants,<br> selective pressure in the environment of deme 1 & 2,<br> relative strength of plastisity in ancestral and derived allele of deme 1 & 2,<br> heterozygosity of target & plastic loci, average fixation time of target & plastic loci,<br> average loss time of target & plastic loci, <br> probability of fixation time of target & plastic loci.<br>
