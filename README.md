# RapidEvolutionNTP

This code simulates evolutionary dynamics at a two-loci, plasticity modifier and a structural locus, in a subdivided population (deme 1 and deme 2) over a patchy (spatially heterogenous) habitat, forward in time as described in Promy et al 2023. The program requires an input file “pars1” which contains a single entry (separated by enter) for each of:
N1 – size of deme1
N2 - size of deme2
TIME - simulation time (100N where N = N1 + N2)
REP - number of replicates (100N)
interval - interval between migration occurs in generation
NMU - mutation rate
MIGS - migrants between the demes
sel - selective pressure
REC - recombination rate
pla - relative strength of plasticity, 0-1
sigma - randomness of selective pressure
A sample pars1 is provided.
Output for the runs is given in “out2” and reports:
size of deme 1 & 2, migrants, selective pressure in the environment of deme 1 & 2, relative strength of plastisity in ancestral and derived allele of deme 1 & 2, heterozygosity of target & plastic loci, average fixation time of target & plastic loci, average loss time of target & plastic loci, probability of fixation time of target & plastic loci.
