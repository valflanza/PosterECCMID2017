# PosterECCMID2017
# MDST (Multi Dimensional Sequence Typing): an Universal Sequence Typing

## Val F Lanza, Fernando Baquero, Teresa M.Coque. 

### Background
One of the major tasks in microbial epidemiology is
the need to type the isolates to allow tracking the outbreaks in
different locations. Up to now the most extensively sequence-based tool
is Multi Locus Sequence Typing (MLST). MLST have some intrinsic
weakness: first, It is a not a correlative scheme (the designation of a
particular ST say nothing about its phylogenetic relation).  Second,
MLST is not always enough precise to define the real clones, and might
overestimate the clone concept. Third, two identical MLSTs in strains of
the same species can include a set of isolates with different
phylogenetic backgrounds. On the other hand, new sequence types are
constantly appearing with the expansion of NGS technology. Approaches as
cgMLST or wgMLST increase the typing accuracy but they do not provide
the possibility of any possible designation to closely phylogenetically
related sequence types. MDST is a new approach that creates a new
continuous scheme classification detect phylogenetically relations
between clones and is applicable to all bacterial organisms.

### Material/Methods
MDST uses a classification based on genetic distance.
Taking a set of sequences  as origin of coordinates  (usually 4 to 6
sequences) MDST calculates the distance from this coordinates to all
other sequences. The distance among sequences are calculated by MASH
software (B. Ondov et al. 2016) providing a normalized way to calculate
phylogenetic distances. About 4000 E. coli strains are classified
respect 6 coordinate sequences is just a few seconds. MDST toolkit
includes an automatic system to select the optimal coordinate sequences,
based on the available sequences. The MDST provide a numerical
classification (coordinates) that allow to represent ensembles of STs
linked by close phylogenetic distances. Particular STs are visualized as
part of these ensembles, easily recognized in a colour scheme in the
coordinates landscape.

### Results: 
MDST has been applied to three different
species: Escherichia coli, Enterococcus faecium and Klebsiella
pneumoniae. MDST grouping and conventional phylogenetic population
structure were closely correlated. Besides, MDST easily identify
sequences that do not correspond with the species classification. MDST
grouping system allows to identify  find the neighbors strains to any
particular clone of interest just calculating the Euclidean distance
among them. Finally, MDST nomenclature allow performing population
structure analysis without the necessity to have all the genomic
sequence data. 

### Conclusions
In the new genomic age, an universal
designation system at supra-specific, specific, and sub-specific levels
is mandatory. The old systems as MLST are not enough precise to
represent the genomic data provided for the NGS and the new
classification system as cgMLST do not have a readable system of
designation for significant ensembles of strains. Other classification
systems such as BAPS (Korander et al 2007) or eBURST (Feil et al 2004)
are dependent of the dataset to analyse.  MDST solve these problems with
a universal system able to classify isolates independently of the
dataset, with a readable classification system based on phylogenetical
relations. 
