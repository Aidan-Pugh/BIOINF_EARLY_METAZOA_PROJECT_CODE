# All the scripts and settings used to construct the trees
- 01_Align_and_trim.sh: Loops through all fasta files. Aligns the sequences within the fasta file using Mafft (Katoh et al., 2002) and then trims the sequences using TrimAl (Capella-Gutierrez, 2009).

- 02_IqItree_concat_data.sh: Script used to construct trees from the aligned and trimmed datasets using IqTree (Kalyaanamoorthy et al., 2017). This script was used for single gene data and multiple genes concatenated data.

- 03_IqTree_constraint_tree.sh: Script used to construct a constraint tree from the 3-gene dataset, using a previously constructed backbone tree.

Concatenation of fasta files was done using the perl script FASconCAT.pl (Kueck, 2009).




_Capella-Gutierrez, S., Silla-Martinez, J.M. and Gabaldon, T. (2009). trimAl: a tool for 
automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), 
pp.1972–1973_

_Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., von Haeseler, A., Jermiin, L.S.. (2017) 
ModelFinder: Fast model selection for accurate phylogenetic estimates. Nat. Methods, 
14:587-589_

_Katoh, K., Misawa, K., Kuma, K.I. and Miyata, T. (2002). MAFFT: a novel method for rapid 
multiple sequence alignment based on fast Fourier transform. Nucleic Acids Research, 
30(14), pp.3059–3066._

_Kueck P., FASconCAT, Version 1.0, Zool. Forschungsmuseum A. Koenig, Germany,
2009_
