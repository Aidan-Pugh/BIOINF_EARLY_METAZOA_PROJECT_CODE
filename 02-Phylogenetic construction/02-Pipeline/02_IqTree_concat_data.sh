#!/bin/bash
# After running FASconCAT.pl, run IqTree on the resulting .phy file

iqtree2 -s FcC_supermatrix_COI_CYTB_18S.phy -m MFP -B 1000 -nt AUTO -mem 20G

exit