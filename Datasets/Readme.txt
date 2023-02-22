In each GRN_x folder, there are Attractors_x.txt, GeneLength_x.txt, and PromoterStrength_x.txt.

Attractors_x.txt contains the input transcriptional profiles. The first column is the indicator which specifies which gene has been knocked out or overexpressed.
If the indicator is -1, no gene has been knocked out or overexpressed.
If the ith gene has been knocked out, the indicator will be (-1 + i). For example, if the 1st gene has been knocked out, the indicator will be 0.
If the ith gene has been overexpressed, the indicator will be (-1 - i). For example, if the 2nd gene has been overexpressed, the indicator will be -3.

GeneLength_x.txt contains the lengths of the genes, in terms of base pairs.

PromoterStrength_x.txt contains the promoter strengths for each gene in each transcriptional profile. If the Attractors_x.txt has n transcriptional profiles, the PromoterStrength_x.txt should also have n rows.