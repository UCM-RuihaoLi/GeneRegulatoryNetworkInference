# An evolutionary algorithm for gene regulatory network inference
## Requirements:

Tested using

Python 3.6.4<br>
pandas (0.22.0)<br>
numpy (1.15.4)<br>
scipy (1.5.2)<br>

## Usage:
```
Usage: EGRNM.py [-h] -r <input_RNAseq> [-c <input_ChIP>] [-i <iteration_num>]
                                   -t <promoter_strengths> [-s <samplesize>] [-n <training_count>]
                                   -p <PerturbationPower> [-m <mRNA_elongation_rate>] [-a <aminoacid_elongation_rate>]
                                   -l <gene_length> -o <output_name> -[d <protein_degradation>]
```
### Example:
```
python EGRN_Multi_Genalg.py -r data/Alpha_Attractors_example.txt -n 250 -i 3 -t data/Alpha_PromoterStrength_example.txt -l data/Alpha_GeneLength_example.txt -o GRN_Test -s 10
```

## How to construct `input_RNAseq.txt`
The `input_RNAseq.txt` file is tab separated. The first row is a comma-separated list of indices for deleted genes. 
Subsequent rows represent the transcriptional profile.

### Example
Sample 1, a wild type strain:
```
Gene1	Gene2	Gene3	Gene4
5	10	20	30
```
Sample 2, a Gene1 deletion strain:
```
Gene1	Gene2	Gene3	Gene4
0	10	20	30
```
Sample 3, a Gene1&2 double deletion strain:
```
Gene1	Gene2	Gene3	Gene4
0	0	20	30
```

The `input_RNAseq.txt` will then be:
```
-1	5	10	20	30
0	0	10	20	30
0,1	0	0	20	30
```

## How to construct `promoter_strengths.txt`
The `promoter_strengths.txt` file stores the relative promoter strengths for each gene.
Note that the order of the genes and the number of rows should be the same as in `input_RNAseq`.

### Example:
If the relative promoter strengths for Gene1-4 are 100, 200, 300, and 400, 
and they are the same for the two attractors in input_RNAseq, the `promoter_strength.txt` file will read:
```
100	200	300	400
100	200	300	400
```
## How to construct `gene_length.txt`
The `gene_length.txt` file contains the length (in base pair) of each gene.
Note that the order of the genes and the number of rows should be the same as in `input_RNAseq`.

### Example:
If Gene1 is 100 bp long, Gene2 and 3 are 200 bp and Gene4 is 500 bp, the `gene_length.txt` file will read:
```
100 200 200 500
```

## How to construct `input_ChIP.txt` (optional)
In `input_ChIP.txt`, `1` represents positive while `0` represents negative for ChIP signal.
Note that the order of the genes should be the same as in `input_RNAseq.txt`.

### Example:
If Gene1 binds to Gene2 and Gene4 binds to Gene3, the `input_ChIP.txt` file will read:
```
0	1	0	0
0	0	0	0
0	0	0	0
0	0	1	0
```
## Additional Notes
Other parameters can be assigned by users. By default, the mRNA elongation rate is $4.8 nt./s$, amino acid elongation rate is $8 aa./s$, 
protein degradation rate is $0.00796/s$. 

The samplesize (by default $100$) determines the number of GRN instances involved in the evolutionary algorithm, 
and the training_count (by default $800$) refers to how many times the evolutionary algorithm will iterate.

The inferred GRN will be stored in `_Result.txt` and `_Result_f0.txt`. 
The `CheckPoint.txt` file stores the parameters for each individual GRN instance. 
The time spent and the decrease of attractor distance on each iteration will be stored in `CheckLoopCounter.txt` and `ResultFlow.txt`, respectively.
