# Exploring batch effects in the LINCS L1000 data

Good afternoon everyone, in this notebook tutorial we're going to explore the batch effects in different levels of the LINCS L1000 data and introduce a few computational pipelines adjust the batch effect. 

We will first get us familiar with the Level 3 data and explore its batch effect.
Next, we will look at the Level 5 data and learn some approaches for visualizing and quantifying the batch effects. 
We will also look at Level 5 signatures generated from different computational pipelines and benchmark these batch effect correction methods.
Finally, we will look into the principles and applications of interactive visualizations of the Level 5 signatures. 

# 1. Explore batch effects in the Level 3 data

As we know, the Level 3 data contains the gene expression profiles of samples from the 384-well plates of the L1000 experiments. It is normalized to the invariant probesets and quantile normalized. It is similar to the normalized expression matrix from RNA-seq or microarrays where one can perform DE analysis to identify gene expression signature for a biological condition versus another.

## 1.2. Examine different batch effects in the Level 3 data with regards to:

Hierarchies of batch factors

- Perturbagen plate (`pert_plate`, example CPC005)
- RNA plate (`rna_plate`, example CPC005_A549_24H_X1)
- Detection plate (`det_plate`, example CPC005_A549_24H_X1_B20)

To do that, we perform PCA on the data and color the sample by their different levels of batch factors.

We can see the RNA plate and detection plate seem to be the source of the batch effect for this particular subset of data.

## 1.3. Mitigate the batch effects by mean centering within the batches

We can try a very simple way to correct for this effect, which is to use the means from each batch to center the expression matrix.

# 2. Batch effects in the Level 5 data

We will compare the Level 5 data generated from different computational pipelines and assess their batch effects. These pipeline include ...
First we will parse the signature-level metadata to select the drug/compound treatment signatures from PC3. We removed signatures generated from only one treatment experiment.

# 2.2. Examination of batch effects in the Level 5 data (signatures)
We'll use two methods to inspect the batch effects.
The intuiation is to perform unsupervised DR to visualize the signatures, and examine the clustering patterns with regards to the batch factors and drugs. 

Higher score suggests more agreements between the clustering of signatures with the drugs. 

MOAs of top 10 drugs in Level 5 signatures
- HDAC inhibitors: trichostatin-a, vorinostat, panobinostat, tubastatin-a, PCI-34051
- wortmannin: PI3K inhibitor
- geldanamycin: HSP90 inhibotor
- LY-294002: mTOR/PI3K inhibitor
- fulvestrant: ER antagonist
- troglitazone: insulin sensitizer, PPAR receptor agonist
- sirolimus: mTOR inhibitor

### 2.2.2. 
Connectivity score is used to quantify the similarity beteen a pair of gene expression signature. 
We would expect the connectivity score for signatures from the same drugs would be high. 
Idealy it should be irrelevant to the batch of the signatures.

# 3. 
The goal of interactive visualization is the allow users to overlap multiple layers of information onto the plot and being able to explore. 

Next, we will first make interactive scatter plot of the signatures.

