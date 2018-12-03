# Exploring batch effects in the LINCS L1000 data

Good afternoon everyone, in this notebook tutorial we're going to explore the batch effects in different levels of the LINCS L1000 data and introduce a few computational pipelines adjust the batch effect. We will also learn how to visualize the data. 

We will first get us familiar with the Level 3 data and explore its batch effect.
Next, we will look at the Level 5 data and learn some approches for visualizing and quantifying the batch effects. 
We will also look at Level 5 signatures generated from different computational pipelines and benchmark these batch effect correction methods.
Finally, we will look into principles and applications of interactive visualizations of the Level 5 signatures. 

# 1. Explore batch effects in the Level 3 data

As we know, the Level 3 data contains the gene expression profiles of samples fomr the 384-well plates of the L1000 experiments. It is normalized to the invariant probesets and quantile normalized. It is similar to the normalized expression matrix from RNA-seq or microarrays.

## 1.2. Examine different batch effects in the Level 3 data with regards to:

- Perturbagen plate (`pert_plate`, example CPC005)
- RNA plate (`rna_plate`, example CPC005_A549_24H_X1)
- Detection plate (`det_plate`, example CPC005_A549_24H_X1_B20)

To do that, we perform PCA on the data and color the sample by their different levels of batch factors.

We can see the RNA plate and detection plate seem to be the source of the batch effect for this particular subset of data.

## 1.3. Mitigate the batch effects by mean centering within the batches

We can try a very simple way to correct for this effect by using the mean from each batch to center the expression matrix.

# 2. Batch effects in the Level 5 data

We will compare the Level 5 data generated from different computational pipelines and assess their batch effects. 
First we will 


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


