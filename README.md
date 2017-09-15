# Act-L-Act-network
Release the codes for network analysis (metrics calculation, randomized network and visualization), activation map files

## Network metrics ##
*   Network efficiency (based on shortest path), and it calls brant (https://github.com/kbxu/brant)
*   Modularity and fractional modularity, adapted from generalized Louvain method by Mucha et al, https://github.com/GenLouvain/GenLouvain
*   Interaction strength, sum-up the connectivty with ends separately locating in the two ROIs
*   Randomized the network for comparing/normalizing modularity index, small-world, et al (call BCT, https://sites.google.com/site/bctnet).

## Matlab plotting ##
*   boxplot overlaid by scattering plot (by notBoxPlot, https://github.com/diazlab/scell/blob/master/mfiles/notBoxPlot.m)
*   Combine the Matalb bar plot with the error plot (by https://github.com/FNNDSC/matlab/blob/master/misc/barwitherr.m)

## Activation map generation ##
*   By FSL/FEAT, input the parameter configuration files
*   Activation map for the 7 tasks in HCP data (slect 30 subjects).

