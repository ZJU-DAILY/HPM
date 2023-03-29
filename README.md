# Host Profit Maximization for Multiple Competing Products


## Introduction
1. This repository contains the extended version of our paper. 
2. This repository contains the codes and datasets used in our paper.
3. **Host Profit Maximization for Multiple Competing Products** is a novel problem of host profit maximization for multiple competing products in a social network, where each merchant is willing to pay a budget if a desired level of influence is achieved. In this problem, the incentivized cost of a user serving as an influence source is treated as a negative part of the host's profit.


## Datasets
We use four publicly available real-world road networks, including NetHEPT, Epinions, DBLP and LiveJournal datasets.
1. NetHEPT can be obtained from [1].
2. Epinions, DBLP and LiveJournal can be obtained from [2].

[1] Chen, Wei, Yajun Wang, and Siyu Yang. "Efficient influence maximization in social networks." Proceedings of the 15th ACM SIGKDD international conference on Knowledge discovery and data mining. 2009.

[2] Leskovec, Jure, and Andrej Krevl. "SNAP Datasets: Stanford large network dataset collection." (2014). http://snap.stanford.edu/data.


## Algorithms
The following files are the codes for our proposed algorithms. We implemented all codes the using C++ with CLion 2021.2.3.

1. First use **genCost** (can be found in Multi-ProMax/genCost directory) to generate propagation probability and cost, specifically,
+ We use the **Weighted-Cascade model** to generate propagation probability of each edge in social graph,
+ We use **Degree-Proportional Cost Model** to generate cost of each node in social graph.

Configuration: `"dataset model type epsilon"`

2. Then use **alg** (can be found in Multi-ProMax/alg directory) to tackle our problem, **alg** includes:
+ **ROI-Greedy**<sup>[3]</sup> and **Simple-Greedy**<sup>[4]</sup> for a single merchant,
+ **HPM**, **SIM**, **OBO**, **ITER** for multiple merchants.

Configuration: `"-dataset xxx -epsilon xxx  -model LT -time xxx -theta xxx -h xxx -scale xxx -type xxx -gammaR xxx -gammaP xxx -algo xxx -batchPS xxx -batchIS xxx"`


[3] Jin, Tianyuan, et al. "Unconstrained submodular maximization with modular costs: Tight approximation and application to profit maximization." Proceedings of the VLDB Endowment 14.10 (2021): 1756-1768.

[4] Lu, Wei, and Laks VS Lakshmanan. "Profit maximization over social networks." 2012 IEEE 12th International Conference on Data Mining. IEEE, 2012.




## Running Environment
A 64-bit Linux-based OS.

Note that the dataset directory should be placed in the upper two layers of the alg directory.

