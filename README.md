# Host Profit Maximization: Leveraging Performance Incentives and User Flexibility (VLDB'24)


## Introduction
1. This repository contains the full version of our paper. 
2. This repository contains the codes and datasets used in our paper.
3. **Host Profit Maximization: Leveraging Performance Incentives and User Flexibility**.

Abstract: The social network host has knowledge of the network structure and user characteristics and can earn a profit by providing merchants with viral marketing campaigns. We investigate the problem of  host profit maximization by leveraging performance incentives and user flexibility. To incentivize the host's performance, we propose setting a desired influence threshold that would allow the host to receive full payment, with the possibility of a small bonus for exceeding the threshold. Unlike existing works that assume a user's choice is frozen once they are activated, we introduce the  Dynamic State Switching model to capture  ``comparative shopping'' behavior from an economic perspective, in which users have the flexibilities to change their minds about which product to adopt based on the accumulated influence and propaganda strength of each product. In addition, the incentivized cost of a user serving as an influence source is treated as a negative part of the host's profit.

The Host Profit Maximization problem is NP-hard, submodular, and non-monotone. To address this challenge, we propose an efficient greedy algorithm and devise a scalable version with an approximation guarantee to select the seed sets. As a side contribution, we develop two seed allocation algorithms to balance the distribution of adoptions among merchants with small profit sacrifice. Through extensive experiments on four real-world social networks, we demonstrate that our methods are effective and scalable.


## Datasets
We use four publicly available real-world road networks, including NetHEPT, Epinions, DBLP and LiveJournal datasets.
1. NetHEPT can be obtained from [1].
2. Epinions, DBLP and LiveJournal can be obtained from [2].

[1] Wei Chen, YajunWang, and Siyu Yang. 2009. Efficient influence maximization in social networks. In Proceedings of the 15th ACM SIGKDD international conference on Knowledge discovery and data mining. 199–208.

[2] Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford large network dataset collection. http://snap.stanford.edu/data.


## Algorithms
The following files are the codes for our proposed algorithms. We implemented all codes the using C++ with CLion 2021.2.3.

1. First use **genCost** (can be found in Multi-ProMax/genCost directory) to generate propagation probability and cost, specifically,
+ We use the **Weighted-Cascade model** to generate propagation probability of each edge in social graph,
+ We use **Degree-Proportional Cost Model** to generate cost of each node in social graph.

Configuration: `"dataset model type epsilon"`

2. Then use **alg** (can be found in Multi-ProMax/alg directory) to tackle our problem, **alg** includes:
+ **ROI-Greedy**<sup>[3]</sup> for a single merchant,
+ **Fill**, **OBO**, **ITER** for multiple merchants.

Configuration: `"-dataset xxx -epsilon xxx  -model LT -time xxx -theta xxx -h xxx -scale xxx -type xxx -gammaR xxx -gammaP xxx -algo xxx -batchPS xxx -batchIS xxx"`


[3] Tianyuan Jin, Yu Yang, Renchi Yang, Jieming Shi, Keke Huang, and Xiaokui Xiao. 2021. Unconstrained submodular maximization with modular costs: Tight approximation and application to profit maximization. Proceedings of the VLDB Endowment 14, 10 (2021), 1756–1768.




## Running Environment
A 64-bit Linux-based OS.

Note that the dataset directory should be placed in the upper two layers of the alg directory.

