# Systematic comparison of ranking aggregation methods for gene lists in experimental results

Code for simulated data generation and the way to get and use investigated algorithms for the study: Systematic comparison of ranking aggregation methods for gene lists in experimental results.

# Authors:
- Bo Wang
- Michael U. Gutmann
- J. Kenneth Baillie

# data generator
A common experimental output in biomedical science is a list of genes implicated in a given biological process or disease. The results of a group of studies answering the same, or similar, questions can be combined by meta-analysis to find a consensus or a more reliable answer. Ranking aggregation methods can be used to combine gene lists from various sources in meta-analyses. Evaluating a ranking aggregation method on a specific type of dataset before using it is required to support the reliability of the result since the property of a dataset can influence the performance of an algorithm. Evaluation of aggregation methods is usually based on a simulated database especially for the algorithms designed for gene lists because of the lack of a known truth for real data. However, 
simulated datasets tend to be too small compared to experimental data and neglect key features, including heterogeneity of quality, relevance and the inclusion of unranked lists. In this study, a group of existing methods and their variations which are suitable for meta-analysis of gene lists are compared using simulated and real data. Simulated data was used to explore the performance of the aggregation methods as a function of emulating the common scenarios of real genomics data, with various heterogeneity of quality, noise level, and a mix of unranked and ranked data using 20000 possible entities. In addition to the evaluation with simulated data, a comparison using real genomic data on the SARS-CoV-2 virus, cancer (NSCLC), and bacteria (macrophage apoptosis) was performed. We summarise our evaluation results in terms of a simple flowchart to select a ranking aggregation method for genomics data.

The stochastic generative model in this study ranks entities by figure Z for each entity generated from Gaussian distribution.
For entity k in experiment(list) i,

<img src="https://render.githubusercontent.com/render/math?math=Z_{ki}\sim \N( {\mu}_{k},{\sigma}_{i}^{2} )">

<img src="https://render.githubusercontent.com/render/math?math=\ln({\sigma}_{i} )\sim \N(\ln(M), D)">

where M is the mean noise level and D is the heterogeneity of quality for sources.

The length of each list is also sampled from specific distribution to emulating real datasets.

The generated data are lists and will be written in a file with name "(mean_noise M)\_(heterogeneity D)\_(dataset_type S).txt"

Examples of running are shown at bottom of this file below "if __name__ == '__main__':" within data_generator.py

# ranking aggregation algorithms
Files in the algorithms folder are for running existing methods, including the algorithm itself, scripts to control the running of each algorithm, where small scripts within ./algorithms for existing algorithms BIRRA and BiG are provided by related authors. Only RepeatChoice and VC (Vote Counting) are newly implemented following the description of related papers. rMixMED, rMixMEAN,
rMixGEO, MixStuart are editted based on R package:[RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg) to include unranked lists.

Investigated algorithms(where to get them):<br />
1: MAIC: [(Li et al., 2020)](https://www.nature.com/articles/s41467-019-13965-x) code: https://github.com/baillielab/maic. <br />
2: Vote Counting: mentioned in [(Li et al., 2020)](https://www.nature.com/articles/s41467-019-13965-x) by simply ranking entities by frequency. code: ./algorithms/vote_counting.py <br />
3: RRA: use R package: [RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg). Run ./algorithms/RRA.R in this folder for running it. <br />
4: RepeatChoice: [(Ailon, 2010)](https://link.springer.com/article/10.1007/s00453-008-9211-1). code: ./algorithms/repeat_choice.py, implemented in this study.<br />
5: BIRRA: [(Badgeley et al., 2015)](https://doi.org/10.1093/bioinformatics/btu518) code: provided by author (http://www.pitt.edu/~mchikina/BIRRA/) and included here (./algorithms/runBIRRA.R), using ./algorithms/useBIRRA.R to run it.<br />
6-9: rMED, rMEAN, rGEO, stuart: use R package: [RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg).<br />
10-13: tMED, tMEAN, tGEO, tL2: use R package:[TopKLists](https://CRAN.R-project.org/package=TopKLists).<br />
14-17: rMixMED, rMixMEAN, rMixGEO, MixStuart: editted from R package:[RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg) to accept unranked lists as input.<br />
18-19: BiGbottom, BiGNA: [(Li et al., 2018)](https://doi.org/10.1002/sim.7920) code: kindly provided by the authors(https://github.com/xuelilyli/BiG) and also included here (./algorithms/BiG_code_platform_changed.R), using ./algorithms/useBiG.R to run it. <br />
20-22: MC1-3: R package: [TopKLists](https://CRAN.R-project.org/package=TopKLists).<br />
23: BARD: [(Deng et al., 2014)](https://doi.org/10.1080/01621459.2013.878660), kindly provided by the authors. 
