# Systematic comparison of rankingaggregation methods for gene lists inexperimental results

Code for simulated data generation and the way to get and use investigated algorithms for the study: Systematic comparison of rankingaggregation methods for gene lists inexperimental results.

# Authors:
- Bo Wang
- Michael Gutmann
- Kenneth Baillie

# data generator
A common experimental output in biomedical science is a list of genes implicated in a given biologicalprocess or disease. 
The results of a group of studies answering the same, or similar, questions can becombined by meta-analysis, aiming at 
finding a more reliable answer. 
Ranking aggregation methods canbe used for combining gene lists from various sources in meta-analyses. 
Evaluating a ranking aggregationmethod on a specific type of dataset before using it can support the reliability of the 
result since theproperty of a dataset can influence the performance of an algorithm. 
Evaluation of aggregation methodsis usually based on a simulated database because of the lack of a known truth for real data. 
Simulated datasets tend to be too small and neglect key features of experimental data, including heterogeneity ofquality, 
relevance and the inclusion of unranked lists. 
This script can generate simulated data to explore the performance of the aggregation methods as a function of emulating
the common scenarios of real genomics data, with various heterogeneity of quality, noise level, and a mix of unranked and 
ranked data using 20000 possible entities.  
It is based on the analysis of real genomic data and the simulated data generation model in the study of MAIC algorithm 
that samples a score from specific distribution for each entity to simulate the rankings. The investigated real data
include SARS-COV-2virus, cancer(NSCLC), and bacteria(macrophage apoptosis).

The model ranks entities by figure Z for each entity generated from Gaussian distribution.
For entity k in experiment(list) i,
Z_ki ~ N (mu_k, sigma_i^2)
log(sigma_i)~N(log(mean_noise_M), heterogeneity_D)
The length of each list is also sampled from specific distribution to emulating real datasets.

The generated data are lists and will be written in a file with name "(mean_noise)_(heterogeneity_D)_(dataset_type(
if unranked lists included)).txt"

Examples of running are shown at bottom of this file below "if __name__ == '__main__':" within data_generator.py

# ranking aggregation algorithms
Files in the algorithms folder are for running existing methods, including the algorithm itself, scripts to control the running of each algorithm, and small scripts for existing algorithms provided by related authors. Only RepeatChoice is newly implemented and only rMixMED, rMixMEAN,
rMixGEO, MixStuart are editted from R package:RobustRankAggregedited to include unranked lists.

Investigated algorithms(where to get them):<br />
1: MAIC: https://github.com/baillielab/maic. (not included here)<br />
2: Vote Counting: vote_counting.py, mentioned in MAIC study by simply ranking entities by frequency.<br />
3: RRA: use R package: RobustRankAggreg. Run RRA.R in this folder for running it.<br />
4: RepeatChoice: repeat_choice.py, implemented in this study.<br />
5: BIRRA: (Badgeley et al., 2015) provided by author.<br />
6-9: rMED, rMEAN, rGEO, stuart: use R package: RobustRankAggreg.<br />
10-13: tMED, tMEAN, tGEO, tL2: R package:TopKLists.<br />
14-17: rMixMED, rMixMEAN, rMixGEO, MixStuart: editted from R package:RobustRankAggregedited to include unranked lists.<br />
18-19: BiGbottom, BiGNA:  (Li et al., 2018) provided by author.<br />
20-22: MC1-3: R package:TopKLists.<br />
23: BARD: (Deng et al., 2014), provided by author. (not included here)
