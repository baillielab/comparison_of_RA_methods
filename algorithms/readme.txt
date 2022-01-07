Files in the algorithms folder are for running existing methods, including the algorithm itself, scripts to control the running of each algorithm, where small scripts within ./algorithms for existing algorithms BIRRA and BiG are provided by related authors. Only RepeatChoice and VC (Vote Counting) are newly implemented following the description of related papers. rMixMED, rMixMEAN,
rMixGEO, MixStuart are editted based on R package:[RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg) to include unranked lists.

Investigated algorithms(where to get them):<br />
1: MAIC: [(Li et al., 2020)](https://www.nature.com/articles/s41467-019-13965-x) code: https://github.com/baillielab/maic. <br />
2: Vote Counting: mentioned in [(Li et al., 2020)](https://www.nature.com/articles/s41467-019-13965-x) by simply ranking entities by frequency. code: ./vote_counting.py <br />
3: RRA: use R package: [RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg). Run ./RRA.R in this folder for running it. <br />
4: RepeatChoice: [(Ailon, 2010)](https://link.springer.com/article/10.1007/s00453-008-9211-1). code: ./repeat_choice.py, implemented in this study.<br />
5: BIRRA: [(Badgeley et al., 2015)](https://doi.org/10.1093/bioinformatics/btu518) code: provided by author (http://www.pitt.edu/~mchikina/BIRRA/) and included here (./runBIRRA.R), using ./useBIRRA.R to run it.<br />
6-9: rMED, rMEAN, rGEO, stuart: use R package: [RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg).<br />
10-13: tMED, tMEAN, tGEO, tL2: use R package:[TopKLists](https://CRAN.R-project.org/package=TopKLists).<br />
14-17: rMixMED, rMixMEAN, rMixGEO, MixStuart: editted from R package:[RobustRankAggreg](https://CRAN.R-project.org/package=RobustRankAggreg) to accept unranked lists as input.<br />
18-19: BiGbottom, BiGNA: [(Li et al., 2018)](https://doi.org/10.1002/sim.7920) coed: kindly provided by the authors (https://github.com/xuelilyli/BiG) and also included here (./BiG_code_platform_changed.R).<br />
20-22: MC1-3: R package: [TopKLists](https://CRAN.R-project.org/package=TopKLists).<br />
23: BARD: [(Deng et al., 2014)](https://doi.org/10.1080/01621459.2013.878660), kindly provided by the authors. 
