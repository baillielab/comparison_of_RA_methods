Files in this folder are for running existing methods, including the algorithm itself and
scripts to control the running of each algorithm. Only RepeatChoice is newly implemented and only rMixMED, rMixMEAN,
rMixGEO, MixStuart are editted from R package:RobustRankAggregedited to include unranked lists.

Investigated algorithms(where to get them):
1. MAIC: https://github.com/baillielab/maic.
2. Vote Counting: vote_counting.py, mentioned in MAIC study by simply ranking entities by frequency.
3. RRA: use R package: RobustRankAggreg. Run RRA.R in this folder for running it.
4. RepeatChoice: repeat_choice.py, implemented in this study.
5. BIRRA: (Badgeley et al., 2015) provided by author.
6-9. rMED, rMEAN, rGEO, stuart: use R package: RobustRankAggreg.
10-13: tMED, tMEAN, tGEO, tL2: R package:TopKLists.
14-17: rMixMED, rMixMEAN, rMixGEO, MixStuart: editted from R package:RobustRankAggregedited to include unranked lists.
18-19: BiGbottom, BiGNA:  (Li et al., 2018) provided by author.
20-22: MC1-3: R package:TopKLists.
23: BARD: (Deng et al., 2014), provided by author.
