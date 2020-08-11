from Pipeline.Concentration.SIRVsuite_concentration import *
from Pipeline.Correlation.ERCC_correlation import *
from Pipeline.Coverage.SIRVsuite_coverage_tool * 

groups = ["c1","c1","c2","c2"]
files = ["/home/tdrozd/TV_estimation/NGS2.54/NGS2.54_0001/transcripts_summary_Aligned.sortedByCoord.out.dat",
"/home/tdrozd/TV_estimation/NGS2.54/NGS2.54_0002/transcripts_summary_Aligned.sortedByCoord.out.dat",
"/home/tdrozd/TV_estimation/NGS2.54/NGS2.54_0003/transcripts_summary_Aligned.sortedByCoord.out.dat",
"/home/tdrozd/TV_estimation/NGS2.54/NGS2.54_0004/transcripts_summary_Aligned.sortedByCoord.out.dat"]
sample_names = ["NGS2.54_0001","NGS2.54_0002","NGS2.54_0003","NGS2.54_0004"]
types = ["MIX2","MIX2","MIX2","MIX2"]

## SIRV concentration graphics
groups_unique = np.unique(groups)
relative_abundance, SIRV_transcripts = get_relative_abundance(files, groups_unique)
boxplot = create_sirvsuite_boxplot(relative_abundance, groups = groups_unique)
heatmap = create_sirvsuite_heatmap(relative_abundance, groups = groups_unique, transcript_names = SIRV_transcripts)

## ERCC correlation
ERCC_correlation(files, types, sample_names, "Mix1")

## TODO: SIRV coverage