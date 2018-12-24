"""
Question
--------
Which reactions, independent of substrate and product, are present in eco01100, but not in a combination of all non-overview metabolic pathways?

Method
------
- Download pathway description as KGML.
- Convert to substance-reaction graph.
- Get set of reactions for each graph.
- Calculate difference of reaction sets.
- Print reactions (reaction ID).

Result
------

::

    161 results
    R00107
    R00149
    R00160
    R00173
    R00185
    R00206
    R00209
    R00229
    R00243
    R00257
    R00287
    R00305
    R00320
    R00346
    R00352
    R00369
    R00396
    R00400
    R00414
    R00431
    R00432
    R00502
    R00530
    R00531
    R00550
    R00616
    R00618
    R00648
    R00688
    R00689
    R00692
    R00703
    R00709
    R00710
    R00717
    R00726
    R00727
    R00729
    R00746
    R00841
    R00859
    R00861
    R00866
    R00883
    R00895
    R00907
    R00943
    R00951
    R00978
    R01022
    R01039
    R01063
    R01067
    R01088
    R01122
    R01128
    R01151
    R01156
    R01171
    R01176
    R01217
    R01218
    R01225
    R01278
    R01287
    R01290
    R01324
    R01326
    R01415
    R01434
    R01447
    R01523
    R01555
    R01667
    R01730
    R01741
    R01769
    R01788
    R01829
    R01867
    R01869
    R01911
    R01967
    R01986
    R02014
    R02020
    R02022
    R02073
    R02089
    R02107
    R02124
    R02189
    R02197
    R02239
    R02282
    R02296
    R02415
    R02474
    R02494
    R02537
    R02596
    R02727
    R02821
    R03222
    R03295
    R03354
    R03438
    R03776
    R03856
    R03919
    R03989
    R04007
    R04027
    R04230
    R04413
    R04432
    R04753
    R04859
    R04985
    R05000
    R05053
    R05081
    R05198
    R05705
    R06613
    R06740
    R06975
    R06985
    R07064
    R07147
    R07276
    R07279
    R07376
    R07407
    R07410
    R07818
    R08211
    R08379
    R08557
    R08558
    R08574
    R08714
    R08768
    R08769
    R08773
    R08774
    R08775
    R08781
    R09084
    R09085
    R09127
    R09254
    R09281
    R09944
    R10161
    R10178
    R10305
    R10343
    R10404
    R10699
    R11308

Conclusion
----------
The overview map eco01100 contains reactions that are not present in any of the single metabolism pathways.

Further research shows the reasons for such strange results:

- There seems to have been some mistake in creating eco01100, because some genes are assigned reactions which actually belong to a different enzyme function.
  Example: R10699. Eco01100 maps this to the gene eco:b0774, which is classified as EC:2.6.1.62, which actually only performs R03231.
  Under this reaction, the gene occurs in pathway 00780 and 01100, the latter just adds rn:R10699, for an unknown reason.
  R10699 actually performs a different enzyme function: 2.6.1.105. This function does not occur in any enzyme associated with an eco gene.

- Some single pathway files in KGML format do not contain a "reaction tag", even though they should, according to the "entry tags" present in the same file and the pathway image.
  Example: R05000 only exists in eco011000, but not in eco00130, where it should, because eco00130 contains all compound and gene entries necessary for this reaction.
  It simply does not contain the "reaction tag", probably due to an error in KEGG's algorithms for creating KGML files.
    
- Maybe more.

The overview map 01100 is unreliable and should, thus, not be used for exact research of reactions. However, using only the single pathways can also be unreliable, if the graph construction was to rely upon the "reaction tag".

One solution could be to combine all single pathways plus the overview pathway 01100.
But this adds reactions unknown to any single pathway. These reactions could possibly represent enzymes not present in any single pathway, but they might not. 
In case no additional enzymes occur in 01100, this solution would at least work on the substance-gene level and onward.

Another solution could be not to rely upon "reaction tags", which represent pre-computed graph data. 
Then, however, we would have to compute the graph data ourselves, requiring a constant factor of more downloads (EC or KO files, plus reaction files) from KEGG (factor of roughly 10) and more computation time needed. 
"""
from FEV_KEGG.Graph.SubstanceGraphs import SubstanceReactionGraph
from FEV_KEGG.KEGG.Organism import Organism


if __name__ == '__main__':
    
    #- Download pathway description as KGML.
    eco = Organism('eco')
    
    eco01100 = eco.getPathway('01100')
    allNonOverviewPathways = eco.getMetabolicPathways(includeOverviewMaps = False)
    
    #- Convert to substance-reaction graph.
    eco01100_reactionGraph = SubstanceReactionGraph.fromPathway(eco01100)
    allNonOverviewPathways_reactionGraph = SubstanceReactionGraph.fromPathway(allNonOverviewPathways)
    
    #- Get set of reactions for each graph.
    eco01100_reactions = eco01100_reactionGraph.getReactions()
    allNonOverviewPathways_reactions = allNonOverviewPathways_reactionGraph.getReactions()
    
    #- Calculate difference of reaction sets.
    difference_reactions = eco01100_reactions.difference(allNonOverviewPathways_reactions)
    
    #- Print reactions (substrate -> reaction ID -> product).
    output = []
    for reaction in difference_reactions:
        output.append(reaction.__str__())
    
    output.sort()
    print(str(len(output)) + ' results')
    for line in output:
        print(line)
