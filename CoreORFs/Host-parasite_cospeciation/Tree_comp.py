import pickle
from skbio.stats.evolve import hommola_cospeciation
import skbio as skb
import numpy as np
import pandas as pd
from skbio import TreeNode
from io import StringIO
import os

Bac_Tree_Dir      = "/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/197_Trees_Rooted_with_OutGroup/"
host_tree_fp      = '/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/Host-parasite_cospeciation/Host.nwk'
bac2host_file     = "/cellulose-degrading-human-gut-bacteria-2023/CoreORFs/Host-parasite_cospeciation/Bac2Host"
Tree_File_Pattern = ".treefile"
Out_Dir           = "/"

def correct_name(x):
    return(x.replace(' ','_').replace('.','_').replace('-','_'))

bac2host     = pd.read_table(bac2host_file)
host_tree    = TreeNode.read(host_tree_fp)
# host_tree.assign_ids()


for tip in host_tree.tips():
    tip.name = correct_name(tip.name)
    
bac2host = bac2host.map(correct_name)

Host_list = {x: 0 for x in bac2host["Host"].unique()}
incidence_dict = {}
for bac in bac2host["Bac"]:
    Host_list_copy = Host_list.copy()
    Host = bac2host.loc[bac2host["Bac"] == bac,"Host"].values[0]
    Host_list_copy[Host] =1 
    incidence_dict[bac] = Host_list_copy

incidence_table = pd.DataFrame.from_dict(incidence_dict,
                                         orient='index')

host_dists = host_tree.tip_tip_distances()
host_dists = host_dists.filter(bac2host["Host"].unique())


reasults = pd.DataFrame()

for file in os.listdir(Bac_Tree_Dir):
    if file.endswith(Tree_File_Pattern):
        print(file)
        bact_tree_fp = os.path.join(Bac_Tree_Dir, file)
        bact_tree = TreeNode.read(bact_tree_fp)
        for tip in bact_tree.tips():
            tip.name = correct_name(tip.name)
        bac_dists  = bact_tree.tip_tip_distances()
        bac_dists  = bac_dists.filter(bac2host["Bac"].values)
        r, p, d = hommola_cospeciation(host_dists,bac_dists, incidence_table,100000)
        reasults.loc[file,"R Correlation"] = r
        reasults.loc[file,"P Value"] = p
        print("Pearson correlation coefficient of host : parasite association:")
        print(r)
        print(" Significance of host : parasite association computed using permutations and a one-sided (greater) alternative hypothesis:")
        print(p)
reasults.to_csv( os.path.join(Out_Dir,"Mantel_Test_Results.csv"))
