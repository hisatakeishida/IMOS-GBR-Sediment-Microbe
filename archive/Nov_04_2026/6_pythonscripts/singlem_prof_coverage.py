import pandas as pd
import os
from pathlib import Path
import shutil

txt_dir = "INDIR with *prof.tsv"
outpath= "singlem.profsum.csv"

os.chdir(txt_dir)

#######################
# data loading   
#######################

sample_otu = {}
for file in os.listdir():
    with open(file) as f:
        otu_dict = {}
        for line in f.readlines():
            items = line.strip()
            if items.startswith("sample"):
                pass
            else:
                coverage = float(items.split("\t")[1])
                sample_name = items.split("\t")[0]
                taxonomy = (items.split("\t")[-1]) 
                # print(coverage, sample_name, taxonomy)
                otu_dict[taxonomy] = coverage

        sample_otu [sample_name] = otu_dict
        # break
print(sample_otu)

#######################
# get all feature
#######################

otu_list = []
for sample, otu_info in sample_otu.items():
    for otu, info in otu_info.items():
        if otu not in otu_list:
            otu_list.append(otu)

print(len(otu_list))

#######################
# calculate coverage of each feature across all samples. Add 0 if absent
#######################

otu_read_list = []
hit_count = 0
miss_count = 0
for sample, otu_info in sample_otu.items():
    sample_otu_read_dict = {}
    sample_otu_read_dict["0_otu"] = sample
    # print(sample_otu_read_dict)
    for otu in otu_list:
        if otu in otu_info.keys():
            sample_otu_read_dict[otu]=otu_info[otu]
        else:
            sample_otu_read_dict[otu]=0
    sorted_dict = dict(sorted(sample_otu_read_dict.items()))
    otu_read_list.append(sorted_dict)

#######################
# convert into dataframe
#######################

df = pd.DataFrame(data=otu_read_list)

df.to_csv(outpath,index=False)
