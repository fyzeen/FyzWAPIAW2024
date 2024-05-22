import pandas as pd

subj_list = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/subjects_list_aging.csv", header=None)
confounds = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/confounds_num_aging.csv")
phen = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/phenotypes_aging.csv")
dkt = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/dkt_yeo_aging.csv")


for column in dkt.columns:
    print(dkt[column].isnull().values.any())
