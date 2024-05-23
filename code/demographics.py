import pandas as pd

confounds_fpath = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/demographics_data/hcp_ya/confounds.csv"
phenotype_fpath = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/demographics_data/hcp_ya/phenotypes.csv"
subjects_list_fpath = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/demographics_data/hcp_ya/subjects_list.csv"

confounds = pd.read_csv(confounds_fpath)
subject_list = pd.read_csv(subjects_list_fpath, header=None)
phenotype = pd.read_csv(phenotype_fpath)

confound_cleaned = confounds.loc[confounds["Subject"].isin(subject_list[0])]
phenotype_cleaned = phenotype.loc[phenotype["Subject"].isin(subject_list[0])]

print(confound_cleaned["age"].describe())
print(confound_cleaned["sex"].describe())
print(phenotype_cleaned["Neuroticism"].describe())