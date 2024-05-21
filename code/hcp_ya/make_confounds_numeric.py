import pandas as pd

confounds_path = "/scratch/l.lexi/WAPIAW2024/Data/HCP_YA/confounds.csv"
out_path = "/scratch/l.lexi/WAPIAW2024/Data/HCP_YA/dummied_confounds.csv"

confounds = pd.read_csv(confounds_path)

fam_id_dummies = pd.get_dummies(confounds["Family_ID"], prefix='dummy', dtype=int) 
dummied_confounds = pd.concat([confounds, fam_id_dummies], axis=1)
out = dummied_confounds.drop("Family_ID", axis=1)

out.to_csv(out_path)
