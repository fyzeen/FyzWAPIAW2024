import pandas as pd

confounds_path = "/scratch/l.lexi/WAPIAW2024/Data/HCP_D/confounds.csv"
out_path = "/scratch/l.lexi/WAPIAW2024/Data/HCP_D/dummied_confounds.csv"

confounds = pd.read_csv(confounds_path)

site_dummies = pd.get_dummies(confounds["site"], prefix='dummy', dtype=int) 
dummied_confounds = pd.concat([confounds, site_dummies], axis=1)
out = dummied_confounds.drop("site", axis=1)

out.to_csv(out_path)
