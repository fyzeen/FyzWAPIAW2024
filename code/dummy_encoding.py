import pandas as pd

confounds_path = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/confounds/HCP_D/confounds.csv"
out_path = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/confounds/HCP_D/dummied_confounds.csv"

confounds = pd.read_csv(confounds_path)

site_dummies = pd.get_dummies(confounds["site"], prefix='dummy', dtype=int) 
dummied_confounds = pd.concat([confounds, site_dummies], axis=1)
out = dummied_confounds.drop("site", axis=1)
out = dummied_confounds.drop("dummy_UCLA", axis=1)


out.to_csv(out_path, index=False)