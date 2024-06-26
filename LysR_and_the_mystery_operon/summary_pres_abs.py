#!/usr/bin/python3

import pandas as pd

# Read files
dfgamma = pd.read_table("operon_pres_abs/gamma.txt", names=["Species", "Operon"])
dfalpha = pd.read_table("operon_pres_abs/alpha.txt", names=["Species", "Operon"])
dfbeta = pd.read_table("operon_pres_abs/beta.txt", names=["Species", "Operon"])

# Store only strain names
dfgamma["Species"] = dfgamma["Species"].str[0:6]
dfalpha["Species"] = dfalpha["Species"].str[0:6]
dfbeta["Species"] = dfbeta["Species"].str[0:6]

# Group by species and operon presence/absence
group_gamma = dfgamma.groupby(["Species", "Operon"]).size().unstack(fill_value=0)
group_alpha = dfalpha.groupby(["Species", "Operon"]).size().unstack(fill_value=0)
group_beta = dfbeta.groupby(["Species", "Operon"]).size().unstack(fill_value=0)

# Rename the columns
group_gamma.columns = ["Absent", "Present"]
group_alpha.columns = ["Absent", "Present"]
group_beta.columns = ["Absent", "Present"]

# Reset the index to make 'Species' a regular column
group_gamma = group_gamma.reset_index()
group_alpha = group_alpha.reset_index()
group_beta = group_beta.reset_index()

# Add percentage column
group_gamma["%Presence"] = (group_gamma["Present"]*100) / (group_gamma["Present"] + group_gamma["Absent"])
group_alpha["%Presence"] = (group_alpha["Present"]*100) / (group_alpha["Present"] + group_alpha["Absent"])
group_beta["%Presence"] = (group_beta["Present"]*100) / (group_beta["Present"] + group_beta["Absent"])

# Save table
group_gamma.to_csv("operon_pres_abs/gamma_summary.txt", index=False, sep="\t")
group_alpha.to_csv("operon_pres_abs/alpha_summary.txt", index=False, sep="\t")
group_beta.to_csv("operon_pres_abs/beta_summary.txt", index=False, sep="\t")
