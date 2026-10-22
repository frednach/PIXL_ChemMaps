### import necessays libraries
import pandas as pd
import py_functions.pixl_edgeDetection as ED
import matplotlib.pyplot as plt
import numpy as np
import py_functions.mineralogy as mineral
import mpltern
from pathlib import Path


#import formatted PIXL data from parquet into DataFrame and collected targets and number of pixels
pixl_data = pd.read_parquet("pixl_cache/processed_data/master_merged.parquet")
target_counts_df = pd.read_csv("pixl_cache/targets/targets_points.csv")

#pick relevant oxides to plot and look for and find similarities in with the edgeDetection algorithm
#Has to include: SiO2,MgO,FeO-T
oxide_columns = [
    "SiO2_wt%",
    "Al2O3_wt%",
    "FeO-T_wt%",
    "MgO_wt%",
    "CaO_wt%",
    "Na2O_wt%",
    "K2O_wt%",
    "TiO2_wt%",
    "SO3_wt%"
]

#define minimum domain size (spot size = 120µm)
minimum_size = 5


all_domains = pd.DataFrame()

#start loop over all relevant locations
# grouping by target name and sol to separate different observations
for i, ((target_name, sol), target_data) in enumerate(
    pixl_data.groupby(["TARGET", "SOL"])
):
    #"best" threshold
    threshold = 30
    #if target_name == "Dourbes" or "Cheyava Falls":

    print(f"Processing {target_name}, sol={sol}")        

    domain_summary, target_data = ED.pixl_edgeDetection(
        pixl_data, 
        oxide_columns, 
        target_name,
        sol, 
        True, 
        False, 
        5, 
        threshold
    )

    #check if the target got skipped i.e. is not a DataFrame
    if not isinstance(domain_summary, pd.DataFrame):
        continue

    #plot ternary composition graphs and write results to domain_summary
    fig, tax1, tax2 = mineral.plot_both_ternary(
    domain_summary,
    threshold,
    target_name="Dourbes"
    )
    plt.savefig(
        Path(f"pixl_cache/target_analysis/{target_name}_{sol}/{target_name}_Ternary_T{threshold}.svg"),
        format="svg"
    )
    #plt.show()
    plt.close()

    #add a column containing the target name and one conatining the SOL for referencing
    domain_summary["target"] = target_name
    domain_summary["SOL"] = target_data["SOL"].iloc[0]

    #initialize dataFrame on first iteration
    if i == 0:
        all_domains = domain_summary
    #combine all domains in one global file
    else:
        all_domains = pd.concat([all_domains, domain_summary], axis=0, ignore_index=True)



    #create Harker Diagram of MgO/SiO2
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )
    #
    ax.errorbar(
        domain_summary["SiO2_wt%"],
        domain_summary["MgO_wt%"],
        xerr=domain_summary["SiO2_wt%_std"],
        yerr=domain_summary["MgO_wt%_std"],
        fmt="o",
        markersize = 8
    )
    ax.set_xlabel("SiO2 [wt%]")
    ax.set_ylabel("MgO [wt%]")
    ax.set_aspect("equal")
    ax.set_title(
        f"{target_name} Harker Diagram T={threshold}"
    )
    plt.savefig(
                Path(f"pixl_cache/target_analysis/{target_name}_{sol}/{target_name}_Harker_SiO2_MgO_T{threshold}.svg"),
                format="svg"
    )
    plt.close()
    #plt.show()



all_domains.to_csv(Path(f"pixl_cache/target_analysis/_AllTargets_Summary_T{threshold}.csv"))

print("Analysis done!")


