import pandas as pd
import numpy as np
import py_functions.pixl_edgeDetection as ED
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

# Thresholds to test
thresholds = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]

sensitivity_results = []

# ---------------------------------------------------------
# Loop over every TARGET + SOL observation
# ---------------------------------------------------------

for (target_name, target_sol), observation in pixl_data.groupby(
    ["TARGET", "SOL"]
):

    # Same basic filtering used by the edge-detection function
    observation = observation.dropna(
        subset=oxide_columns + ["x", "y"]
    )

    # Skip small observations
    if len(observation) <= 800:
        print(
            f"{target_name}, SOL {target_sol}: skipped "
            f"({len(observation)} points)"
        )
        continue

    print(
        f"\n{target_name}, SOL {target_sol}"
    )

    # -----------------------------------------------------
    # Run the edge detection for every threshold
    # -----------------------------------------------------

    for T in thresholds:

        domain_summary, target_result = ED.pixl_edgeDetection(
            pixl_data,
            oxide_columns,
            target_name,
            target_sol,
            do_save=False,
            do_show=False,
            minimum_size=5,
            thres=T
        )

        # The function returns 0, 0 for skipped observations
        if isinstance(domain_summary, int):
            continue

        # -------------------------------------------------
        # Number and size of domains
        # -------------------------------------------------

        n_domains = len(domain_summary)

        domain_size = domain_summary["mean_width_mm"].mean()

        domain_area = domain_summary["area_mm2"].mean()

        # -------------------------------------------------
        # Classified pixels
        # -------------------------------------------------

        classified = target_result[
            target_result["domain_filtered"] != -1
        ]

        n_classified = len(classified)

        classified_fraction = (
            n_classified / len(target_result)
        )

        # -------------------------------------------------
        # Probable olivine
        # -------------------------------------------------

        olivine_domains = domain_summary[
            domain_summary["is_olivine_candidate"]
        ]

        n_olivine_domains = len(
            olivine_domains
        )

        if n_domains > 0:
            olivine_domain_fraction = (
                n_olivine_domains / n_domains
            )
        else:
            olivine_domain_fraction = np.nan

        # Domain IDs of probable olivine
        olivine_ids = olivine_domains.index

        # Pixels belonging to probable-olivine domains
        olivine_pixels = classified[
            classified["domain_filtered"].isin(
                olivine_ids
            )
        ]

        n_olivine_pixels = len(
            olivine_pixels
        )

        olivine_pixel_fraction = (
            n_olivine_pixels / len(target_result)
        )

        # -------------------------------------------------
        # Store one result row
        # -------------------------------------------------

        sensitivity_results.append({
            "TARGET": target_name,
            "SOL": target_sol,
            "T": T,
            "n_points": len(target_result),
            "area_mm2": domain_area,
            "domainSize_mm": domain_size,
            "n_domains": n_domains,
            "n_classified": n_classified,
            "classified_fraction": classified_fraction,
            "n_olivine_domains": n_olivine_domains,
            "olivine_domain_fraction": olivine_domain_fraction,
            "n_olivine_pixels": n_olivine_pixels,
            "olivine_pixel_fraction": olivine_pixel_fraction
        })


# ---------------------------------------------------------
# Convert to DataFrame
# ---------------------------------------------------------

sensitivity_df = pd.DataFrame(sensitivity_results)
sensitivity_df.to_csv(Path("pixl_cache/target_analysis/T_SensitivityAnalysis.csv"))



## Calculate mean values for each T value
mean_domains = (
    sensitivity_df
    .groupby("T")["n_domains"]
    .mean()
)
mean_classified = (
    sensitivity_df
    .groupby("T")["classified_fraction"]
    .mean()
)
mean_domainsize = (
    sensitivity_df
    .groupby("T")["domainSize_mm"]
    .mean()
)
mean_domainArea = (
    sensitivity_df
    .groupby("T")["area_mm2"]
    .mean()
)
mean_olivine_px = (
    sensitivity_df
    .groupby("T")["olivine_pixel_fraction"]
    .mean()
)
mean_olivine_domain = (
    sensitivity_df
    .groupby("T")["olivine_domain_fraction"]
    .mean()
)

#combine calculated values into one DataFrame
sensitivity_summary = pd.DataFrame({
    "mean_n_domains": mean_domains,
    "mean_domain_size_mm": mean_domainsize,
    "mean_domain_Area_mm2": mean_domainArea,
    "mean_classified_fraction": mean_classified,
    "mean_olivine_pixel_fraction": mean_olivine_px,
    "mean_olivine_domain_fraction": mean_olivine_domain
})
#Write DataFrame to csv file
sensitivity_summary.to_csv(Path("pixl_cache/target_analysis/T_SensitivityAnalysis_Summary.csv"))