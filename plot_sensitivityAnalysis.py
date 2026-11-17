import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

#import data created from T_sensitivityAnalysis.py
T_summary = pd.read_csv(Path("pixl_cache/target_analysis/T_SensitivityAnalysis_Summary.csv"))


# Create Sensitivity Analysis plot for all output parameters
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True)
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flatten()
#
ax0.scatter(
    T_summary["T"],
    T_summary["mean_n_domains"],
    marker="o",
    linestyle="solid",
    s = 2
)
ax0.set_ylabel("n_Domains")
#
ax1.scatter(
    T_summary["T"],
    T_summary["mean_domain_size_mm"],
    marker="o",
    linestyle="solid",
    s = 2
)
ax1.set_ylabel("Domain Size [mm]")
#
ax2.scatter(
    T_summary["T"],
    T_summary["mean_domain_Area_mm2"],
    marker="o",
    linestyle="solid",
    s = 2
)
ax2.set_ylabel("Domain Area [mm2]")
#
ax3.scatter(
    T_summary["T"],
    T_summary["mean_classified_fraction"]*100,
    marker="o",
    linestyle="solid",
    s = 2
)
ax3.set_ylabel("Classified [%]")
#
ax4.scatter(
    T_summary["T"],
    T_summary["mean_olivine_pixel_fraction"]*100,
    marker="o",
    linestyle="solid",
    s = 2
)
ax4.set_ylabel("Ol px [%]")
#
ax5.scatter(
    T_summary["T"],
    T_summary["mean_olivine_domain_fraction"]*100,
    marker="o",
    linestyle="solid",
    s = 2
)
ax5.set_ylabel("Ol Domains [%]")
ax4.set_xlabel("Threshold")
ax5.set_xlabel("Threshold")

plt.subplots_adjust(bottom=0.1, right=1.5, top=0.9)
plt.savefig(
            Path(f"pixl_cache/target_analysis/T_Sensitivity.svg"),
            format="svg",
            bbox_inches = "tight"
)
plt.close()


