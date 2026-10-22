#
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import mpltern
from adjustText import adjust_text


def calculate_fo_fa(
    summary_df,
    mg_column="MgO_wt%",
    fe_column="FeO-T_wt%",
    si_column="SiO2_wt%",
    al_column="Al2O3_wt%",
    ca_column="CaO_wt%",
    ratio_low=1.8,
    ratio_high=2.2,
    max_al=5.0,
    max_ca=5.0
):

    # Molar masses [g/mol]
    MW_MgO = 40.304
    MW_FeO = 71.844
    MW_SiO2 = 60.084

    # Convert oxide concentrations to cation amounts
    Mg = summary_df[mg_column] / MW_MgO
    Fe = summary_df[fe_column] / MW_FeO
    Si = summary_df[si_column] / MW_SiO2

    MgFe = Mg + Fe

    valid = (
        (MgFe > 0) &
        (Si > 0)
    )

    # Initialize
    summary_df["Fo"] = np.nan
    summary_df["Fa"] = np.nan
    summary_df["MgFe_Si_ratio"] = np.nan
    summary_df["is_olivine_candidate"] = False

    # Fo
    summary_df.loc[valid, "Fo"] = (
        100 * Mg[valid] / MgFe[valid]
    )

    # Fa
    summary_df.loc[valid, "Fa"] = (
        100 * Fe[valid] / MgFe[valid]
    )

    # Stoichiometric ratio
    summary_df.loc[valid, "MgFe_Si_ratio"] = (
        MgFe[valid] / Si[valid]
    )

    # Individual criteria
    ratio_ok = (
        summary_df["MgFe_Si_ratio"]
        .between(ratio_low, ratio_high)
    )

    al_ok = (
        summary_df[al_column] <= max_al
    )

    ca_ok = (
        summary_df[ca_column] <= max_ca
    )

    # Combined olivine candidate criterion
    summary_df["is_olivine_candidate"] = (
        valid &
        ratio_ok &
        al_ok &
        ca_ok
    )

    return summary_df








def plot_pixl_ternary(domain_summary, target_name=None):

    # --------------------------------------------------
    # Molar masses [g/mol]
    # --------------------------------------------------

    MW_Al2O3 = 101.961
    MW_CaO = 56.077
    MW_Na2O = 61.979
    MW_K2O = 94.196
    MW_FeO = 71.844
    MW_MgO = 40.304

    # --------------------------------------------------
    # Convert oxide wt% to molar abundances
    # --------------------------------------------------

    Al = domain_summary["Al2O3_wt%"] / MW_Al2O3

    Ca = domain_summary["CaO_wt%"] / MW_CaO
    Na = domain_summary["Na2O_wt%"] / MW_Na2O
    K = domain_summary["K2O_wt%"] / MW_K2O

    Fe = domain_summary["FeO-T_wt%"] / MW_FeO
    Mg = domain_summary["MgO_wt%"] / MW_MgO

    # --------------------------------------------------
    # Combine components to plotting axis t,l,r
    # --------------------------------------------------

    Al_component = Al

    CNK_component = (
        Ca + Na + K
    )

    FM_component = (
        Fe + Mg
    )

    # --------------------------------------------------
    # Normalize to 100%
    # --------------------------------------------------

    total = (
        Al_component
        + CNK_component
        + FM_component
    )

    Al_norm = Al_component / total
    CNK_norm = CNK_component / total
    FM_norm = FM_component / total

    # --------------------------------------------------
    # Create ternary figure
    # --------------------------------------------------

    fig = plt.figure(
        figsize=(8, 8)
    )

    tax = fig.add_subplot(
        111,
        projection="ternary"
    )

    # --------------------------------------------------
    # Plot all domains
    #
    # mpltern:
    #     top, left, right
    #
    # Therefore:
    #     Al2O3       -> top
    #     Ca+Na+K     -> left
    #     Fe+Mg       -> right
    # --------------------------------------------------

    tax.scatter(
        Al_norm,
        CNK_norm,
        FM_norm,
        s=40,
        alpha=0.6,
        label="Chemical domains"
    )

    # --------------------------------------------------
    # Highlight probable olivine domains
    # --------------------------------------------------

    if "is_olivine_candidate" in domain_summary.columns:

        olivine = domain_summary[
            domain_summary["is_olivine_candidate"]
        ].copy()

        if len(olivine) > 0:

            Al_ol = (
                olivine["Al2O3_wt%"]
                / MW_Al2O3
            )

            CNK_ol = (
                olivine["CaO_wt%"]
                / MW_CaO
                +
                olivine["Na2O_wt%"]
                / MW_Na2O
                +
                olivine["K2O_wt%"]
                / MW_K2O
            )

            FM_ol = (
                olivine["FeO-T_wt%"]
                / MW_FeO
                +
                olivine["MgO_wt%"]
                / MW_MgO
            )

            total_ol = (
                Al_ol
                + CNK_ol
                + FM_ol
            )

            Al_ol /= total_ol
            CNK_ol /= total_ol
            FM_ol /= total_ol

            tax.scatter(
                Al_ol,
                CNK_ol,
                FM_ol,
                s=70,
                marker="o",
                edgecolor="black",
                linewidth=0.8,
                label="Probable olivine"
            )

    # --------------------------------------------------
    # Axis labels
    # --------------------------------------------------
    tax.set_tlabel(
        "Al$_2$O$_3$ [mol]"
    )

    tax.set_llabel(
        "CaO + Na$_2$O + K$_2$O [mol]"
    )

    tax.set_rlabel(
        "FeO-T + MgO [mol]"
    )

    # --------------------------------------------------
    # Grid
    # --------------------------------------------------
    tax.grid(
        linewidth=0.8,
        alpha=0.5
    )

    # --------------------------------------------------
    # Title
    # --------------------------------------------------
    if target_name == None:
        tax.set_title(
                f"PIXL mineral chemistry"
        )
    else:
        tax.set_title(
            f"{target_name}: PIXL mineral chemistry"
        )

    # --------------------------------------------------
    # Legend
    # --------------------------------------------------
    tax.legend(
        loc="upper right"
    )
    plt.tight_layout()

    return fig, tax









def plot_both_ternary(domain_summary, threshold, target_name=None):

    # ==================================================
    # MOLAR MASSES [g/mol]
    # ==================================================

    MW_Al2O3 = 101.961
    MW_CaO = 56.077
    MW_Na2O = 61.979
    MW_K2O = 94.196
    MW_FeO = 71.844
    MW_MgO = 40.304
    MW_SiO2 = 60.084

    # ==================================================
    # CREATE FIGURE
    # ==================================================

    fig = plt.figure(figsize=(15, 7))

    # Two ternary axes
    tax1 = fig.add_subplot(
        1, 2, 1,
        projection="ternary"
    )

    tax2 = fig.add_subplot(
        1, 2, 2,
        projection="ternary"
    )

    # ==================================================
    # GENERAL PIXL MINERAL CHEMISTRY
    # ==================================================

    Al = (
        domain_summary["Al2O3_wt%"]
        / MW_Al2O3
    )

    Ca = (
        domain_summary["CaO_wt%"]
        / MW_CaO
    )

    Na = (
        domain_summary["Na2O_wt%"]
        / MW_Na2O
    )

    K = (
        domain_summary["K2O_wt%"]
        / MW_K2O
    )

    Fe = (
        domain_summary["FeO-T_wt%"]
        / MW_FeO
    )

    Mg = (
        domain_summary["MgO_wt%"]
        / MW_MgO
    )

    # Combined components
    CNK = Ca + Na + K
    FM = Fe + Mg

    total = Al + CNK + FM

    Al_norm = Al / total
    CNK_norm = CNK / total
    FM_norm = FM / total

    # --------------------------------------------------
    # Plot all domains
    #
    # mpltern order:
    # top, left, right
    #
    # top  = Al2O3
    # left = CaO + Na2O + K2O
    # right = FeO-T + MgO
    # --------------------------------------------------

    tax1.scatter(
        Al_norm,
        CNK_norm,
        FM_norm,
        s=40,
        alpha=0.5,
        label="All domains"
    )

    # --------------------------------------------------
    # Highlight probable olivine
    # --------------------------------------------------

    if "is_olivine_candidate" in domain_summary.columns:

        olivine = domain_summary[
            domain_summary["is_olivine_candidate"]
        ].copy()

        if len(olivine) > 0:

            Al_ol = (
                olivine["Al2O3_wt%"]
                / MW_Al2O3
            )

            CNK_ol = (
                olivine["CaO_wt%"]
                / MW_CaO
                +
                olivine["Na2O_wt%"]
                / MW_Na2O
                +
                olivine["K2O_wt%"]
                / MW_K2O
            )

            FM_ol = (
                olivine["FeO-T_wt%"]
                / MW_FeO
                +
                olivine["MgO_wt%"]
                / MW_MgO
            )

            total_ol = (
                Al_ol
                + CNK_ol
                + FM_ol
            )

            Al_ol /= total_ol
            CNK_ol /= total_ol
            FM_ol /= total_ol

            tax1.scatter(
                Al_ol,
                CNK_ol,
                FM_ol,
                s=70,
                marker="o",
                edgecolor="black",
                linewidth=0.8,
                label="Probable olivine"
            )

        #Add mineral labels as schematic reference
        mineral_labels = {
        "Kaolinite": (0.88, 0.07, 0.05),
        "Muscovite": (0.72, 0.23, 0.05),
        "Clay minerals": (0.58, 0.25, 0.17),
        "Feldspar": (0.48, 0.48, 0.04),
        "Ca-sulfate /\nCarbonate":
            (0.04, 0.91, 0.05),
        "Chlorite": (0.20, 0.05, 0.75),
        "Cpx": (0.07, 0.06, 0.87),
        "Opx": (0.03, 0.02, 0.95),
    }

    texts = []

    for label, coords in mineral_labels.items():

        tax1.scatter(
            *coords,
            s=60,
            marker="D"
        )

        text = tax1.text(
            *coords,
            label,
            fontsize=14
        )

        texts.append(text)

    
    # Axis labels
    tax1.set_tlabel(
        "Al$_2$O$_3$ [mol%]"
    )

    tax1.set_llabel(
        "CaO + Na$_2$O + K$_2$O [mol%]"
    )

    tax1.set_rlabel(
        "FeO-T + MgO [mol%]"
    )

    tax1.grid(
        linewidth=0.8,
        alpha=0.5
    )

    tax1.set_title(
        "General PIXL mineral chemistry"
    )

    tax1.legend(
        loc="upper right"
    )

    #Adjust Text in Figure to not overlap
    plt.subplots_adjust(
    left=0.04,
    right=0.96,
    top=0.88,
    bottom=0.08,
    wspace=0.40
    )
    adjust_text(
        texts,
        expand=(1.5, 2),
        ax=tax1,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            linewidth=2,
            alpha=1
        )
    )





    # ==================================================
    # OLIVINE-FOCUSED MgO-FeO-T-SiO2 TERNARY
    # ==================================================

    Si = (
        domain_summary["SiO2_wt%"]
        / MW_SiO2
    )

    total_olivine_ternary = (
        Mg + Fe + Si
    )

    Mg_norm = (
        Mg / total_olivine_ternary
    )

    Fe_norm = (
        Fe / total_olivine_ternary
    )

    Si_norm = (
        Si / total_olivine_ternary
    )

    # --------------------------------------------------
    # Plot all domains
    #
    # top  = SiO2
    # left = MgO
    # right = FeO-T
    # --------------------------------------------------

    tax2.scatter(
        Si_norm,
        Mg_norm,
        Fe_norm,
        s=40,
        alpha=0.5,
        label="All domains"
    )

    # --------------------------------------------------
    # Highlight probable olivine
    # --------------------------------------------------

    if "is_olivine_candidate" in domain_summary.columns:

        olivine = domain_summary[
            domain_summary["is_olivine_candidate"]
        ].copy()

        if len(olivine) > 0:

            Mg_ol = (
                olivine["MgO_wt%"]
                / MW_MgO
            )

            Fe_ol = (
                olivine["FeO-T_wt%"]
                / MW_FeO
            )

            Si_ol = (
                olivine["SiO2_wt%"]
                / MW_SiO2
            )

            total = (
                Mg_ol
                + Fe_ol
                + Si_ol
            )

            Mg_ol /= total
            Fe_ol /= total
            Si_ol /= total

            tax2.scatter(
                Si_ol,
                Mg_ol,
                Fe_ol,
                s=70,
                marker="o",
                edgecolor="black",
                linewidth=0.8,
                label="Probable olivine"
            )

    # --------------------------------------------------
    # Ideal olivine series
    # --------------------------------------------------

    Fo = np.linspace(0, 1, 200)

    Mg_ideal = 2 * Fo
    Fe_ideal = 2 * (1 - Fo)
    Si_ideal = np.ones_like(Fo)

    total = (
        Mg_ideal
        + Fe_ideal
        + Si_ideal
    )

    Mg_ideal /= total
    Fe_ideal /= total
    Si_ideal /= total

    tax2.plot(
        Si_ideal,
        Mg_ideal,
        Fe_ideal,
        linestyle="--",
        linewidth=2,
        label="Ideal olivine"
    )

    # Axis labels
    tax2.set_tlabel(
        "SiO$_2$ [mol%]"
    )

    tax2.set_llabel(
        "MgO [mol%]"
    )

    tax2.set_rlabel(
        "FeO-T [mol%]"
    )

    tax2.grid(
        linewidth=0.8,
        alpha=0.5
    )

    tax2.set_title(
        "Olivine-focused chemistry"
    )

    tax2.legend(
        loc="upper right"
    )

    # ==================================================
    # OVERALL TITLE
    # ==================================================

    if target_name is not None:

        fig.suptitle(
            f"{target_name} PIXL chemical domains T={threshold}",
            fontsize=16
        )

    else:

        fig.suptitle(
            "PIXL chemical domains",
            fontsize=16
        )

    # ==================================================
    # LAYOUT
    # ==================================================

    plt.subplots_adjust(
        top=0.88,
        wspace=0.35
    )

    return fig, tax1, tax2