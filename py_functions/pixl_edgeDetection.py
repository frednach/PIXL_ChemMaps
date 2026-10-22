import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.preprocessing import StandardScaler
from pathlib import Path
#import local functions
from py_functions import geometry, normalize, mineralogy



def pixl_edgeDetection (pixl_data, oxide_columns, target_name, target_sol,
                         do_save=False, do_show=True, minimum_size=5, thres=40):

    
    target_data = pixl_data[pixl_data["TARGET"] == target_name]
    target_data = target_data[target_data["SOL"] == target_sol]
    #drop points without chemical composition
    target_data = target_data.dropna(
        subset=oxide_columns + ["x", "y"]
    ).copy()

    #skip small values
    if len(target_data["x"]) <= 800:
        print(f"{target_name}: Skipped due to short length")
        return 0, 0


    #set filepaths and filenames
    analysis_dir = f"pixl_cache/target_analysis/{target_name}_{target_sol}"
    plotname_gradient = f"{target_name}_GradientMap.svg"
    plotname_domains = f"{target_name}_ChemDomains_T{thres}.svg"
    plotname_RGB = f"{target_name}_RGB.svg"
    summary_filename = f"{target_name}_Summary_T{thres}.csv"
    target_filename = f"{target_name}_Target_T{thres}.csv"

    # Create directory for target if needed
    directory_path = Path(analysis_dir)
    try:
        directory_path.mkdir()
    except FileExistsError:
        pass
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


    #scale oxide values
    scaler = StandardScaler()
    X = target_data[
        oxide_columns
    ].to_numpy(dtype=float)
    X_scaled = scaler.fit_transform(X)

    #get positions as xy
    xy = target_data[
        ["x", "y"]
    ].to_numpy(dtype=float)

    #Construct neighbor network with Delaunay
    tri = Delaunay(xy)
    edges = np.vstack([
        tri.simplices[:, [0, 1]],
        tri.simplices[:, [1, 2]],
        tri.simplices[:, [2, 0]]
    ])
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    #remove very long connections
    tree = cKDTree(xy)
    distances, _ = tree.query(
        xy,
        k=2
    )
    nearest_distance = distances[:, 1]
    typical_spacing = np.median(
        nearest_distance
    )
    edge_lengths = np.linalg.norm(
        xy[edges[:, 0]] -
        xy[edges[:, 1]],
        axis=1
    )
    max_edge_length = 3 * typical_spacing
    keep = edge_lengths <= max_edge_length
    edges = edges[keep]

    #Calculate chemical differences
    chemical_difference = np.linalg.norm(
        X_scaled[edges[:, 0]]-
        X_scaled[edges[:, 1]],
        axis=1
    )

    #Plot chemical gradient
    edge_midpoints = (
        xy[edges[:, 0]]
        +
        xy[edges[:, 1]]
    ) / 2
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )
    scatter = ax.scatter(
        edge_midpoints[:, 0],
        edge_midpoints[:, 1],
        c=chemical_difference,
        cmap="inferno",
        s=5
    )
    fig.colorbar(
        scatter,
        ax=ax,
        label="Multivariate chemical difference"
    )
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_aspect("equal")
    ax.set_title(
        f"{target_name} Chemical Gradient"
    )
    if do_save:
        plt.savefig(
            Path(analysis_dir,plotname_gradient),
            format="svg"
        )
    if do_show:
        plt.show()
    plt.close()

    #Extract boundaries above the threshold percentile of gradient
    threshold = np.percentile(
        chemical_difference,
        thres
    )
    boundary_edges = edges[
        chemical_difference >= threshold
    ]
    #convert to line segments
    segments = [
        [
            xy[i],
            xy[j]
        ]
        for i, j in boundary_edges
    ]


    #plot chemical map in RGB
    #plot boundaries on top of the chemical map for the locality
    #create figure
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )
    # Original chemical map
    ax.scatter(
        target_data["x"],
        target_data["y"],
        c=np.column_stack([
            normalize.normalize_percentile(
                target_data["FeO-T_wt%"]
            ),
            normalize.normalize_percentile(
                target_data["MgO_wt%"]
            ),
            normalize.normalize_percentile(
                target_data["SiO2_wt%"]
            )
        ]),
        marker="s",
        s=5
    )
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_aspect("equal")
    ax.set_title(
        f"{target_name} Chemical Map R,G,B = FeO-T, MgO, SiO2"
    )
    if do_save:
        plt.savefig(
            Path(analysis_dir,plotname_RGB),
            format="svg"
        )
    if do_show:
        plt.show()
    plt.close()



    #plot boundaries on top of the chemical map for the locality
    #create figure
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )
    # Original chemical map
    ax.scatter(
        target_data["x"],
        target_data["y"],
        c=np.column_stack([
            normalize.normalize_percentile(
                target_data["FeO-T_wt%"]
            ),
            normalize.normalize_percentile(
                target_data["MgO_wt%"]
            ),
            normalize.normalize_percentile(
                target_data["SiO2_wt%"]
            )
        ]),
        marker="s",
        s=5
    )
    # Candidate boundaries
    lc = LineCollection(
        segments,
        linewidths=1.0
    )
    ax.add_collection(lc)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_aspect("equal")
    ax.set_title(
        f"{target_name} Candidate Chemical Boundaries T={thres}"
    )
    plotname_Boundaries = f"{target_name}_BoundaryMap_T{thres}.svg"
    if do_save:
        plt.savefig(
            Path(analysis_dir,plotname_Boundaries),
            format="svg"
        )
    if do_show:
        plt.show()
    plt.close()

    #Construct the within-domain network
    #
    domain_edges = edges[
        chemical_difference < threshold
    ]

    #convert network to sparse graph
    #create undirected graph to identify neighboring domains of similarity
    n_points = len(target_data)
    row = np.concatenate([
        domain_edges[:, 0],
        domain_edges[:, 1]
    ])
    col = np.concatenate([
        domain_edges[:, 1],
        domain_edges[:, 0]
    ])
    data = np.ones(len(row))
    graph = coo_matrix(
        (data, (row, col)),
        shape=(n_points, n_points)
    )

    #Find connected components
    n_domains, domain_labels = connected_components(
        graph,
        directed=False
    )
    #add to domains in DataFrame
    target_data["domain"] = domain_labels

    #get domain sizes
    domain_sizes = (
        target_data
        .groupby("domain")
        .size()
        .sort_values(ascending=False)
    )
    #force a minimum size of chemical domains
    valid_domains = domain_sizes[
        domain_sizes >= minimum_size
    ].index

    #create new classification in DataFrame for valid domains (-1 = too small or unclassified)
    target_data["domain_filtered"] = np.where(
        target_data["domain"].isin(valid_domains),
        target_data["domain"],
        -1
    )



    

    # Number of actual domains
    valid_domain_ids = sorted(
        target_data.loc[
            target_data["domain_filtered"] >= 0,
            "domain_filtered"
        ].unique()
    )
    n_valid = len(valid_domain_ids)

    # Generate enough colors
    base_colors = (
        plt.cm.gist_ncar(
            np.linspace(0, 1, n_valid)
        )
    )
    # Map domain IDs to consecutive color indices
    domain_to_color = {
        domain_id: i
        for i, domain_id in enumerate(valid_domain_ids)
    }
    plot_colors = target_data["domain_filtered"].map(
        domain_to_color
    )
    # Give unclassified points a separate value
    plot_colors = plot_colors.fillna(-1)


    #create plot
    fig, ax = plt.subplots(
        figsize=(10, 6)
    )
    # Plot unclassified points first as gray
    unclassified = (
        target_data["domain_filtered"] == -1
    )
    ax.scatter(
        target_data.loc[unclassified, "x"],
        target_data.loc[unclassified, "y"],
        color="lightgray",
        marker="s",
        s=8,
        label="Unclassified"
    )
    # Plot actual domains
    classified = ~unclassified
    scatter = ax.scatter(
        target_data.loc[classified, "x"],
        target_data.loc[classified, "y"],
        c=plot_colors[classified],
        cmap=mcolors.ListedColormap(base_colors),
        marker="s",
        s=8
    )

    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_aspect("equal")

    ax.set_title(
        f"{target_name} Chemically Defined Domains T={thres}"
    )
    if do_save:
        plt.savefig(
                    Path(analysis_dir,plotname_domains),
                    format="svg"
                )
    if do_show:
        plt.show()
    plt.close()






    #Calculate each domains chemistry
    domain_chemistry = (
        target_data[
            target_data["domain_filtered"] != -1
        ]
        .groupby("domain_filtered")[oxide_columns]
        .mean()
    )

    #Calculate chemical domain size
    domain_summary = (
        target_data[
            target_data["domain_filtered"] != -1
        ]
        .groupby("domain_filtered")
        .agg(
            n_points=("domain_filtered", "size"),
            x_min=("x", "min"),
            x_max=("x", "max"),
            y_min=("y", "min"),
            y_max=("y", "max")
        )
    )
    #Add domain chemistry
    domain_summary = domain_summary.join(
        domain_chemistry
    )

    #Add domain geometry
    domain_geometry_df = (
    target_data[
        target_data["domain_filtered"] != -1
    ]
    .groupby("domain_filtered")
    .apply(geometry.domain_geometry)
    )
    #Join on summary csv file
    domain_summary = domain_summary.join(
    domain_geometry_df
    )

    #Create columns with standard deviation for every oxide
    oxide_std = (
    target_data[
        target_data["domain_filtered"] != -1
    ]
    .groupby("domain_filtered")[oxide_columns]
    .std()
    )
    #Give names to the std columns as oxide_std
    oxide_std = oxide_std.rename(
        columns={
            oxide: f"{oxide}_std"
            for oxide in oxide_columns
        }
    )
    #Join on existing summary file
    domain_summary = domain_summary.join(oxide_std)

    #Add olivine classification to domain summary
    domain_summary = mineralogy.calculate_fo_fa(domain_summary)


    if do_save:
        #Save data to .csv files
        target_data.to_csv(
            Path(analysis_dir,target_filename),
            index=False
        )
        domain_summary.to_csv(
            Path(analysis_dir,summary_filename)
        )

    return domain_summary, target_data
