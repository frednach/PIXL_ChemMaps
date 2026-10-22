import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

def domain_geometry(group):

        points = group[["x", "y"]].drop_duplicates().to_numpy()

        # If there are enough points to define an area, return NAN
        if len(points) < 3:
            return pd.Series({
                "area_mm2": np.nan,
                "mean_width_mm": np.nan
            })

        try:
            hull = ConvexHull(points)

        # If ConvexHull does not work, return NAN
        except Exception:
            return pd.Series({
                "area_mm2": np.nan,
                "mean_width_mm": np.nan
            })

        # ConvexHull.volume is the area for a 2D hull
        area = hull.volume

        # Coordinates of hull vertices
        hull_points = points[hull.vertices]

        widths = []

        # Examine the width perpendicular to every hull edge
        for i in range(len(hull_points)):

            p1 = hull_points[i]
            p2 = hull_points[
                (i + 1) % len(hull_points)
            ]
            edge = p2 - p1

            # Normal vector to the edge
            normal = np.array([
                -edge[1],
                edge[0]
            ])
            # Normalize
            length = np.linalg.norm(normal)

            if length == 0:
                continue

            normal = normal / length
            # Project all points onto the normal
            projections = points @ normal
            width = (
                projections.max()
                - projections.min()
            )

            widths.append(width)

        if len(widths) == 0:
            max_width = np.nan
        else:
            max_width = max(widths)

        #return values as mm2 and mm
        return pd.Series({
            "area_mm2": area*1e+6,
            "mean_width_mm": max_width*1e+3
        })