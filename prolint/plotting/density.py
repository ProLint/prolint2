"""Density map plotters for spatial distribution visualization.

This module provides 2D density and radial density visualizations
for database molecule distributions around query atoms.
"""

from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import COLORS, get_colormap, apply_prolint_style


class DensityMapPlotter(BasePlotter):
    """Plotter for 2D spatial density maps.

    Visualizes the spatial distribution of database molecules
    around the query center of mass.

    See Also
    --------
    DensityMapAnalysis : Generates density map data
    RadialDensityPlotter : Radial density profiles
    """

    name = "density_map"
    required_analysis = "density_map"
    description = "2D spatial density of database molecules around query"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required density map data."""
        required_keys = ["density", "x_edges", "y_edges"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        colorscheme: str = "viridis",
        xlabel: str = "X (Å)",
        ylabel: str = "Y (Å)",
        title: str = "Density Map",
        figsize: Tuple[float, float] = (8, 8),
        ax: Optional[Axes] = None,
        show_colorbar: bool = True,
        cbar_label: str = "Density",
        log_scale: bool = False,
        aspect: str = "equal",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        show_query_contours: bool = True,
        highlight_query_residues: Optional[List[int]] = None,
        highlight_database_ids: Optional[List[int]] = None,
        universe=None,
        frame: int = 0,
    ) -> Tuple[Figure, Axes]:
        """Create 2D density map visualization.

        Parameters
        ----------
        result : AnalysisResult
            Result from density_map analysis.
        colorscheme : str, default="viridis"
            Color scale name.
        log_scale : bool, default=False
            Whether to use logarithmic color scale.
        show_query_contours : bool, default=True
            Whether to overlay query density as contours.
        highlight_query_residues : list of int, optional
            Query residues to highlight on the map.
        highlight_database_ids : list of int, optional
            Database molecules to highlight.
        universe : Universe, optional
            Required for highlighting residues.
        ax : Axes, optional
            Existing axes to plot on.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data from result
        density = np.asarray(result.data["density"])
        x_edges = result.data["x_edges"]
        y_edges = result.data["y_edges"]

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        cmap = get_colormap(colorscheme)

        x_min, x_max = x_edges[0], x_edges[-1]
        y_min, y_max = y_edges[0], y_edges[-1]

        if vmin is None:
            vmin = np.nanmin(density[density > 0]) if log_scale else np.nanmin(density)
        if vmax is None:
            vmax = np.nanmax(density)

        if log_scale:
            density_plot = np.where(density > 0, density, vmin)
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            density_plot = density
            norm = Normalize(vmin=vmin, vmax=vmax)

        im = ax.imshow(
            density_plot.T,
            origin="lower",
            extent=[x_min, x_max, y_min, y_max],
            cmap=cmap,
            norm=norm,
            aspect=aspect,
            interpolation="bilinear",
        )

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")

        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label(cbar_label, fontsize=10)
            if log_scale:
                cbar.ax.set_ylabel(f"{cbar_label} (log scale)", fontsize=10)

        # Overlay query density as contours using pre-computed centers
        if show_query_contours and "query_density" in result.data:
            query_density = np.asarray(result.data["query_density"])
            x_centers = result.data["x_centers"]
            y_centers = result.data["y_centers"]

            X, Y = np.meshgrid(x_centers, y_centers)

            ax.contour(
                X,
                Y,
                query_density.T,
                levels=5,
                colors=COLORS["data"]["highlight"],
                linewidths=1.5,
                alpha=0.8,
            )

            ax.annotate(
                "Query (contours)",
                xy=(0.02, 0.98),
                xycoords="axes fraction",
                fontsize=9,
                color=COLORS["data"]["highlight"],
                va="top",
            )

        # Highlight specific residues if requested
        if universe is not None and (
            highlight_query_residues or highlight_database_ids
        ):
            universe.trajectory[frame]
            query_com = universe.query.atoms.center_of_mass()[:2]

            legend_handles = []

            if highlight_query_residues:
                for resid in highlight_query_residues:
                    residue_atoms = universe.query.select_atoms(f"resid {resid}")
                    if len(residue_atoms) > 0:
                        pos = residue_atoms.center_of_mass()[:2] - query_com
                        scatter = ax.scatter(
                            pos[0],
                            pos[1],
                            s=150,
                            c=COLORS["data"]["query"],
                            edgecolors="white",
                            linewidths=2,
                            marker="o",
                            zorder=10,
                        )
                        ax.annotate(
                            str(resid),
                            (pos[0], pos[1]),
                            xytext=(5, 5),
                            textcoords="offset points",
                            fontsize=9,
                            fontweight="bold",
                            color=COLORS["data"]["query"],
                        )
                if highlight_query_residues:
                    legend_handles.append(
                        plt.Line2D(
                            [0],
                            [0],
                            marker="o",
                            color="w",
                            markerfacecolor=COLORS["data"]["query"],
                            markersize=10,
                            label="Query residues",
                        )
                    )

            if highlight_database_ids:
                for resid in highlight_database_ids:
                    residue_atoms = universe.database.select_atoms(f"resid {resid}")
                    if len(residue_atoms) > 0:
                        pos = residue_atoms.center_of_mass()[:2] - query_com
                        scatter = ax.scatter(
                            pos[0],
                            pos[1],
                            s=150,
                            c=COLORS["data"]["database"],
                            edgecolors="white",
                            linewidths=2,
                            marker="s",
                            zorder=10,
                        )
                        ax.annotate(
                            str(resid),
                            (pos[0], pos[1]),
                            xytext=(5, 5),
                            textcoords="offset points",
                            fontsize=9,
                            fontweight="bold",
                            color=COLORS["data"]["database"],
                        )
                if highlight_database_ids:
                    legend_handles.append(
                        plt.Line2D(
                            [0],
                            [0],
                            marker="s",
                            color="w",
                            markerfacecolor=COLORS["data"]["database"],
                            markersize=10,
                            label="Database residues",
                        )
                    )

            if legend_handles:
                ax.legend(handles=legend_handles, loc="lower right", fontsize=9)

        plt.tight_layout()
        return fig, ax


class RadialDensityPlotter(BasePlotter):
    """Plotter for radial density profiles.

    Visualizes radially-averaged density as a function of distance
    from the query center.

    See Also
    --------
    RadialDensityAnalysis : Generates radial density data
    DensityMapPlotter : 2D density maps
    """

    name = "radial_density"
    required_analysis = "radial_density"
    description = "Radial density profile from 2D density map"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required radial density data."""
        required_keys = ["r_centers", "radial_density"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        xlabel: str = "Distance from Center (Å)",
        ylabel: str = "Radial Density",
        title: str = "Radial Density Profile",
        figsize: Tuple[float, float] = (8, 4),
        ax: Optional[Axes] = None,
    ) -> Tuple[Figure, Axes]:
        """Create radial density profile plot.

        Parameters
        ----------
        result : AnalysisResult
            Result from radial_density analysis.
        ax : Axes, optional
            Existing axes to plot on.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data from result
        r_centers = result.data["r_centers"]
        radial_density = result.data["radial_density"]
        r_max = result.data.get("r_max", max(r_centers) if r_centers else 1.0)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        ax.plot(r_centers, radial_density, color=COLORS["data"]["query"], linewidth=2)
        ax.fill_between(
            r_centers, radial_density, alpha=0.3, color=COLORS["data"]["query"]
        )

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")

        ax.set_xlim(0, r_max)
        ax.set_ylim(0, None)
        ax.grid(True, alpha=0.3, linestyle="--")

        plt.tight_layout()
        return fig, ax


# Register plotters
PlottingRegistry.register("density_map", DensityMapPlotter)
PlottingRegistry.register("radial_density", RadialDensityPlotter)
