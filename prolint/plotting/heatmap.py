"""Heatmap plotters for contact visualization.

This module provides heatmap-style visualizations for contact matrices,
distance matrices, and database contact patterns.
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import get_colormap, COLORS, apply_prolint_style


class HeatmapPlotter(BasePlotter):
    """Heatmap plotter for contact counts or correlation matrices.

    Works with timeseries results (residue × frame heatmaps) or
    shared_contacts results (residue × residue correlation matrices).

    See Also
    --------
    TimeSeriesAnalysis : Generates timeseries data
    SharedContactsAnalysis : Generates correlation matrices
    """

    name = "heatmap"
    required_analysis = "timeseries"  # Primary, but also works with shared_contacts
    description = "2D heatmap for contact counts or correlation matrices"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        # Check for timeseries format
        if "contact_counts" in result.data and "frames" in result.data:
            return  # Valid timeseries result

        # Check for shared_contacts format
        if "matrix" in result.data and "labels" in result.data:
            return  # Valid shared_contacts result

        raise ValueError(
            f"AnalysisResult not suitable for {cls.name}. "
            f"Expected either timeseries result (with frames, contact_counts) "
            f"or shared_contacts result (with matrix, labels)."
        )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        colorscheme: str = "viridis",
        figsize: Optional[Tuple[float, float]] = None,
        show_row_labels: bool = True,
        show_col_labels: bool = True,
        show_colorbar: bool = True,
        xlabel: str = "",
        ylabel: str = "",
        title: str = "",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        aspect: str = "auto",
        ax: Optional[Axes] = None,
        cbar_label: str = "",
        max_display_rows: int = 40,
        max_display_cols: int = 200,
        origin: str = "upper",
    ) -> Tuple[Figure, Axes]:
        """Create a 2D heatmap visualization.

        Parameters
        ----------
        result : AnalysisResult
            Result from timeseries or shared_contacts analysis.
        colorscheme : str, default="viridis"
            Color scale name.
        figsize : tuple of float, optional
            Figure size (width, height). Auto-calculated if None.
        vmin, vmax : float, optional
            Color scale limits.
        ax : Axes, optional
            Existing axes to plot on.
        max_display_rows, max_display_cols : int
            Maximum rows/columns before sampling.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data based on result type
        if "contact_counts" in result.data:
            # Timeseries format
            contact_counts = result.data["contact_counts"]
            query_residues = result.data.get(
                "query_residues", list(contact_counts.keys())
            )
            frames = result.data["frames"]

            row_labels = [str(r) for r in query_residues]
            col_labels = [str(f) for f in frames]
            data = np.array([contact_counts[r] for r in query_residues])

            if not xlabel:
                xlabel = "Frame"
            if not ylabel:
                ylabel = "Residue"
            if not cbar_label:
                cbar_label = "Contact Count"
        else:
            # Shared contacts format
            data = np.asarray(result.data["matrix"])
            labels = result.data["labels"]
            row_labels = [str(l) for l in labels]
            col_labels = [str(l) for l in labels]

            if not xlabel:
                xlabel = "Residue"
            if not ylabel:
                ylabel = "Residue"
            if not cbar_label:
                cbar_label = "Shared Frames"

        n_rows, n_cols = data.shape

        # Sample columns if too many
        sample_step = 1
        if n_cols > max_display_cols:
            sample_step = n_cols // max_display_cols
            sampled_indices = np.arange(0, n_cols, sample_step)
            data = data[:, sampled_indices]
            col_labels = [col_labels[i] for i in sampled_indices]
            n_cols = len(sampled_indices)

        # Calculate figure size
        if figsize is None:
            row_height = 0.3 if n_rows <= 30 else (0.2 if n_rows <= 60 else 0.12)
            height = max(4, min(n_rows * row_height + 2, 12))
            width = max(8, min(n_cols * 0.05 + 3, 16))
            figsize = (width, height)

        # Create figure/axes
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Get colormap
        cmap = get_colormap(colorscheme)

        # Determine color range
        if vmin is None:
            vmin = np.nanmin(data)
        if vmax is None:
            vmax = np.nanmax(data)

        # Plot heatmap
        im = ax.imshow(
            data,
            aspect=aspect,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            origin=origin,
        )

        # Configure axes
        ax.set_xlabel(xlabel, fontsize=11, fontweight="medium")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="medium")
        if title:
            ax.set_title(title, fontsize=12, fontweight="semibold")

        # Row labels (Y-axis)
        if show_row_labels and n_rows <= max_display_rows:
            ax.set_yticks(np.arange(n_rows))
            ax.set_yticklabels(row_labels, fontsize=9)
        else:
            tick_step = max(1, n_rows // 10)
            ax.set_yticks(np.arange(0, n_rows, tick_step))
            ax.set_yticklabels(
                [row_labels[i] for i in range(0, n_rows, tick_step)], fontsize=9
            )

        # Column labels (X-axis)
        if show_col_labels and n_cols <= 20:
            ax.set_xticks(np.arange(n_cols))
            ax.set_xticklabels(col_labels, fontsize=9, rotation=45, ha="right")
        else:
            tick_step = max(1, n_cols // 10)
            ax.set_xticks(np.arange(0, n_cols, tick_step))
            ax.set_xticklabels(
                [col_labels[i] for i in range(0, n_cols, tick_step)], fontsize=9
            )

        # Colorbar
        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
            if cbar_label:
                cbar.set_label(cbar_label, fontsize=10)
            cbar.ax.tick_params(labelsize=9)

        plt.tight_layout()

        # Add sampling note
        if sample_step > 1:
            ax.annotate(
                f"Showing every {sample_step}th column",
                xy=(0.01, -0.12),
                xycoords="axes fraction",
                ha="center",
                fontsize=8,
                color=COLORS["text"]["secondary"],
            )

        return fig, ax


class DistanceHeatmapPlotter(BasePlotter):
    """Heatmap plotter for atom-atom distance matrices.

    Visualizes pairwise distances between atoms of two residues
    at a specific frame.

    See Also
    --------
    AtomDistancesAnalysis : Generates distance matrix data
    """

    name = "distance_heatmap"
    required_analysis = "atom_distances"
    description = "Atom-atom distance matrix heatmap"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains distance matrix data."""
        required_keys = ["distance_matrix", "query_atoms", "database_atoms"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from 'atom_distances' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        colorscheme: str = "viridis",
        figsize: Optional[Tuple[float, float]] = None,
        show_colorbar: bool = True,
        xlabel: str = "Database Atoms",
        ylabel: str = "Query Atoms",
        title: str = "",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        ax: Optional[Axes] = None,
        cbar_label: str = "Distance (Å)",
        annotate: bool = False,
        annotation_fontsize: int = 8,
    ) -> Tuple[Figure, Axes]:
        """Create atom-atom distance matrix heatmap.

        Parameters
        ----------
        result : AnalysisResult
            Result from atom_distances analysis.
        colorscheme : str, default="viridis"
            Color scale name.
        annotate : bool, default=False
            Whether to annotate cells with distance values.
        ax : Axes, optional
            Existing axes to plot on.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data
        distance_matrix = np.asarray(result.data["distance_matrix"])
        query_atoms = result.data["query_atoms"]
        database_atoms = result.data["database_atoms"]

        n_query, n_db = distance_matrix.shape

        # Auto title
        if not title:
            query_res = result.metadata.get("query_residue", "?")
            db_res = result.metadata.get("database_residue", "?")
            frame = result.data.get("frame", "?")
            title = f"Atom Distances: Res {query_res} ↔ Res {db_res} (Frame {frame})"

        # Calculate figure size
        if figsize is None:
            width = max(6, min(n_db * 0.5 + 2, 14))
            height = max(5, min(n_query * 0.4 + 2, 12))
            figsize = (width, height)

        # Create figure/axes
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Get colormap
        cmap = get_colormap(colorscheme)

        # Determine color range
        if vmin is None:
            vmin = np.nanmin(distance_matrix)
        if vmax is None:
            vmax = np.nanmax(distance_matrix)

        # Plot heatmap
        im = ax.imshow(
            distance_matrix,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )

        # Configure axes
        ax.set_xlabel(xlabel, fontsize=11, fontweight="medium")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="medium")
        ax.set_title(title, fontsize=12, fontweight="semibold")

        # Set tick labels
        ax.set_xticks(np.arange(n_db))
        ax.set_xticklabels(database_atoms, fontsize=9, rotation=45, ha="right")
        ax.set_yticks(np.arange(n_query))
        ax.set_yticklabels(query_atoms, fontsize=9)

        # Annotate cells with values
        if annotate and n_query <= 20 and n_db <= 20:
            for i in range(n_query):
                for j in range(n_db):
                    value = distance_matrix[i, j]
                    # Choose text color based on background
                    text_color = "white" if value < (vmin + vmax) / 2 else "black"
                    ax.text(
                        j,
                        i,
                        f"{value:.1f}",
                        ha="center",
                        va="center",
                        fontsize=annotation_fontsize,
                        color=text_color,
                    )

        # Colorbar
        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label(cbar_label, fontsize=10)
            cbar.ax.tick_params(labelsize=9)

        plt.tight_layout()
        return fig, ax


class DatabaseContactsHeatmapPlotter(BasePlotter):
    """Heatmap plotter for per-molecule contact timelines.

    Shows binary contact states (on/off) for each database molecule
    across trajectory frames for a single query residue.

    See Also
    --------
    DatabaseContactsAnalysis : Generates contact timeline data
    """

    name = "database_contacts_heatmap"
    required_analysis = "database_contacts"
    description = "Per-residue database contact timeline heatmap"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains database contacts data."""
        required_keys = ["contact_matrix", "database_ids", "frames"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from 'database_contacts' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        colorscheme: str = "blues",
        figsize: Optional[Tuple[float, float]] = None,
        show_colorbar: bool = True,
        xlabel: str = "Frame",
        ylabel: str = "Database Molecule ID",
        title: str = "",
        ax: Optional[Axes] = None,
        cbar_label: str = "Contact",
        max_display_rows: int = 50,
        max_display_cols: int = 500,
    ) -> Tuple[Figure, Axes]:
        """Create binary contact timeline heatmap.

        Parameters
        ----------
        result : AnalysisResult
            Result from database_contacts analysis.
        colorscheme : str, default="blues"
            Color scale name.
        max_display_rows, max_display_cols : int
            Maximum rows/columns before sampling.
        ax : Axes, optional
            Existing axes to plot on.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data
        contact_matrix = result.data["contact_matrix"]
        database_ids = result.data["database_ids"]
        frames = result.data["frames"]

        # Build matrix from dict
        sorted_ids = sorted(database_ids)
        data = np.array([contact_matrix[db_id] for db_id in sorted_ids])
        row_labels = [str(db_id) for db_id in sorted_ids]
        col_labels = [str(f) for f in frames]

        n_rows, n_cols = data.shape

        # Auto title
        if not title:
            query_res = result.metadata.get("query_residue", "?")
            db_type = result.metadata.get("database_type", "")
            title = f"Database Contacts: Res {query_res}" + (
                f" - {db_type}" if db_type else ""
            )

        # Sample columns if too many
        sample_step = 1
        if n_cols > max_display_cols:
            sample_step = n_cols // max_display_cols
            sampled_indices = np.arange(0, n_cols, sample_step)
            data = data[:, sampled_indices]
            col_labels = [col_labels[i] for i in sampled_indices]
            n_cols = len(sampled_indices)

        # Calculate figure size
        if figsize is None:
            row_height = 0.25 if n_rows <= 30 else (0.15 if n_rows <= 60 else 0.1)
            height = max(4, min(n_rows * row_height + 2, 12))
            width = max(8, min(n_cols * 0.03 + 3, 16))
            figsize = (width, height)

        # Create figure/axes
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Get colormap
        cmap = get_colormap(colorscheme)

        # Plot heatmap
        im = ax.imshow(
            data,
            aspect="auto",
            cmap=cmap,
            vmin=0,
            vmax=1,
            interpolation="nearest",
            origin="upper",
        )

        # Configure axes
        ax.set_xlabel(xlabel, fontsize=11, fontweight="medium")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="medium")
        ax.set_title(title, fontsize=12, fontweight="semibold")

        # Row labels (Y-axis)
        if n_rows <= max_display_rows:
            ax.set_yticks(np.arange(n_rows))
            ax.set_yticklabels(row_labels, fontsize=8)
        else:
            tick_step = max(1, n_rows // 10)
            ax.set_yticks(np.arange(0, n_rows, tick_step))
            ax.set_yticklabels(
                [row_labels[i] for i in range(0, n_rows, tick_step)], fontsize=8
            )

        # Column labels (X-axis)
        tick_step = max(1, n_cols // 10)
        ax.set_xticks(np.arange(0, n_cols, tick_step))
        ax.set_xticklabels(
            [col_labels[i] for i in range(0, n_cols, tick_step)], fontsize=9
        )

        # Colorbar
        if show_colorbar:
            cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label(cbar_label, fontsize=10)
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(["No", "Yes"])
            cbar.ax.tick_params(labelsize=9)

        plt.tight_layout()

        # Add sampling note
        if sample_step > 1:
            ax.annotate(
                f"Showing every {sample_step}th frame",
                xy=(0.5, -0.12),
                xycoords="axes fraction",
                ha="center",
                fontsize=8,
                color=COLORS["text"]["secondary"],
            )

        return fig, ax


# Register plotters
PlottingRegistry.register("heatmap", HeatmapPlotter)
PlottingRegistry.register("distance_heatmap", DistanceHeatmapPlotter)
PlottingRegistry.register("database_contacts_heatmap", DatabaseContactsHeatmapPlotter)
