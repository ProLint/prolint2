"""Residue metrics plotters for per-residue visualization.

This module provides bar charts and logo grids for displaying
per-residue contact metrics.
"""

from typing import List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import (
    COLORS,
    AMINO_ACID_COLORS,
    AMINO_ACID_ONE_LETTER,
    get_colormap,
    get_color_for_value,
    apply_prolint_style,
)


class ResidueMetricsPlotter(BasePlotter):
    """Plotter for per-residue contact metrics.

    Visualizes metrics as bar charts or scatter plots with
    amino acid coloring and highlighting options.

    See Also
    --------
    MetricsAnalysis : Generates per-residue metric data
    LogoGridPlotter : Grid-based residue visualization
    """

    name = "residue_metrics"
    required_analysis = "metrics"
    description = "Per-residue metrics visualization (bar/scatter)"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains residue metrics data."""
        required_keys = ["residues", "values"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis with "
                f"keys: residues (list of dicts with resid, resname), values (list of floats)."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        style: str = "bar",
        colorscheme: str = "prolint",
        xlabel: str = "Residue",
        ylabel: str = "Value",
        title: str = "Per-Residue Metrics",
        figsize: Optional[Tuple[float, float]] = None,
        ax: Optional[Axes] = None,
        show_aa_labels: bool = True,
        highlight_residues: Optional[List[int]] = None,
        bar_width: float = 0.8,
        sort_by_value: bool = False,
        marker_size: int = 50,
    ) -> Tuple[Figure, Axes]:
        """Create per-residue metrics visualization.

        Parameters
        ----------
        result : AnalysisResult
            Result from metrics analysis.
        style : {"bar", "scatter"}, default="bar"
            Plot style.
        colorscheme : str, default="prolint"
            Color scheme ("prolint", "amino_acid", or scale name).
        show_aa_labels : bool, default=True
            Whether to show amino acid labels on x-axis.
        highlight_residues : list of int, optional
            Residue IDs to highlight.
        sort_by_value : bool, default=False
            Whether to sort residues by value.
        ax : Axes, optional
            Existing axes to plot on.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        residues = result.data["residues"]
        values = result.data["values"]

        # Pair residues with values for sorting
        paired = list(zip(residues, values))
        if sort_by_value:
            paired = sorted(paired, key=lambda x: x[1], reverse=True)
            residues, values = zip(*paired) if paired else ([], [])
            residues = list(residues)
            values = list(values)

        n_residues = len(residues)
        resids = [r["resid"] for r in residues]
        resnames = [r.get("resname", "X") for r in residues]

        # Calculate figure size
        if figsize is None:
            width = max(8, min(n_residues * 0.3, 20))
            figsize = (width, 4)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Determine colors
        min_val = min(values) if values else 0
        max_val = max(values) if values else 1

        if colorscheme == "amino_acid":
            colors = [
                AMINO_ACID_COLORS.get(rn, COLORS["neutral"]["400"]) for rn in resnames
            ]
        else:
            colors = [
                get_color_for_value(v, colorscheme, min_val, max_val) for v in values
            ]

        highlight_set = set(highlight_residues or [])

        if style == "scatter":
            edgecolors = [
                COLORS["data"]["highlight"] if r in highlight_set else "white"
                for r in resids
            ]
            linewidths = [2 if r in highlight_set else 0.5 for r in resids]

            ax.scatter(
                resids,
                values,
                c=colors,
                s=marker_size,
                edgecolors=edgecolors,
                linewidths=linewidths,
            )
            ax.set_xlim(min(resids) - 1, max(resids) + 1)
        else:
            # Bar chart style (default)
            edgecolors = []
            linewidths = []
            for resid in resids:
                if resid in highlight_set:
                    edgecolors.append(COLORS["data"]["highlight"])
                    linewidths.append(2)
                else:
                    edgecolors.append("white")
                    linewidths.append(0.5)

            x = np.arange(n_residues)
            ax.bar(
                x,
                values,
                width=bar_width,
                color=colors,
                edgecolor=edgecolors,
                linewidth=linewidths,
            )

            # X-axis ticks for bar chart
            if n_residues <= 50:
                ax.set_xticks(x)
                if show_aa_labels:
                    labels = [
                        f"{AMINO_ACID_ONE_LETTER.get(rn, rn[0])}{resids[i]}"
                        for i, rn in enumerate(resnames)
                    ]
                    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
                else:
                    ax.set_xticklabels(
                        [str(r) for r in resids], fontsize=8, rotation=45, ha="right"
                    )
            else:
                step = max(1, n_residues // 20)
                ax.set_xticks(x[::step])
                if show_aa_labels:
                    labels = [
                        f"{AMINO_ACID_ONE_LETTER.get(resnames[i], resnames[i][0])}{resids[i]}"
                        for i in range(0, n_residues, step)
                    ]
                else:
                    labels = [str(resids[i]) for i in range(0, n_residues, step)]
                ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")

            ax.set_xlim(-0.5, n_residues - 0.5)

        # Common axis settings
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")
        ax.set_ylim(0, max_val * 1.1 if max_val > 0 else 1)
        ax.grid(True, alpha=0.3, linestyle="--", axis="y")

        # Add colorbar for value-based coloring
        if colorscheme != "amino_acid":
            from matplotlib.cm import ScalarMappable

            cmap = get_colormap(colorscheme)
            norm = Normalize(vmin=min_val, vmax=max_val)
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label(ylabel, fontsize=9)

        plt.tight_layout()
        return fig, ax


class LogoGridPlotter(BasePlotter):
    """Plotter for grid-based residue logo visualization.

    Displays residues as colored cells arranged in rows with
    one-letter amino acid codes and residue numbers.

    See Also
    --------
    MetricsAnalysis : Generates per-residue metric data
    ResidueMetricsPlotter : Bar/scatter visualization
    """

    name = "logo_grid"
    required_analysis = "metrics"
    description = "Grid-based residue logo plot with amino acid annotations"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains residue metrics data."""
        required_keys = ["residues", "values"]
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
        colorscheme: str = "prolint",
        residues_per_row: int = 80,
        cell_size: float = 0.3,
        title: str = "Residue Logo Plot",
        figsize: Optional[Tuple[float, float]] = None,
        highlight_residues: Optional[List[int]] = None,
    ) -> Tuple[Figure, Axes]:
        """Create grid-based residue logo visualization.

        Parameters
        ----------
        result : AnalysisResult
            Result from metrics analysis.
        colorscheme : str, default="prolint"
            Color scale name for value-based coloring.
        residues_per_row : int, default=80
            Number of residue cells per row.
        cell_size : float, default=0.3
            Size of each cell in inches.
        title : str, default="Residue Logo Plot"
            Plot title.
        figsize : tuple of (float, float), optional
            Figure dimensions (width, height). Auto-calculated if None.
        highlight_residues : list of int, optional
            Residue IDs to highlight with colored borders.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        residues = result.data["residues"]
        values = result.data["values"]

        n_residues = len(residues)
        n_rows = (n_residues + residues_per_row - 1) // residues_per_row

        # Calculate figure size
        if figsize is None:
            width = min(residues_per_row, n_residues) * cell_size + 1
            height = n_rows * cell_size * 1.5 + 1
            figsize = (width, height)

        fig, ax = plt.subplots(figsize=figsize)

        min_val = min(values) if values else 0
        max_val = max(values) if values else 1
        highlight_set = set(highlight_residues or [])

        for i, (r, value) in enumerate(zip(residues, values)):
            row = i // residues_per_row
            col = i % residues_per_row

            x = col
            y = n_rows - 1 - row  # Flip y so row 0 is at top

            resname = r.get("resname", "X")
            resid = r["resid"]

            # Get color
            color = get_color_for_value(value, colorscheme, min_val, max_val)

            # Determine edge color for highlighting
            if resid in highlight_set:
                edgecolor = COLORS["data"]["highlight"]
                linewidth = 2
            else:
                edgecolor = COLORS["border"]["default"]
                linewidth = 0.5

            # Draw rectangle
            rect = Rectangle(
                (x - 0.4, y - 0.4),
                0.8,
                0.8,
                facecolor=color,
                edgecolor=edgecolor,
                linewidth=linewidth,
            )
            ax.add_patch(rect)

            # Determine text color based on background brightness
            normalized = (
                (value - min_val) / (max_val - min_val) if max_val > min_val else 0.5
            )
            text_color = "black" if normalized > 0.5 else "white"

            # One-letter code
            one_letter = AMINO_ACID_ONE_LETTER.get(
                resname, resname[0] if resname else "?"
            )
            ax.text(
                x,
                y + 0.1,
                one_letter,
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
                color=text_color,
                family="monospace",
            )

            # Residue ID (smaller)
            ax.text(
                x,
                y - 0.2,
                str(resid),
                ha="center",
                va="center",
                fontsize=5,
                color=text_color,
                family="monospace",
            )

        # Configure axes
        ax.set_xlim(-0.5, min(residues_per_row, n_residues) - 0.5)
        ax.set_ylim(-0.5, n_rows - 0.5)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, fontsize=12, fontweight="semibold", pad=10)

        # Add colorbar
        from matplotlib.cm import ScalarMappable

        cmap = get_colormap(colorscheme)
        norm = Normalize(vmin=min_val, vmax=max_val)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", shrink=0.5, pad=0.05)
        cbar.ax.tick_params(labelsize=8)

        plt.tight_layout()
        return fig, ax


# Register plotters
PlottingRegistry.register("residue_metrics", ResidueMetricsPlotter)
PlottingRegistry.register("logo_grid", LogoGridPlotter)
