"""Time series plotters for temporal contact visualization.

This module provides time series plots for contact dynamics
and distance evolution over trajectory frames.
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import (
    COLORS,
    COLOR_SCALES,
    apply_prolint_style,
    get_unit_label,
)


class TimeSeriesPlotter(BasePlotter):
    """Plotter for contact time series.

    Displays contact counts over trajectory frames for multiple
    residues as overlaid line plots.

    See Also
    --------
    TimeSeriesAnalysis : Generates time series data
    DistanceTimeSeriesPlotter : Distance evolution plots
    """

    name = "timeseries"
    required_analysis = "timeseries"
    description = "Contact counts over time for multiple residues"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required time series data."""
        required_keys = ["frames", "contact_counts"]
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
        xlabel: str = "Frame",
        ylabel: str = "Contact Count",
        title: str = "Contacts Over Time",
        figsize: Tuple[float, float] = (10, 4),
        ax: Optional[Axes] = None,
        time_units: Optional[str] = None,
        dt: float = 1.0,
        legend: bool = True,
        max_series: int = 10,
    ) -> Tuple[Figure, Axes]:
        """Create contact count time series plot.

        Parameters
        ----------
        result : AnalysisResult
            Result from timeseries analysis.
        xlabel : str, default="Frame"
            X-axis label.
        ylabel : str, default="Contact Count"
            Y-axis label.
        title : str, default="Contacts Over Time"
            Plot title.
        figsize : tuple of (float, float), default=(10, 4)
            Figure dimensions (width, height).
        ax : Axes, optional
            Existing axes to plot on.
        time_units : str, optional
            Time unit for x-axis (e.g., "ns", "us").
        dt : float, default=1.0
            Time step multiplier when using time_units.
        legend : bool, default=True
            Whether to show legend.
        max_series : int, default=10
            Maximum number of residue series to plot.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data from result
        frames = result.data["frames"]
        contact_counts = result.data["contact_counts"]

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Convert frames to time if units specified
        if time_units:
            x_values = np.array(frames) * dt
            xlabel = f"Time ({get_unit_label(time_units)})"
        else:
            x_values = np.array(frames)

        # Plot multiple series
        colors = COLOR_SCALES["categorical"]
        series_items = list(contact_counts.items())[:max_series]

        for i, (label, series_values) in enumerate(series_items):
            line_color = colors[i % len(colors)]
            ax.plot(
                x_values,
                series_values,
                color=line_color,
                linewidth=1.2,
                label=str(label),
            )

        if legend and len(series_items) <= 10:
            ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")

        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_xlim(x_values[0], x_values[-1])
        ax.set_ylim(0, None)

        plt.tight_layout()
        return fig, ax


class DistanceTimeSeriesPlotter(BasePlotter):
    """Plotter for distance time series.

    Displays distance evolution between query and database residues
    over trajectory frames with optional cutoff line and contact highlighting.

    See Also
    --------
    DistanceAnalysis : Generates distance data
    TimeSeriesPlotter : Contact count time series
    """

    name = "distance_timeseries"
    required_analysis = "distances"
    description = "Distance over time between query and database residue"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required distance data."""
        required_keys = ["frames", "distances"]
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
        cutoff: Optional[float] = None,
        xlabel: str = "Frame",
        ylabel: str = "Distance (Å)",
        title: str = "Distance Over Time",
        figsize: Tuple[float, float] = (12, 4),
        ax: Optional[Axes] = None,
        time_units: Optional[str] = None,
        dt: float = 1.0,
    ) -> Tuple[Figure, Axes]:
        """Create distance time series plot.

        Parameters
        ----------
        result : AnalysisResult
            Result from distances analysis.
        cutoff : float, optional
            Distance cutoff to draw as horizontal line.
        xlabel : str, default="Frame"
            X-axis label.
        ylabel : str, default="Distance (Å)"
            Y-axis label.
        title : str, default="Distance Over Time"
            Plot title.
        figsize : tuple of (float, float), default=(12, 4)
            Figure dimensions (width, height).
        ax : Axes, optional
            Existing axes to plot on.
        time_units : str, optional
            Time unit for x-axis (e.g., "ns", "us").
        dt : float, default=1.0
            Time step multiplier when using time_units.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        # Extract data from result
        frames = result.data["frames"]
        distances = result.data["distances"]
        min_distances = result.data.get("min_distances")
        contact_frames = result.data.get("contact_frames", [])

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Convert to time if specified
        if time_units:
            x_values = np.array(frames) * dt
            xlabel = f"Time ({get_unit_label(time_units)})"
        else:
            x_values = np.array(frames)

        # Plot COM distance
        ax.plot(
            x_values,
            distances,
            color=COLORS["data"]["query"],
            linewidth=1.2,
            label="COM Distance",
            alpha=0.8,
        )

        # Plot min distance if provided
        if min_distances is not None:
            ax.plot(
                x_values,
                min_distances,
                color=COLORS["data"]["database"],
                linewidth=1.2,
                label="Min Distance",
                alpha=0.8,
            )

        # Highlight contact frames
        if contact_frames:
            frame_set = set(contact_frames)
            contact_indices = [i for i, f in enumerate(frames) if f in frame_set]
            contact_x = x_values[contact_indices]
            contact_y = [distances[i] for i in contact_indices]

            ax.scatter(
                contact_x,
                contact_y,
                color=COLORS["data"]["highlight"],
                s=10,
                alpha=0.5,
                zorder=5,
                label="Contact",
            )

        # Draw cutoff line
        if cutoff is not None:
            ax.axhline(
                y=cutoff,
                color=COLORS["error"]["main"],
                linestyle="--",
                linewidth=1,
                alpha=0.7,
                label=f"Cutoff ({cutoff} Å)",
            )

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_xlim(x_values[0], x_values[-1])
        ax.set_ylim(0, None)

        plt.tight_layout()
        return fig, ax


class ContactEventsPlotter(BasePlotter):
    """Plotter for contact event timelines.

    Displays contact events as horizontal bars on a timeline,
    showing when contacts occur during the trajectory.

    See Also
    --------
    KineticsAnalysis : Generates contact event data
    SurvivalCurvePlotter : Survival probability curves
    """

    name = "contact_events"
    required_analysis = "kinetics"
    description = "Contact events timeline showing when contacts occur"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains contact_frames data."""
        if "contact_frames" not in result.data:
            raise ValueError(
                f"AnalysisResult missing 'contact_frames' key for {cls.name}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        xlabel: str = "Frame",
        title: str = "Contact Events",
        figsize: Tuple[float, float] = (12, 2),
        ax: Optional[Axes] = None,
        time_units: Optional[str] = None,
        dt: float = 1.0,
    ) -> Tuple[Figure, Axes]:
        """Create contact events timeline plot.

        Parameters
        ----------
        result : AnalysisResult
            Result from kinetics analysis.
        xlabel : str, default="Frame"
            X-axis label.
        title : str, default="Contact Events"
            Plot title.
        figsize : tuple of (float, float), default=(12, 2)
            Figure dimensions (width, height).
        ax : Axes, optional
            Existing axes to plot on.
        time_units : str, optional
            Time unit for x-axis (e.g., "ns", "us").
        dt : float, default=1.0
            Time step multiplier when using time_units.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.
        """
        cls.validate_result(result)
        apply_prolint_style()

        contact_frames = result.data["contact_frames"]
        n_frames = result.metadata.get("n_frames", 0)
        kinetics = result.data.get("kinetics", {})
        occupancy = kinetics.get("occupancy", 0)
        n_events = kinetics.get("n_events", 0)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        if not contact_frames:
            ax.text(
                0.5,
                0.5,
                "No contact events",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=10,
                color=COLORS["text"]["secondary"],
            )
            ax.set_xlim(0, n_frames if n_frames > 0 else 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return fig, ax

        # Convert to time if units specified
        scale = dt if time_units else 1.0
        if time_units:
            xlabel = f"Time ({get_unit_label(time_units)})"

        # Find continuous segments and draw bars
        sorted_frames = sorted(contact_frames)
        start = sorted_frames[0]
        prev = sorted_frames[0]

        for frame in sorted_frames[1:]:
            if frame != prev + 1:
                # Draw bar for completed segment
                ax.barh(
                    0.5,
                    (prev - start + 1) * scale,
                    left=start * scale,
                    height=0.6,
                    color=COLORS["data"]["query"],
                    edgecolor=COLORS["data"].get("query_dark", COLORS["data"]["query"]),
                    linewidth=0.5,
                )
                start = frame
            prev = frame

        # Draw last segment
        ax.barh(
            0.5,
            (prev - start + 1) * scale,
            left=start * scale,
            height=0.6,
            color=COLORS["data"]["query"],
            edgecolor=COLORS["data"].get("query_dark", COLORS["data"]["query"]),
            linewidth=0.5,
        )

        ax.set_xlim(0, n_frames * scale)
        ax.set_ylim(0, 1)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_yticks([])
        ax.set_title(title, fontsize=11, fontweight="medium")

        # Add stats
        stats_text = f"{n_events} events, {occupancy:.1%} occupancy"
        ax.text(
            0.98,
            0.98,
            stats_text,
            transform=ax.transAxes,
            fontsize=8,
            ha="right",
            va="top",
            color=COLORS["text"]["secondary"],
        )

        plt.tight_layout()
        return fig, ax


# Register plotters
PlottingRegistry.register("timeseries", TimeSeriesPlotter)
PlottingRegistry.register("distance_timeseries", DistanceTimeSeriesPlotter)
PlottingRegistry.register("contact_events", ContactEventsPlotter)
