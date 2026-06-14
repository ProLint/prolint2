"""Kinetics plotters for binding dynamics visualization.

This module provides survival curves and residence time distributions
for contact kinetics analysis results.
"""

from typing import Optional, Tuple
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import COLORS, apply_prolint_style, get_unit_label


class SurvivalCurvePlotter(BasePlotter):
    """Plotter for survival curves with exponential fits.

    Visualizes contact survival probability over lag time with
    optional mono- and bi-exponential model fits.

    See Also
    --------
    KineticsAnalysis : Generates survival curve data
    ResidenceDistributionPlotter : Residence time histograms
    """

    name = "survival_curve"
    required_analysis = "kinetics"
    description = "Survival curve with mono/bi-exponential fits"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required survival curve data."""
        if "survival_curve" not in result.data:
            raise ValueError(
                f"AnalysisResult missing 'survival_curve' key for {cls.name}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )
        survival = result.data["survival_curve"]
        required_keys = ["lag_times", "survival_probability"]
        missing = [k for k in required_keys if k not in survival]
        if missing:
            raise ValueError(
                f"survival_curve missing required keys: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        xlabel: str = "Lag Time (frames)",
        ylabel: str = "Survival Probability",
        title: str = "Survival Curve",
        figsize: Tuple[float, float] = (8, 5),
        ax: Optional[Axes] = None,
        show_legend: bool = True,
        time_units: Optional[str] = None,
        dt: float = 1.0,
    ) -> Tuple[Figure, Axes]:
        """Create survival curve plot with exponential fits.

        Parameters
        ----------
        result : AnalysisResult
            Result from kinetics analysis.
        time_units : str, optional
            Time unit for x-axis (e.g., "ns", "us").
        dt : float, default=1.0
            Time step multiplier when using time_units.
        show_legend : bool, default=True
            Whether to show fit parameters in legend.
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
        survival = result.data["survival_curve"]
        kinetics = result.data.get("kinetics", {})
        n_events = kinetics.get("n_events", 0)
        min_events_mono = survival.get("min_events_mono", 5)

        # Check for insufficient data (same limit as frontend)
        if n_events < min_events_mono:
            warnings.warn(
                f"Insufficient data for survival curve plotting ({n_events} events, "
                f"need >= {min_events_mono}). Skipping plot.",
                UserWarning,
            )
            if ax is None:
                fig, ax = plt.subplots(figsize=figsize)
            else:
                fig = ax.figure
            ax.text(
                0.5,
                0.5,
                f"Insufficient data\n({n_events} events, need >= {min_events_mono})",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=11,
                color=COLORS["text"]["secondary"],
            )
            ax.set_xlabel(xlabel, fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_title(title, fontsize=12, fontweight="semibold")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return fig, ax

        lag_times = survival["lag_times"]
        survival_probability = survival["survival_probability"]
        mono_fit = survival.get("mono_fit")
        bi_fit = survival.get("bi_fit")
        selected_model = survival.get("selected_model")

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Convert to time if specified
        if time_units:
            x_values = np.array(lag_times) * dt
            xlabel = f"Lag Time ({get_unit_label(time_units)})"
        else:
            x_values = np.array(lag_times)

        # Plot raw survival curve data points
        ax.scatter(
            x_values,
            survival_probability,
            color=COLORS["neutral"]["500"],
            s=30,
            alpha=0.6,
            label="Data",
            zorder=5,
        )

        # Plot monoexponential fit if available
        if mono_fit is not None and "fitted_curve" in mono_fit:
            fit_y = mono_fit["fitted_curve"]
            is_selected = selected_model == "monoexponential"
            linewidth = 2.5 if is_selected else 1.5
            alpha = 1.0 if is_selected else 0.5
            linestyle = "-" if is_selected else "--"

            k_off = mono_fit.get("k_off", 0)
            r2 = mono_fit.get("r_squared", 0)
            half_life = mono_fit.get("half_life")
            label = f"Mono: k={k_off:.3f}, R²={r2:.3f}"
            if half_life:
                label += f", t½={half_life:.1f}"

            ax.plot(
                x_values,
                fit_y,
                color=COLORS["data"]["query"],
                linewidth=linewidth,
                alpha=alpha,
                linestyle=linestyle,
                label=label,
            )

        # Plot biexponential fit if available
        if bi_fit is not None and "fitted_curve" in bi_fit:
            fit_y = bi_fit["fitted_curve"]
            is_selected = selected_model == "biexponential"
            linewidth = 2.5 if is_selected else 1.5
            alpha = 1.0 if is_selected else 0.5
            linestyle = "-" if is_selected else "--"

            k_fast = bi_fit.get("k_fast", 0)
            k_slow = bi_fit.get("k_slow", 0)
            r2 = bi_fit.get("r_squared", 0)
            label = f"Bi: k_fast={k_fast:.3f}, k_slow={k_slow:.3f}, R²={r2:.3f}"

            ax.plot(
                x_values,
                fit_y,
                color=COLORS["data"]["database"],
                linewidth=linewidth,
                alpha=alpha,
                linestyle=linestyle,
                label=label,
            )

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")

        ax.set_ylim(0, 1.05)
        ax.set_xlim(x_values[0], x_values[-1])
        ax.grid(True, alpha=0.3, linestyle="--")

        if show_legend:
            ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

        plt.tight_layout()
        return fig, ax


class ResidenceDistributionPlotter(BasePlotter):
    """Plotter for residence time distributions.

    Visualizes the distribution of contact durations as a histogram.

    See Also
    --------
    KineticsAnalysis : Generates residence distribution data
    SurvivalCurvePlotter : Survival curves
    """

    name = "residence_distribution"
    required_analysis = "kinetics"
    description = "Histogram of residence time durations"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required residence distribution data."""
        if "residence_distribution" not in result.data:
            raise ValueError(
                f"AnalysisResult missing 'residence_distribution' key for {cls.name}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )
        residence = result.data["residence_distribution"]
        required_keys = ["bins", "counts"]
        missing = [k for k in required_keys if k not in residence]
        if missing:
            raise ValueError(
                f"residence_distribution missing required keys: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        xlabel: str = "Residence Time (frames)",
        ylabel: str = "Count",
        title: str = "Residence Time Distribution",
        figsize: Tuple[float, float] = (8, 4),
        ax: Optional[Axes] = None,
        time_units: Optional[str] = None,
        dt: float = 1.0,
        log_scale: bool = False,
    ) -> Tuple[Figure, Axes]:
        """Create residence time histogram.

        Parameters
        ----------
        result : AnalysisResult
            Result from kinetics analysis.
        time_units : str, optional
            Time unit for x-axis.
        dt : float, default=1.0
            Time step multiplier.
        log_scale : bool, default=False
            Whether to use log scale for y-axis.
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
        residence = result.data["residence_distribution"]
        bins = residence["bins"]
        counts = residence["counts"]
        kinetics = result.data.get("kinetics", {})
        mean_residence = kinetics.get("mean_residence_time", 0)
        n_events = kinetics.get("n_events", sum(counts))

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # Convert bins to time if specified
        if time_units:
            x_values = np.array(bins) * dt
            xlabel = f"Residence Time ({get_unit_label(time_units)})"
        else:
            x_values = np.array(bins)

        # Calculate bar width
        if len(x_values) > 1:
            width = (x_values[1] - x_values[0]) * 0.8
        else:
            width = 0.8

        ax.bar(
            x_values,
            counts,
            width=width,
            color=COLORS["data"]["query"],
            edgecolor=COLORS["data"]["query_dark"],
            alpha=0.8,
        )

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="semibold")

        if log_scale and max(counts) > 0:
            ax.set_yscale("log")

        ax.set_xlim(0, x_values[-1] + width if len(x_values) > 0 else 10)
        ax.grid(True, alpha=0.3, linestyle="--", axis="y")

        # Add stats using pre-computed mean from kinetics
        stats_text = f"N={n_events}, Mean={mean_residence:.1f} frames"
        ax.text(
            0.98,
            0.98,
            stats_text,
            transform=ax.transAxes,
            fontsize=9,
            ha="right",
            va="top",
            color=COLORS["text"]["secondary"],
        )

        plt.tight_layout()
        return fig, ax


# Register plotters
PlottingRegistry.register("survival_curve", SurvivalCurvePlotter)
PlottingRegistry.register("residence_distribution", ResidenceDistributionPlotter)
