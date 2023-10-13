r""":mod:`prolint2.plotting.multiples`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os
import math
import pandas as pd
import numpy as np
import networkx as nx
import logomaker as lm
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FormatStrFormatter
from prolint2.computers.distances import SerialDistances
from prolint2.computers.distances import TwoPointDistances
from matplotlib.cm import ScalarMappable
import inspect
from scipy import interpolate
from prolint2.server.chord_utils import contact_chord
from .utils import *
from .plotting import Plotter

## seaborn config for paper quality plots
import seaborn as sns

sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")

__all__ = [
    "MultiplePointDistribution",
    "MultipleRadar",
    "MultipleMosaics",
]


class MultiplePointDistribution(Plotter):
    """
    Initialize the MultiplePointDistribution instance.

    Parameters
    ----------
    universe : Universe
        The universe containing the data for the plot.
    metric : dict
        A dictionary containing the metric data.
    xlabel : str, optional
        The label for the x-axis. Default is None.
    ylabel : str, optional
        The label for the y-axis. Default is None.
    fn : str, optional
        The filename to save the plot. Default is None.
    title : str, optional
        The title for the plot. Default is None.
    fig_size : tuple, optional
        The figure size as a tuple (width, height). Default is (8, 8).

    Returns
    -------
    None
    """

    def __init__(
        self,
        universe,
        metric,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        # Initialize the PointDistribution instance.
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def _plot_single(self, lipid_type, metric_name, axis, res_ids=None, **kwargs):
        """
        Plot a scatter plot for a specific lipid type and metric.

        Parameters
        ----------
        lipid_type : str
            The lipid type to plot.
        metric_name : str
            The name of the metric to plot.
        axis : matplotlib axis
            The axis for the scatter plot.
        res_ids : list, optional
            List of residue IDs to include. Default is None.
        **kwargs
            Additional keyword arguments for the scatterplot.

        Returns
        -------
        axis : matplotlib axis
            The updated axis containing the scatter plot.
        """
        # Get the list of metric values for specified lipid_type and metric_name
        res_list, metric_list = get_metric_list_by_residues(
            self.universe, self.metric, lipid_type, metric_name, res_list=res_ids
        )

        # Create scatter plot without palette
        axis = sns.scatterplot(x=res_list, y=metric_list, ax=axis, **kwargs)

        axis.tick_params(axis="both", which="major", labelsize=10)
        axis.set_title(
            "{} | {}".format(lipid_type, metric_name),
            fontsize=11,
            fontfamily="Arial Unicode MS",
        )
        # add axes labels
        axis.set_xlabel(self.xlabel, fontsize=12, fontfamily="Arial Unicode MS")

        return axis

    def create_plot(self, rows="metrics", res_ids=None, decimals=1, **kwargs):
        """
        Create the multiple point distribution plot.

        Parameters
        ----------
        rows : str
            Specify whether to arrange plots by 'metrics' or 'lipids'.
        res_ids : list, optional
            List of residue IDs to include in the plot. Default is None.
        decimals : int, optional
            Number of decimal places to format y-axis labels. Default is 1.
        **kwargs
            Additional keyword arguments for the plot.

        Returns
        -------
        None
        """
        if rows not in ["metrics", "lipids"]:
            raise ValueError("rows must be either 'metrics' or 'lipids'.")

        n_lipids = len(self.universe.database.unique_resnames)
        n_metrics = len(
            self.metric[list(self.metric.keys())[0]][
                list(self.metric[list(self.metric.keys())[0]].keys())[0]
            ]
        )
        if n_lipids == 1 and n_metrics == 1:
            raise ValueError(
                "You only have one lipid type and one metric. We encourage you to use a single plot instead."
            )

        metric_names = []
        # Collect all unique metric names from the metric dictionary
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names:
                        metric_names.append(key)

        if rows == "metrics":
            n_rows = n_metrics
            n_cols = n_lipids
        else:
            n_rows = n_lipids
            n_cols = n_metrics

        # set default labels
        if self.xlabel is None:
            self.xlabel = "Residue ID"
        if self.ylabel is None:
            self.ylabel = "Metric Value"

        if rows == "metrics":
            fig, axes = plt.subplots(
                n_rows, n_cols, figsize=self.fig_size, sharex=True, sharey="row"
            )
            if n_rows == 1 and n_cols > 1:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    axes[i] = self._plot_single(
                        lipid_type, metric_names[0], axes[i], res_ids=res_ids, **kwargs
                    )
            elif n_rows > 1 and n_cols == 1:
                for i, metric_name in enumerate(metric_names):
                    axes[i] = self._plot_single(
                        self.universe.database.unique_resnames[0],
                        metric_name,
                        axes[i],
                        res_ids=res_ids,
                        **kwargs
                    )
            else:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    for j, metric_name in enumerate(metric_names):
                        axes[j, i] = self._plot_single(
                            lipid_type,
                            metric_name,
                            axes[j, i],
                            res_ids=res_ids,
                            **kwargs
                        )
        elif rows == "lipids":
            fig, axes = plt.subplots(
                n_rows, n_cols, figsize=self.fig_size, sharex=True, sharey="col"
            )
            if n_rows == 1 and n_cols > 1:
                for i, metric_name in enumerate(metric_names):
                    axes[i] = self._plot_single(
                        self.universe.database.unique_resnames[0],
                        metric_name,
                        axes[i],
                        res_ids=res_ids,
                        **kwargs
                    )
            elif n_rows > 1 and n_cols == 1:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    axes[i] = self._plot_single(
                        lipid_type, metric_names[0], axes[i], res_ids=res_ids, **kwargs
                    )
            else:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    for j, metric_name in enumerate(metric_names):
                        axes[i, j] = self._plot_single(
                            lipid_type,
                            metric_name,
                            axes[i, j],
                            res_ids=res_ids,
                            **kwargs
                        )

        # format decimal places in yaxis labels
        for ax in axes.flat:
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.{}f".format(decimals)))

        # adding ylabel to the first column
        if n_rows == 1:
            axes[0].set_ylabel(self.ylabel, fontsize=12, fontfamily="Arial Unicode MS")
        else:
            for ax in axes[:, 0]:
                ax.set_ylabel(self.ylabel, fontsize=12, fontfamily="Arial Unicode MS")

        # Set default filename, and title if not provided
        if self.fn is None:
            self.fn = os.path.join(
                os.getcwd(),
                "multiple_point_distribution.pdf",
            )
        if self.title is None:
            if rows == "metrics":
                self.title = "Point distributions - Metrics x Lipid types"
            else:
                self.title = "Point distributions - Lipid types x Metrics"

        # # Add title to the plot
        fig.suptitle(
            self.title,
            fontsize=14,
            weight="bold",
            fontfamily="Arial Rounded MT Bold",
        )
        plt.tight_layout()


class MultipleRadar(Plotter):
    """
    Initialize a MultipleRadar instance.

    Parameters
    ----------
    universe : Universe
        The molecular dynamics universe.
    metric : dict
        A dictionary containing metrics for lipids.
    xlabel : str, optional
        The label for the x-axis.
    ylabel : str, optional
        The label for the y-axis.
    fn : str, optional
        The filename to save the plot.
    title : str, optional
        The title of the radar chart.
    fig_size : tuple, optional
        The figure size in inches (width, height).

    Inherits from Plotter and sets plot labels, title, and figure size.
    """

    def __init__(
        self,
        universe,
        metric,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        # Initialize the PointDistribution instance.
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def _plot_single(
        self, lipid_type, metric_name, axis, theta, res_ids=None, **kwargs
    ):
        """
        Plot a radar chart for a single lipid type and metric.

        Parameters
        ----------
        lipid_type : str
            The lipid type.
        metric_name : str
            The metric name.
        axis : matplotlib.axes._subplots.PolarAxesSubplot
            The polar axes for plotting.
        theta : numpy.ndarray
            The angles for the radar chart.
        res_ids : list, optional
            A list of residue IDs to include in the plot.
        kwargs : dict
            Additional keyword arguments for plotting.

        Returns
        -------
        matplotlib.axes._subplots.PolarAxesSubplot
            The polar axes with the radar chart.
        """
        # Get a list of metrics for the given lipid and metric_name
        res_list, metric_list = get_metric_list_by_residues(
            self.universe, self.metric, lipid_type, metric_name, res_list=res_ids
        )

        axis.set_theta_zero_location("S")
        axis.set_theta_direction(-1)

        # Plot the radar chart bars
        axis.bar(theta, metric_list, width=0.02, alpha=0.5, **kwargs)

        axis.set_xlim(0, 2 * np.pi)
        axis.xaxis.grid(True, linestyle="dashed", alpha=0.5)
        axis.set_yticks([max(metric_list) / 2, max(metric_list)])
        axis.set_yticklabels([])

        axis.tick_params(axis="x", pad=10, labelsize=8)
        axis.yaxis.grid(True, linestyle="dashed", alpha=0.5)
        axis.set_rlabel_position(0)
        axis.set_thetamin(15)
        axis.set_thetamax(345)

        axis.set_title(
            "{} | {}".format(lipid_type, metric_name),
            fontsize=11,
            pad=10,
            fontfamily="Arial Rounded MT Bold",
        )

        return axis

    def create_plot(
        self, rows="metrics", res_ids=None, decimals=1, magic_number=None, **kwargs
    ):
        """
        Create a radar chart plot.

        Parameters
        ----------
        rows : str, optional
            Determines the layout of the plot. Must be 'metrics' or 'lipids'.
        res_ids : list, optional
            A list of residue IDs to include in the plot.
        decimals : int, optional
            The number of decimal places for labels.
        magic_number : int, optional
            A factor to determine the number of variables to display.
        kwargs : dict
            Additional keyword arguments for plotting.

        Raises
        ------
        ValueError
            If rows is neither 'metrics' nor 'lipids'.

        """
        if rows not in ["metrics", "lipids"]:
            raise ValueError("rows must be either 'metrics' or 'lipids'.")

        n_lipids = len(self.universe.database.unique_resnames)
        n_metrics = len(
            self.metric[list(self.metric.keys())[0]][
                list(self.metric[list(self.metric.keys())[0]].keys())[0]
            ]
        )
        if n_lipids == 1 and n_metrics == 1:
            raise ValueError(
                "You only have one lipid type and one metric. We encourage you to use a single plot instead."
            )

        metric_names = []
        # Collect all unique metric names from the metric dictionary
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names:
                        metric_names.append(key)

        if rows == "metrics":
            n_rows = n_metrics
            n_cols = n_lipids
        else:
            n_rows = n_lipids
            n_cols = n_metrics

        resids_total_list = self.universe.query.residues.resids.tolist()
        if res_ids is None:
            res_list = resids_total_list
        else:
            res_list = res_ids

        resnames = [
            self.universe.query.residues.resnames[j]
            for j in [resids_total_list.index(i) for i in res_list]
        ]

        # Calculate the number of variables and the corresponding angles for the radar chart
        num_vars = len(res_list)
        theta = np.linspace(0, 2 * np.pi - np.pi / 6, num_vars, endpoint=False)
        theta += np.pi / 12

        if magic_number is None:
            magic_number = len(res_list) // 15 + 1

        if rows == "metrics":
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=self.fig_size,
                sharey="row",
                subplot_kw={"polar": True},
            )
            if n_rows == 1 and n_cols > 1:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    axes[i] = self._plot_single(
                        lipid_type,
                        metric_names[0],
                        axes[i],
                        theta,
                        res_ids=res_ids,
                        **kwargs
                    )
            elif n_rows > 1 and n_cols == 1:
                for i, metric_name in enumerate(metric_names):
                    axes[i] = self._plot_single(
                        self.universe.database.unique_resnames[0],
                        metric_name,
                        axes[i],
                        theta,
                        res_ids=res_ids,
                        **kwargs
                    )
            else:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    for j, metric_name in enumerate(metric_names):
                        axes[j, i] = self._plot_single(
                            lipid_type,
                            metric_name,
                            axes[j, i],
                            theta,
                            res_ids=res_ids,
                            **kwargs
                        )
        elif rows == "lipids":
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=self.fig_size,
                sharey="col",
                subplot_kw={"polar": True},
            )
            if n_rows == 1 and n_cols > 1:
                for i, metric_name in enumerate(metric_names):
                    axes[i] = self._plot_single(
                        self.universe.database.unique_resnames[0],
                        metric_name,
                        axes[i],
                        theta,
                        res_ids=res_ids,
                        **kwargs
                    )
            elif n_rows > 1 and n_cols == 1:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    axes[i] = self._plot_single(
                        lipid_type,
                        metric_names[0],
                        axes[i],
                        theta,
                        res_ids=res_ids,
                        **kwargs
                    )
            else:
                for i, lipid_type in enumerate(self.universe.database.unique_resnames):
                    for j, metric_name in enumerate(metric_names):
                        axes[i, j] = self._plot_single(
                            lipid_type,
                            metric_name,
                            axes[i, j],
                            theta,
                            res_ids=res_ids,
                            **kwargs
                        )

        for ax in axes.flat:
            ax.set_xticks(theta[::magic_number])
            ax.set_xticklabels(
                [
                    "{} {}".format(
                        resnames[i],
                        res_list[i],
                    )
                    for i in range(0, len(res_list), magic_number)
                ],
                fontfamily="Arial Unicode MS",
            )

        for ax in axes.flat:
            yticks = ax.get_yticks()
            for loc in yticks:
                ax.text(
                    0.0,
                    loc,
                    str(round(float(loc), decimals)),
                    ha="center",
                    va="top",
                    fontsize=8,
                    fontfamily="Arial Unicode MS",
                )

        # Set default filename, and title if not provided
        if self.fn is None:
            self.fn = os.path.join(
                os.getcwd(),
                "multiple_radars.pdf",
            )
        if self.title is None:
            if rows == "metrics":
                self.title = "Radar - Metrics x Lipid types"
            else:
                self.title = "Radar - Lipid types x Metrics"

        # # Add title to the plot
        fig.suptitle(
            self.title,
            fontsize=14,
            weight="bold",
            fontfamily="Arial Rounded MT Bold",
        )
        plt.tight_layout()


class MultipleMosaics(Plotter):
    """
    Initialize the MultipleMosaics object with specified plot attributes.

    Parameters
    ----------
    list_of_grids : list of str
        List of filenames containing grid data.
    xlabel : str, optional
        Label for the x-axis of the plot.
    ylabel : str, optional
        Label for the y-axis of the plot.
    fn : str, optional
        Filename to save the plot. If None, the plot won't be saved.
    title : str, optional
        Title for the plot.
    fig_size : tuple, optional
        Figure size in inches (width, height).

    Returns
    -------
    None
    """

    def __init__(
        self,
        list_of_grids,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(10, 10),
    ):
        # Initialize the MosaicsGridData object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        # load grid data from file
        self.grid_data_list = [np.loadtxt(file) for file in list_of_grids]

    def _plot_single(
        self, grid_data, axis, frame=None, title=None, v_min=None, v_max=None, **kwargs
    ):
        """
        Create a single heatmap plot.

        Parameters
        ----------
        grid_data : numpy.ndarray
            Grid data to be plotted.
        axis : matplotlib.axes.Axes
            The axis on which the heatmap will be plotted.
        frame : int, optional
            Frame parameter for the grid data.
        title : str, optional
            Title for the plot.
        v_min : float, optional
            Minimum value for color scaling.
        v_max : float, optional
            Maximum value for color scaling.
        **kwargs : dict
            Additional keyword arguments for the imshow function.

        Returns
        -------
        None
        """
        # Create heatmap using the grid data
        if frame is None:
            if v_min is not None:
                im = axis.imshow(
                    grid_data, origin="lower", vmin=v_min, vmax=v_max, **kwargs
                )
            else:
                im = axis.imshow(grid_data, origin="lower", **kwargs)
        else:
            if v_min is not None:
                im = axis.imshow(
                    grid_data[frame:-frame, frame:-frame],
                    origin="lower",
                    vmin=v_min,
                    vmax=v_max,
                    **kwargs
                )
            else:
                im = axis.imshow(
                    grid_data[frame:-frame, frame:-frame], origin="lower", **kwargs
                )

        axis.grid(False)

        # Create a colorbar of the same size as the plot
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        # Set colorbar label and font size
        cbar.ax.tick_params(labelsize=12)

        # Remove labels and ticks from the axes
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_xticklabels([])
        axis.set_yticklabels([])

        # Add title to the plot
        if title is not None:
            axis.set_title(
                title,
                fontsize=12,
                weight="bold",
                pad=15,
                fontfamily="Arial Rounded MT Bold",
            )

    def create_plot(
        self,
        shape: tuple,
        frame=None,
        nan_to_zero=False,
        titles=None,
        share_limits=False,
        **kwargs
    ):
        """
        Create a mosaic of heatmap plots.

        Parameters
        ----------
        shape : tuple
            Shape of the mosaic as (rows, columns).
        frame : int, optional
            Frame parameter for the grid data.
        nan_to_zero : bool, optional
            Convert NaN values to zero.
        titles : list of str, optional
            Titles for individual subplots. If None, no titles are used.
        share_limits : bool, optional
            Share color scale limits across subplots.
        **kwargs : dict
            Additional keyword arguments for the _plot_single function.

        Returns
        -------
        None
        """

        if len(shape) != 2:
            raise ValueError("shape must be a tuple of length 2.")
        elif shape[0] * shape[1] != len(self.grid_data_list):
            raise ValueError(
                "The number of grids must be equal to the product of the shape."
            )
        elif len(self.grid_data_list) == 1:
            raise ValueError(
                "You only have one grid. We encourage you to use a single plot instead."
            )
        if titles is not None:
            if len(titles) != len(self.grid_data_list):
                raise ValueError(
                    "titles must have the same length as the number of grids."
                )
        if nan_to_zero:
            self.grid_data_list = [np.nan_to_num(arr) for arr in self.grid_data_list]
        if share_limits:
            v_min = 0
            v_max = 0
            for grid_data in self.grid_data_list:
                v_min = min(v_min, np.nanmin(grid_data))
                v_max = max(v_max, np.nanmax(grid_data))

        n_rows = shape[0]
        n_cols = shape[1]

        # Generate a figure and axis for the plot
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=self.fig_size,
        )

        if n_rows == 1 or n_cols == 1:
            for i, grid_data in enumerate(self.grid_data_list):
                if titles is None:
                    if share_limits:
                        self._plot_single(
                            grid_data,
                            axes[i],
                            frame=frame,
                            v_min=v_min,
                            v_max=v_max,
                            **kwargs
                        )
                    else:
                        self._plot_single(grid_data, axes[i], frame=frame, **kwargs)
                else:
                    if share_limits:
                        self._plot_single(
                            grid_data,
                            axes[i],
                            frame=frame,
                            title=titles[i],
                            v_min=v_min,
                            v_max=v_max,
                            **kwargs
                        )
                    else:
                        self._plot_single(
                            grid_data, axes[i], frame=frame, title=titles[i], **kwargs
                        )
        else:
            for i, grid_data in enumerate(self.grid_data_list):
                if titles is None:
                    if share_limits:
                        self._plot_single(
                            grid_data,
                            axes[i // n_cols, i % n_cols],
                            frame=frame,
                            v_min=v_min,
                            v_max=v_max,
                            **kwargs
                        )
                    else:
                        self._plot_single(
                            grid_data,
                            axes[i // n_cols, i % n_cols],
                            frame=frame,
                            **kwargs
                        )
                else:
                    if share_limits:
                        self._plot_single(
                            grid_data,
                            axes[i // n_cols, i % n_cols],
                            frame=frame,
                            title=titles[i],
                            v_min=v_min,
                            v_max=v_max,
                            **kwargs
                        )
                    else:
                        self._plot_single(
                            grid_data,
                            axes[i // n_cols, i % n_cols],
                            frame=frame,
                            title=titles[i] ** kwargs,
                        )

        # # Set default filename and title if not provided
        # if self.fn is None:
        #     self.fn = os.path.join(
        #         os.getcwd(), "mosaics_grid_data_{}.pdf".format(prop_label)
        #     )
        # if self.title is None:
        #     self.title = "Mosaics Grid Data - {}".format(prop_label)

        # # Add title to the plot
        # ax.set_title(
        #     self.title,
        #     fontsize=14,
        #     weight="bold",
        #     pad=15,
        #     fontfamily="Arial Rounded MT Bold",
        # )
        plt.tight_layout()
