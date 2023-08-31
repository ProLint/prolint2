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
]


class MultiplePointDistribution(Plotter):
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
                "point_distribution_{}_{}.pdf".format(lipid_type, metric_name),
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
        # Get a list of metrics for the given lipid and metric_name
        res_list, metric_list = get_metric_list_by_residues(
            self.universe, self.metric, lipid_type, metric_name, res_list=res_ids
        )

        # Calculate the number of variables and the corresponding angles for the radar chart
        num_vars = len(metric_list)
        theta = np.linspace(0, 2 * np.pi - np.pi / 6, num_vars, endpoint=False)
        theta += np.pi / 12

        axis.set_theta_zero_location("S")
        axis.set_theta_direction(-1)
        magic_number = len(metrics) // 32 + 1
        cmap = mpl.cm.get_cmap(palette)
        rescale = lambda metrics: (metrics - np.min(metrics)) / (
            np.max(metrics) - np.min(metrics)
        )

        # Plot the radar chart bars
        axis = axis.bar(
            theta, metrics, width=0.01, alpha=0.5, color=cmap(rescale(metrics))
        )

        # Customize x-axis ticks and labels
        resnames = [
            self.universe.query.residues.resnames[j]
            for j in [
                self.universe.query.residues.resids.tolist().index(i)
                for i in resids
            ]
        ]
        axis.set_xticks(theta[::magic_number])
        axis.set_xticklabels(
            [
                "{} {}".format(
                    resnames[i],
                    resids[i],
                )
                for i in range(0, len(resids), magic_number)
            ],
            fontfamily="Arial Rounded MT Bold",
        )
        axis.tick_params(axis="x", pad=20, labelsize=12)
        axis.xaxis.grid(True, linestyle="dashed", alpha=0.5)

        return axis

    def create_plot(self, rows="metrics", res_ids=None, decimals=1, **kwargs):
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
                n_rows, n_cols, figsize=self.fig_size, sharex=True, sharey="row", subplot_kw={"polar": True}
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
                n_rows, n_cols, figsize=self.fig_size, sharex=True, sharey="col", subplot_kw={"polar": True}
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

        # # format decimal places in yaxis labels
        # for ax in axes.flat:
        #     ax.yaxis.set_major_formatter(FormatStrFormatter("%.{}f".format(decimals)))

        # # adding ylabel to the first column
        # if n_rows == 1:
        #     axes[0].set_ylabel(self.ylabel, fontsize=12, fontfamily="Arial Unicode MS")
        # else:
        #     for ax in axes[:, 0]:
        #         ax.set_ylabel(self.ylabel, fontsize=12, fontfamily="Arial Unicode MS")

        # # Set default filename, and title if not provided
        # if self.fn is None:
        #     self.fn = os.path.join(
        #         os.getcwd(),
        #         "point_distribution_{}_{}.pdf".format(lipid_type, metric_name),
        #     )
        # if self.title is None:
        #     if rows == "metrics":
        #         self.title = "Point distributions - Metrics x Lipid types"
        #     else:
        #         self.title = "Point distributions - Lipid types x Metrics"

        # # # Add title to the plot
        # fig.suptitle(
        #     self.title,
        #     fontsize=14,
        #     weight="bold",
        #     fontfamily="Arial Rounded MT Bold",
        # )
        plt.tight_layout()