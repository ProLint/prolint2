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
