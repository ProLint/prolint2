import os
from typing import Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import inspect
from .utils import *

## seaborn config for paper quality plots
import seaborn as sns

sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")

__all__ = ["Plotter", "PointDistribution", "Radar", "DensityMap"]


class Plotter:
    def __init__(self, xlabel: str = None, ylabel: str = None, fn: str = None, title: str =None, fig_size: tuple = (8, 8)):

        self.xlabel = xlabel
        self.ylabel = ylabel
        self.fn = fn
        self.title = title
        self.fig_size = fig_size

    def save_plot(self, **kwargs):
        self.create_plot(**kwargs)
        plt.savefig(self.fn, dpi=300)

    def generate_script(self, class_code, script_filename):
        # create a python script that generates the plot
        plotting_function_source = inspect.getsource(class_code)
        if self.__class__.__name__ in ["PointDistribution"]:
            tail = """
    mean_instance = MeanMetric()
    metric_instance = Metric(contacts, mean_instance)
    mean_contacts = metric_instance.compute()

    # Generate the plot
    PLOT = {}(u, mean_contacts, lipid='CHOL', metric_name='MeanMetric')
    PLOT.save_plot()
            """.format(
                self.__class__.__name__
            )
        elif self.__class__.__name__ in ["Radar"]:
            tail = """
    metric_instances_list = [MeanMetric(), SumMetric(), MaxMetric()]
    metric_instance = Metric(contacts, metric_instances_list) 
    contacts_out = metric_instance.compute()

    # Generate the plot
    PLOT = {}(contacts_out, resIDs=[2, 3, 5], lipid='POPS', metric_names=['MeanMetric', 'SumMetric', 'MaxMetric'])
    PLOT.save_plot()
            """.format(
                self.__class__.__name__
            )
        elif self.__class__.__name__ in ["DensityMap"]:
            tail = """
    # Generate the plot
    PLOT = {}(u, lipid='CHOL')
    PLOT.save_plot()
            """.format(
                self.__class__.__name__
            )

        script_code = use_1d_script_template(plotting_function_source, tail)

        with open(script_filename, "w") as f:
            f.write(script_code)


class PointDistribution(Plotter):
    def __init__(
        self,
        universe,
        metric,
        lipid=None,
        metric_name=None,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.resids = universe.query.residues.resids
        self.metric = get_metric_list_by_residues(universe, metric, lipid, metric_name)
        self.metric_name = metric_name

    def create_plot(self, **kwargs):
        """Plot the distribution of a metric for each residue."""
        fig, ax = plt.subplots(figsize=self.fig_size)

        if 'palette' in kwargs:
            ax = sns.scatterplot(
            x=self.resids,
            y=self.metric,
            hue=self.metric,
            **kwargs
            )
            norm = plt.Normalize(self.metric.min(), self.metric.max())
            sm = plt.cm.ScalarMappable(cmap=kwargs['palette'], norm=norm)
            sm.set_array([])

            # remove legend and add color bar
            ax.get_legend().remove()
            ax.figure.colorbar(sm)
        else:
            ax = sns.scatterplot(
            x=self.resids,
            y=self.metric,
            **kwargs
            )

        if self.xlabel is None:
            self.xlabel = "Residue ID"
        if self.ylabel is None:
            self.ylabel = self.metric_name
        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "point_distribution_{}.pdf".format(self.metric_name))
        if self.title is None:
            self.title = "Point distribution - {}".format(self.metric_name)

        # add labels and title
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title, fontsize=12, weight='bold')
        plt.tight_layout()

class Radar(Plotter):
    def __init__(
        self,
        metrics,
        resIDs=None,
        lipid=None,
        metric_names=None,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8)
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.metrics, self.metric_names = get_metrics_for_radar(
            metrics, metric_names, resIDs=resIDs, lipid=lipid
        )

    def create_plot(self, **kwargs):
        num_vars = len(self.metric_names)
        theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)
        theta += np.pi / 2

        fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw={"polar": True})
        ax.set_xticks(theta)
        ax.set_xticklabels(self.metric_names, fontsize=8)
        ax.tick_params(axis="x", pad=15)
        ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)

        for resi in self.metrics.keys():
            values = self.metrics[resi]
            xs = np.concatenate((theta, [theta[0]]))
            values = np.concatenate((values, [values[0]]))
            ax.plot(xs, values, label=resi, **kwargs)
            ax.fill(xs, values, alpha=0.1)

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        ax.legend(loc="upper right", bbox_to_anchor=(1.1, 1.1), title="Residue ID")

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "metrics_comparison.pdf")
        if self.title is None:
            self.title = "Metrics comparison"
        
        # add labels and title
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title, fontsize=12, weight='bold', pad=30)
        plt.tight_layout()


class DensityMap(Plotter):
    def __init__(
        self,
        universe,
        lipid=None,
        bins=150,
        size_in_mb=50000,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.lipid = lipid
        self.bins = bins
        self.size_in_mb = size_in_mb

    def create_plot(self, **kwargs):
        """Plot the preferential localization of lipids using 2D density maps."""

        # Compute the lipid coordinates
        computed_coords = compute_density(self.universe, self.lipid, self.size_in_mb)

        # Compute the lipid density
        H, xe, ye = np.histogram2d(
            computed_coords[:, 0], computed_coords[:, 1], bins=self.bins, density=True
        )

        # Generate the 2D histogram (density plot)
        fig, ax = plt.subplots(figsize=self.fig_size)

        # Plot the density map
        im = ax.imshow(
            H,
            origin="lower",
            extent=[xe[0], xe[-1], ye[0], ye[-1]],
            **kwargs
        )
        ax.grid(False)

        # Create colorbar of the same size as the plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        # increase the font size
        cbar.set_label(label="Density distribution", size=12)
        cbar.ax.tick_params(labelsize=10)
        
        # Remove labels and ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "density_map_{}.pdf".format(self.lipid))
        if self.title is None:
            self.title = "Density map - {}".format(self.lipid)
        
        # add title
        ax.set_title(self.title, fontsize=12, weight='bold', pad=20)
        plt.tight_layout()
