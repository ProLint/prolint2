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


class Plotter:
    def __init__(self, xlabel=None, ylabel=None, fn=None, title=None, fig_close=False):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.fn = fn
        self.title = title
        self.fig_close = fig_close

    def save_plot(self, show=True):
        self.create_plot()
        plt.savefig(self.fn, dpi=300)

        if show:
            plt.show()

    def generate_script(self, class_code, script_filename):
        # create a python script that generates the plot
        plotting_function_source = inspect.getsource(class_code)
        if self.__class__.__name__ in [
            "ResiduePlot",
            "ResidueLogo",
            "PointDistribution",
        ]:
            tail = """
    mean_instance = MeanMetric()
    metric_instance = Metric(contacts, mean_instance)
    mean_contacts = metric_instance.compute()

    # Generate the plot
    PLOT = {}(u, mean_contacts, lipid='CHOL', metric_name='MeanMetric')
    PLOT.save_plot(show=False)
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
    PLOT.save_plot(show=False)
            """.format(
                self.__class__.__name__
            )
        elif self.__class__.__name__ in ["DensityMap"]:
            tail = """
    # Generate the plot
    PLOT = {}(u, lipid='CHOL')
    PLOT.save_plot(show=False)
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
        fig_close=False,
        generate_script=False,
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_close)
        self.resids = universe.query.residues.resids
        self.metric = get_metric_list_by_residues(universe, metric, lipid, metric_name)

    def create_plot(self):
        ax = sns.scatterplot(
            x=self.resids,
            y=self.metric,
            hue=self.metric,
            palette="flare",
            linewidth=0.24,
        )
        norm = plt.Normalize(self.metric.min(), self.metric.max())
        sm = plt.cm.ScalarMappable(cmap="flare", norm=norm)
        sm.set_array([])

        # remove legend and add color bar
        ax.get_legend().remove()
        ax.figure.colorbar(sm)

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "Figure_interactions_points.pdf")

        # add labels and title
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)


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
        fig_close=False,
        generate_script=False,
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_close)
        self.metrics, self.metric_names = get_metrics_for_radar(
            metrics, metric_names, resIDs=resIDs, lipid=lipid
        )

    def create_plot(self):
        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "Figure_interactions_radar.pdf")

        num_vars = len(self.metric_names)
        theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)
        theta += np.pi / 2

        fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={"polar": True})
        ax.set_xticks(theta)
        ax.set_xticklabels(self.metric_names, fontsize=8)
        ax.tick_params(axis="x", pad=15)
        ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)

        for resi in self.metrics.keys():
            if sum(self.metrics[resi]) == 0:
                ax.scatter(0, 0, label=resi)
            else:
                values = self.metrics[resi]
                values = np.concatenate((values, [values[0]]))
                ax.plot(np.concatenate((theta, [theta[0]])), values, label=resi)
                ax.fill(np.concatenate((theta, [theta[0]])), values, alpha=0.1)

        ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), title="Residue ID")
        plt.tight_layout()


class DensityMap(Plotter):
    def __init__(
        self,
        universe,
        lipid=None,
        bins=150,
        size_in_mb=50000,
        colormap='viridis',
        interpolation='nearest',
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_close=False,
        generate_script=False,
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_close)
        self.universe = universe
        self.lipid = lipid
        self.bins = bins
        self.size_in_mb = size_in_mb
        self.colormap = colormap
        self.interpolation = interpolation

    def create_plot(self):
        """Plot the preferential localization of lipids using 2D density maps."""

        # Compute the lipid coordinates
        computed_coords = compute_density(self.universe, [self.lipid], self.size_in_mb) 
                    
        # Compute the lipid density
        H, xe, ye = np.histogram2d(computed_coords[:, 0], computed_coords[:, 1], bins=self.bins, density=True)

        # Generate the 2D histogram (density plot)
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot the density map
        im = ax.imshow(H, interpolation=self.interpolation, origin='lower',
                   extent=[xe[0], xe[-1], ye[0], ye[-1]],
                   cmap=self.colormap)
        ax.grid(False)

        # Add colorbar
        # Create colorbar of the same size as the plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax, label="Density")

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "Figure_density_map.pdf")

        # Remove labels and ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        plt.tight_layout()





        
