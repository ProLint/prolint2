import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import inspect
from .utils import *

## seaborn config for paper quality plots
import seaborn as sns
sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")


##### ALL THE PLOTS ARE GOING TO BE GENERATED DIRECTLY FROM A UNIVERSE AND A CONTACTS INSTANCE
# THAT CAN BE LOADED FROM A CSV FILE OR CALCULATED FROM SCRATCH #####


class Plotter1D:
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
        if self.__class__.__name__ in ["ResiduePlot", "ResidueLogo", "PointDistribution"]:
            tail = """
    # Generate the plot
    PLOT = {}(u, mean_contacts, lipid='CHOL', metric_name='MeanMetric')
    PLOT.save_plot(show=False)
            """.format(
                self.__class__.__name__
            )
        script_code = use_1d_script_template(plotting_function_source, tail)

        with open(script_filename, "w") as f:
            f.write(script_code)


class PointDistribution(Plotter1D):
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
        ax = sns.scatterplot(x=self.resids, y=self.metric, hue=self.metric, palette="flare", linewidth=0.24)
        norm = plt.Normalize(self.metric.min(), self.metric.max())
        sm = plt.cm.ScalarMappable(cmap='flare', norm=norm)
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

