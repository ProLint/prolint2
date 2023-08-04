import os
import numpy as np
import matplotlib.pyplot as plt
import inspect
from .utils import use_1d_script_template, get_metric_list_by_residues


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

        if self.fig_close:
            plt.close(fig)

    def generate_script(self, class_code, script_filename):
        # create a python script that generates the plot
        plotting_function_source = inspect.getsource(class_code)
        if self.__class__.__name__ == 'ResiduePlot':
            tail = """
    # Generate the plot
    PLOT = {}(u, mean_contacts, lipid='CHOL', metric_name='MeanMetric')
    PLOT.save_plot(show=False)
            """.format(self.__class__.__name__)
        script_code = use_1d_script_template(plotting_function_source, tail)

        with open(script_filename, 'w') as f:
            f.write(script_code)


class ResiduePlot(Plotter1D):
    def __init__(self, universe, metric, lipid=None, metric_name=None, xlabel=None, ylabel=None, fn=None, title=None, fig_close=False, generate_script=False):
        super().__init__(xlabel, ylabel, fn, title, fig_close)
        self.resids = universe.query.residues.resids
        self.metric = get_metric_list_by_residues(universe, metric, lipid, metric_name)

        
    def create_plot(self):    
        """Plot interactions as a function of residue index"""
        
        bar_color = "#176BA0"
        if self.ylabel is None:
            self.ylabel = "Interactions"

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "Figure_interactions.pdf")

        # plot
        size_mapping = {
            20: (2.8, 1.5),
            50: (3.2, 1.5),
            300: (3.8, 1.8),
            1000: (4.5, 1.8),
            2000: (6.0, 1.8),
            np.inf: (7.5, 1.8)
        }

        size = next(size for limit, size in size_mapping.items() if len(self.resids) <= limit)
        figsize = size_mapping.get(size, (7.5, 1.8))

        self.figure, self.axis = plt.subplots(1, 1, figsize=figsize)

        self.axis.bar(self.resids, self.metric, linewidth=0, color=bar_color)

        self.axis.set_ylim(0, self.metric.max() * 1.05)
        self.axis.set_xlim(self.resids.min() - 1, self.resids.max() + 1)
        self.axis.set_ylabel(self.ylabel, fontsize=8, weight="bold")
        self.axis.set_xlabel("Residue Index", fontsize=8, weight="bold")

        for label in self.axis.xaxis.get_ticklabels() + self.axis.yaxis.get_ticklabels():
            label.set_fontsize(8)
            label.set_weight("bold")

        if self.title is not None:
            self.axis.set_title(self.title, fontsize=8, weight="bold")

        plt.tight_layout()


    



