import os
import numpy as np
import matplotlib.pyplot as plt


##### ALL THE PLOTS ARE GOING TO BE GENERATED DIRECTLY FROM A UNIVERSE AND A CONTACTS INSTANCE 
# THAT CAN BE LOADED FROM A CSV FILE OR CALCULATED FROM SCRATCH #####

def use_1d_script_template(xs, ys, code_body, fn=None):
    template = """
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-

    import os
    import numpy as np
    import matplotlib.pyplot as plt

    {} 
    """.format(code_body)

def get_metric_list_by_residues(universe, metric, lipid=None, metric_name=None):
    """Get a list of metric values for a given lipid and metric name"""
    resids = universe.residues.resids
    metric_list = []
    for resi in resids:
        if resi in metric.keys():
            if metric_name:
                if lipid:
                    try:
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        metric_list.append(0)
                else:
                    # get the lipid names available
                    lipids_dict = {i: lipid for i, lipid in enumerate(universe.database.unique_resnames)}
                    lipid = int(input("Please enter the lipid name that you want to analyze: \n {}.".format(lipids_dict)))
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available".format(lipid))
                    if metric_name not in metric[resi][lipid].keys():
                        # break the loops and print error message
                        print("The metric {} is not available for the lipid {}.".format(metric_name, lipid))
                        break
                    try:
                        metric_list.append(metric[resi][lipid][metric_name])
                    except KeyError:
                        print("The metric {} is not available for the lipid {}.".format(metric_name, lipid))
            else:
                if lipid:
                    if lipid not in universe.database.unique_resnames:
                        # break the loops and print error message
                        print("The lipid {} is not available.".format(lipid))
                        break
                else:
                    # get the lipid names available
                    lipids_dict = {i: lipid for i, lipid in enumerate(universe.database.unique_resnames)}
                    lipid = int(input("Please enter the lipid name that you want to analyze: \n {}.".format(lipids_dict)))
                    try:
                        lipid = lipids_dict[lipid]
                    except KeyError:
                        print("The lipid {} is not available.".format(lipid))
                # get the metric names available
                metric_names = list(metric[resi][lipid].keys())
                inc = 1
                while len(metric_names) == 0:
                    metric_names = list(metric[resi+inc][lipid].keys())
                    inc += 1
                metric_names_dict = {i: metric_name for i, metric_name in enumerate(metric_names)}
                metric_name = int(input("Please enter the metric name that you want to analyze: \n {}.".format(metric_names_dict)))
                try:
                    metric_name = metric_names_dict[metric_name]
                except KeyError:
                    print("The metric {} is not available.".format(metric_name))
                try:
                    metric_list.append(metric[resi][lipid][metric_name])
                except KeyError:
                    metric_list.append(0)
        else:
            metric_list.append(0)
    return np.array(metric_list)

class ResiduePlotter:
    def __init__(self, universe, metric, ylabel=None, fn=None, title=None, fig_close=False):
        self.resids = universe.query.residues.resids
        self.metric = metric
        self.ylabel = ylabel
        self.fn = fn
        self.title = title
        self.fig_close = fig_close

    def plot_residue_data(self):
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

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        ax.bar(self.resids, self.metric, linewidth=0, color=bar_color)

        ax.set_ylim(0, self.metric.max() * 1.05)
        ax.set_xlim(self.resids.min() - 1, self.resids.max() + 1)
        ax.set_ylabel(self.ylabel, fontsize=8, weight="bold")
        ax.set_xlabel("Residue Index", fontsize=8, weight="bold")

        for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
            label.set_fontsize(8)
            label.set_weight("bold")

        if self.title is not None:
            ax.set_title(self.title, fontsize=8, weight="bold")

        plt.tight_layout()

        fig.savefig(self.fn, dpi=300)

        if self.fig_close:
            plt.close(fig)

        return

    def generate_plot_residue_data(self):
        # create a python script that generates the plot
        body = """
        

        """







