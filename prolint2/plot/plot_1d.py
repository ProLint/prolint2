import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import inspect
from .utils import *
import logomaker


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
        if self.__class__.__name__ in ["ResiduePlot", "ResidueLogo"]:
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


class ResiduePlot(Plotter1D):
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
            np.inf: (7.5, 1.8),
        }

        size = next(
            size for limit, size in size_mapping.items() if len(self.resids) <= limit
        )
        figsize = size_mapping.get(size, (7.5, 1.8))

        self.figure, self.axis = plt.subplots(1, 1, figsize=figsize)

        self.axis.bar(self.resids, self.metric, linewidth=0, color=bar_color)

        self.axis.set_ylim(0, self.metric.max() * 1.05)
        self.axis.set_xlim(self.resids.min() - 1, self.resids.max() + 1)
        self.axis.set_ylabel(self.ylabel, fontsize=8, weight="bold")
        self.axis.set_xlabel("Residue Index", fontsize=8, weight="bold")

        for label in (
            self.axis.xaxis.get_ticklabels() + self.axis.yaxis.get_ticklabels()
        ):
            label.set_fontsize(8)
            label.set_weight("bold")

        if self.title is not None:
            self.axis.set_title(self.title, fontsize=8, weight="bold")

        plt.tight_layout()


class ResidueLogo(Plotter1D):
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
        self.logos = universe.query.residues.resnames
        self.metric = get_metric_list_by_residues(universe, metric, lipid, metric_name)

    def create_plot(self):
        letter_map = {
            "CYS": "C",
            "ASP": "D",
            "SER": "S",
            "GLN": "Q",
            "LYS": "K",
            "ILE": "I",
            "PRO": "P",
            "THR": "T",
            "PHE": "F",
            "ASN": "N",
            "GLY": "G",
            "HIS": "H",
            "LEU": "L",
            "ARG": "R",
            "TRP": "W",
            "ALA": "A",
            "VAL": "V",
            "GLU": "E",
            "TYR": "Y",
            "MET": "M",
        }

        color_scheme = "chemistry"

        logos_checked = []
        for name in self.logos:
            if len(name) == 1:
                logos_checked.append(name)
            else:
                logos_checked.append(letter_map[name])
        if self.ylabel is None:
            self.ylabel = "Interactions"
        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "Figure_interactions_logo.pdf")

        length = 100
        gap = 1000
        # check for chain breaks, gray_areas and axis breaks
        axis_obj = AxisIndex(self.resids, logos_checked, self.metric, length, gap)
        axis_obj.sort()
        # plot
        for page_idx in axis_obj.breaks.keys():
            n_rows = len(axis_obj.breaks[page_idx])
            fig, axes = plt.subplots(
                n_rows, 1, figsize=(4.5, 1.3 * n_rows), sharey=True
            )
            plt.subplots_adjust(hspace=0.5, left=0.2)
            ymax = []
            for ax_idx, ax in enumerate(np.atleast_1d(axes)):
                resi_selected = [item[0] for item in axis_obj.breaks[page_idx][ax_idx]]
                logos_selected = [item[1] for item in axis_obj.breaks[page_idx][ax_idx]]
                interaction_selected = [
                    item[2] for item in axis_obj.breaks[page_idx][ax_idx]
                ]
                ymax.append(np.max(interaction_selected))
                if np.sum(interaction_selected) > 0:
                    df = pd.DataFrame(
                        {
                            "Resid": resi_selected,
                            "Resn": logos_selected,
                            "Data": interaction_selected,
                        }
                    )
                    matrix = df.pivot(
                        index="Resid", columns="Resn", values="Data"
                    ).fillna(0)
                    logomaker.Logo(matrix, color_scheme=color_scheme, ax=ax)
                if ax_idx == (n_rows - 1):
                    ax.set_xlabel("Residue Index", fontsize=8, weight="bold")
                ax.xaxis.set_major_locator(MultipleLocator(20))
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                ax.set_xlim(resi_selected[0] - 0.5, resi_selected[-1] + 0.5)
                ax.set_ylabel(self.ylabel, fontsize=8, weight="bold", va="center")
                for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
                    plt.setp(label, fontsize=8, weight="bold")
            np.atleast_1d(axes)[-1].set_ylim(0, np.max(ymax) * 1.05)
            # plot missing areas
            if page_idx in axis_obj.gray_areas.keys():
                for item in axis_obj.gray_areas[page_idx]:
                    np.atleast_1d(axes)[item[0]].axvspan(
                        item[1], item[2], facecolor="#c0c0c0", alpha=0.3
                    )
            if self.title is not None:
                np.atleast_1d(axes)[0].set_title(self.title, fontsize=10, weight="bold")
            plt.tight_layout()
            if len(axis_obj.breaks.keys()) == 1:
                fig.savefig(self.fn, dpi=300)
            else:
                name, ext = os.path.splitext(self.fn)
                fig.savefig("{}_{}{}".format(name, page_idx, ext), dpi=300)
