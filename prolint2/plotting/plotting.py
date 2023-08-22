import os
import math
import numpy as np
import logomaker as lm
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import FancyBboxPatch
import inspect
from .utils import *

## seaborn config for paper quality plots
import seaborn as sns

sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")

__all__ = ["Plotter", "PointDistribution", "Radar", "DensityMap", "DurationGantt"]


class Plotter:
    def __init__(
        self,
        xlabel: str = None,
        ylabel: str = None,
        fn: str = None,
        title: str = None,
        fig_size: tuple = (8, 8),
    ):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.fn = fn
        self.title = title
        self.fig_size = fig_size

    def save_plot(self, **kwargs):
        self.create_plot(**kwargs)
        plt.savefig(self.fn, dpi=300, bbox_inches='tight')

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

        if "palette" in kwargs:
            ax = sns.scatterplot(
                x=self.resids, y=self.metric, hue=self.metric, **kwargs
            )
            norm = plt.Normalize(self.metric.min(), self.metric.max())
            sm = plt.cm.ScalarMappable(cmap=kwargs["palette"], norm=norm)
            sm.set_array([])

            # remove legend and add color bar
            ax.get_legend().remove()
            cbar = ax.figure.colorbar(sm)
            cbar.ax.tick_params(labelsize=12)
        else:
            ax = sns.scatterplot(x=self.resids, y=self.metric, **kwargs)

        ax.tick_params(axis="both", which="major", labelsize=12)

        if self.xlabel is None:
            self.xlabel = "Residue ID"
        if self.ylabel is None:
            self.ylabel = self.metric_name
        if self.fn is None:
            self.fn = os.path.join(
                os.getcwd(), "point_distribution_{}.pdf".format(self.metric_name)
            )
        if self.title is None:
            self.title = "Point distribution - {}".format(self.metric_name)

        # add labels and title
        plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
        plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
        plt.title(
            self.title, fontsize=14, weight="bold", fontfamily="Arial Rounded MT Bold"
        )
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
        fig_size=(8, 8),
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
        ax.set_xticklabels(self.metric_names)
        ax.tick_params(axis="x", pad=20)
        ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)

        for resi in self.metrics.keys():
            values = self.metrics[resi]
            xs = np.concatenate((theta, [theta[0]]))
            values = np.concatenate((values, [values[0]]))
            ax.plot(xs, values, label=resi, **kwargs)
            ax.fill(xs, values, alpha=0.1)

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        ax.tick_params(axis="both", which="major", labelsize=10)
        tt = ax.legend(loc="upper right", bbox_to_anchor=(1.15, 1.1), title="Residue ID", title_fontsize=12, prop={'family': "Arial Unicode MS", 'size': 12})
        tt.set_title(title="Residue ID", prop={'family': "Arial Rounded MT Bold", 'size': 12})

        if self.fn is None:
            self.fn = os.path.join(os.getcwd(), "metrics_comparison.pdf")
        if self.title is None:
            self.title = "Metrics comparison"

        # add labels and title
        plt.yticks(fontname = "Arial Unicode MS")
        plt.xticks(fontname = "Arial Rounded MT Bold")
        plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
        plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
        plt.title(self.title, fontsize=14, weight="bold", pad=30, fontfamily="Arial Rounded MT Bold")
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
            H, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]], **kwargs
        )
        ax.grid(False)

        # Create colorbar of the same size as the plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        # increase the font size
        cbar.set_label(label="Density distribution", size=12, fontfamily='Arial Rounded MT Bold')
        cbar.ax.tick_params(labelsize=12)

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
        ax.set_title(self.title, fontsize=14, weight="bold", pad=15, fontfamily="Arial Rounded MT Bold")
        plt.tight_layout()


class DurationGantt(Plotter):
    def __init__(
        self,
        universe,
        contacts,
        lipid_type,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts
        self.lipid_type = lipid_type
        self.lipid_frequencies = None

    def get_lipid_contact_durations(self, frequency_filter=20):
        self.lipid_frequencies = get_lipid_contact_frequencies(
            self.universe, self.contacts, self.lipid_type, frequency_filter
        )

    def create_plot(self, lipid_id=None, top_filter=10, **kwargs):
        if self.lipid_frequencies is None:
            raise ValueError("Please run get_lipid_contact_durations() first.")
        elif lipid_id is None:
            raise ValueError("Please specify a lipid_id.")
        else:
            residues = self.contacts.get_residues_by_lipid_id(lipid_id=lipid_id)
            per_res_freq = {}
            for res in residues:
                per_res_freq[res] = len(self.contacts.contact_frames[res][lipid_id])
            # sort by values
            temp = sorted(per_res_freq.items(), key=lambda x: x[1], reverse=True)
            top_residues = temp[:top_filter]

            res_chunks = {}
            for res in top_residues:
                res_chunks[res[0]] = find_continuous_chunks(
                    self.contacts.contact_frames[res[0]][lipid_id]
                )

            fig, ax = plt.subplots(figsize=(10, 6))

            # Create Gantt chart
            for res, chunks in inverse_dict_keys(res_chunks).items():
                for chunk in chunks:
                    ax.barh(str(res), chunk[1] - chunk[0], left=chunk[0])

            new_patches = []
            for patch in reversed(ax.patches):
                bb = patch.get_bbox()
                color = patch.get_facecolor()
                p_bbox = FancyBboxPatch(
                    (bb.xmin, bb.ymin),
                    abs(bb.width),
                    abs(bb.height),
                    boxstyle="round,pad=-0.003,rounding_size=4",
                    # boxstyle="round",
                    ec="none",
                    fc=color,
                    mutation_aspect=0.2,
                )
                patch.remove()
                new_patches.append(p_bbox)
            for patch in new_patches:
                ax.add_patch(patch)

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "durations_gantt_{}_{}.pdf".format(self.lipid_type, lipid_id),
                )
            if self.title is None:
                self.title = "Lipid Contact Durations - {} {}".format(
                    self.lipid_type, lipid_id
                )
            if self.xlabel is None:
                self.xlabel = "Trajectory Frames"
            if self.ylabel is None:
                self.ylabel = "Residue ID"

            # Set labels and title
            ax.tick_params(axis='both', which='major', labelsize=12)
            plt.xticks(fontname = "Arial Rounded MT Bold")
            plt.yticks(fontname = "Arial Rounded MT Bold")
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            ax.set_title(self.title, fontsize=14, weight="bold", pad=20, fontfamily="Arial Rounded MT Bold")
            plt.tight_layout()


class LogoResidues(Plotter):
    def __init__(
            self,
            universe,
            metric,
            lipid_type=None,
            metric_name=None,
            xlabel=None,
            ylabel=None,
            fn=None,
            title=None,
            fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.lipid_type = lipid_type
        self.metric_name = metric_name
        self.df = create_logo_df(universe, metric, lipid=lipid_type, metric_name=metric_name)

    def create_plot(self, **kwargs):
        mat_df = lm.sequence_to_matrix(''.join(self.df['Resname'].to_list()))

        def ceildiv(number):
            if number % 75 == 0:
                return number // 75
            else:
                return number // 75 + 1
            
        n_rows = ceildiv(len(self.df['Resname']))

        fig, axs = plt.subplots(n_rows, figsize=[10,(n_rows*0.75) + 0.4])

        # create colormap based on Metric
        cmap = mpl.cm.get_cmap('Blues')
        norm = mpl.colors.Normalize(vmin=self.df['Metric'].min(), vmax=self.df['Metric'].max())
        colors = [cmap(norm(value)) for value in self.df['Metric']]
        # colors
        if n_rows == 1:
            ww_logo = lm.Logo(mat_df, ax=axs, color_scheme='silver', vpad=.4, font_name='Arial Rounded MT Bold') 
            for residx, metric in enumerate(self.df['Metric']):
                ww_logo.highlight_position(p=residx, color=list(colors[residx][:3]), alpha=1)
            ww_logo.style_spines(visible=False)

            plt.xticks(range(0, len(self.df['Resname']), 5), self.df['ResID'][::5], rotation=0)
            axs.tick_params(axis='x', which='major', labelsize=12)
            axs.set_yticks([])
            axs.grid(False)
        else:
            magic_number = math.ceil(len(self.df['Resname'])/n_rows)
            for i in range(n_rows):
                if i == n_rows - 1:
                    ww_logo = lm.Logo(mat_df[i*magic_number:], ax=axs[i], color_scheme='silver', vpad=.4, font_name='Arial Rounded MT Bold')  
                    ww_logo.ax.set_xticks(range(i*magic_number, len(self.df['Resname']), 5), self.df['ResID'][i*magic_number::5], rotation=0)
                else:
                    ww_logo = lm.Logo(mat_df[i*magic_number:(i+1)*magic_number], ax=axs[i], color_scheme='silver', vpad=.4, font_name='Arial Rounded MT Bold')  
                    ww_logo.ax.set_xticks(range(i*magic_number, (i+1)*magic_number, 5), self.df['ResID'][i*magic_number:(i+1)*magic_number:5], rotation=0)

                for residx, metric in enumerate(self.df['Metric']):
                    ww_logo.highlight_position(p=residx, color=list(colors[residx][:3]), alpha=1)
                ww_logo.style_spines(visible=False)

            for ax in axs.flat:
                ax.tick_params(axis='x', which='major', labelsize=12)
                ax.set_yticks([])
                ax.grid(False)

        if self.fn is None:
            self.fn = os.path.join(
                os.getcwd(),
                "logo_{}_{}.pdf".format(self.lipid_type, self.metric_name),
            )
        if self.title is None:
            self.title = 'Logo based on {}'.format(self.metric_name)
        if self.xlabel is None:
            self.xlabel = "Residue ID"

        # add unique xlabel
        fig.text(0.5, -0.1/(n_rows - n_rows/2), self.xlabel, ha='center', fontsize=12, fontfamily='Arial Rounded MT Bold')

        # # add color bar to the bottom
        cax = plt.axes([0.1, -0.3/(n_rows - n_rows/2), 0.8, 0.12/n_rows])
        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', ticks=[self.df['Metric'].min(), self.df['Metric'].max()/2, self.df['Metric'].max()])
        cbar.set_label(label=self.metric_name, size=12, fontfamily='Arial Rounded MT Bold')
        cbar.ax.tick_params(labelsize=12)

        fig.suptitle(self.title, fontsize=14, weight="bold", fontfamily='Arial Rounded MT Bold')
        plt.tight_layout()
