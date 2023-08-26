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
from prolint2.computers.distances import SerialDistances
from matplotlib.cm import ScalarMappable
import inspect
import configparser
from prolint2.server.chord_utils import contact_chord
from .utils import *

## seaborn config for paper quality plots
import seaborn as sns

sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_palette("colorblind")

__all__ = [
    "Plotter",
    "PointDistribution",
    "Radar",
    "DensityMap",
    "DurationGantt",
    "LogoResidues",
    "InteractionHeatMap",
    "RadarMetrics",
    "SharedContacts",
]

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini"))
parameters_config = config["Parameters"]


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
        plt.savefig(self.fn, dpi=300, bbox_inches="tight")

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
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def create_plot(self, lipid_type=None, metric_name=None, **kwargs):
        """Plot the distribution of a metric for each residue."""
        metric_names = []
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names:
                        metric_names.append(key)
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        elif metric_name is None or metric_name not in metric_names:
            raise ValueError(
                "Please specify a valid metric name: {}.".format(metric_names)
            )
        else:
            metric_list = get_metric_list_by_residues(
                self.universe, self.metric, lipid_type, metric_name
            )

            fig, ax = plt.subplots(figsize=self.fig_size)

            if "palette" in kwargs:
                ax = sns.scatterplot(
                    x=self.universe.query.residues.resids,
                    y=metric_list,
                    hue=metric_list,
                    **kwargs
                )
                norm = plt.Normalize(metric_list.min(), metric_list.max())
                sm = plt.cm.ScalarMappable(cmap=kwargs["palette"], norm=norm)
                sm.set_array([])

                # remove legend and add color bar
                ax.get_legend().remove()
                cbar = ax.figure.colorbar(sm)
                cbar.ax.tick_params(labelsize=12)
            else:
                ax = sns.scatterplot(
                    x=self.universe.query.residues.resids, y=metric_list, **kwargs
                )

            ax.tick_params(axis="both", which="major", labelsize=12)

            if self.xlabel is None:
                self.xlabel = "Residue ID"
            if self.ylabel is None:
                self.ylabel = metric_name
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "point_distribution_{}_{}.pdf".format(lipid_type, metric_name),
                )
            if self.title is None:
                self.title = "Point distribution - {} - {}".format(
                    lipid_type, metric_name
                )

            # add labels and title
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.title(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class Radar(Plotter):
    def __init__(
        self,
        metric,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.metric = metric

    def create_plot(self, resIDs=None, lipid_type=None, metric_names=None, **kwargs):
        metric_names_aux = []
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names_aux:
                        metric_names_aux.append(key)
        if resIDs is None:
            raise ValueError("Please specify resIDs.")
        elif lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        elif metric_names is None or not all(
            item in metric_names_aux for item in metric_names
        ):
            raise ValueError(
                "Please specify valid metric names: {}.".format(metric_names_aux)
            )
        else:
            metric_dict, metric_names = get_metrics_for_radar(
                self.metric, metric_names, resIDs=resIDs, lipid=lipid_type
            )
            num_vars = len(metric_names)
            theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)
            theta += np.pi / 2

            fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw={"polar": True})
            ax.set_xticks(theta)
            ax.set_xticklabels(metric_names)
            ax.tick_params(axis="x", pad=20)
            ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)

            for resi in metric_dict.keys():
                values = metric_dict[resi]
                xs = np.concatenate((theta, [theta[0]]))
                values = np.concatenate((values, [values[0]]))
                ax.plot(xs, values, label=resi, **kwargs)
                ax.fill(xs, values, alpha=0.1)

            ax.set_theta_offset(np.pi / 2)
            ax.set_theta_direction(-1)

            ax.tick_params(axis="both", which="major", labelsize=10)
            tt = ax.legend(
                loc="upper right",
                bbox_to_anchor=(1.15, 1.1),
                title="Residue ID",
                title_fontsize=12,
                prop={"family": "Arial Unicode MS", "size": 12},
            )
            tt.set_title(
                title="Residue ID", prop={"family": "Arial Rounded MT Bold", "size": 12}
            )

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(), "metrics_comparison_{}.pdf".format(lipid_type)
                )
            if self.title is None:
                self.title = "Metrics comparison - {}".format(lipid_type)

            # add labels and title
            plt.yticks(fontname="Arial Unicode MS")
            plt.xticks(fontname="Arial Rounded MT Bold")
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=30,
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class DensityMap(Plotter):
    def __init__(
        self,
        universe,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe

    def create_plot(self, lipid_type=None, bins=150, size_in_mb=50000, **kwargs):
        """Plot the preferential localization of lipids using 2D density maps."""
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        else:
            # Compute the lipid coordinates
            computed_coords = compute_density(self.universe, lipid_type, size_in_mb)

            # Compute the lipid density
            H, xe, ye = np.histogram2d(
                computed_coords[:, 0], computed_coords[:, 1], bins=bins, density=True
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
            cbar.set_label(
                label="Density distribution",
                size=12,
                fontfamily="Arial Rounded MT Bold",
            )
            cbar.ax.tick_params(labelsize=12)

            # Remove labels and ticks
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(), "density_map_{}.pdf".format(lipid_type)
                )
            if self.title is None:
                self.title = "Density map - {}".format(lipid_type)

            # add title
            ax.set_title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=15,
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class DurationGantt(Plotter):
    def __init__(
        self,
        universe,
        contacts,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts

    def get_contact_durations(self, lipid_type, frequency_filter=20):
        return get_lipid_contact_durations(
            self.universe, self.contacts, lipid_type, frequency_filter
        )

    def create_plot(
        self,
        lipid_id=None,
        top_filter=10,
        continuity_filter=int(parameters_config["intervals_to_filter_out"]),
        tolerance=int(parameters_config["tolerance"]),
        **kwargs
    ):
        if lipid_id is None:
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
                res_chunks[res[0]] = get_frame_contact_intervals(
                    self.contacts.contact_frames[res[0]][lipid_id],
                    continuity_filter,
                    tolerance,
                )

            fig, ax = plt.subplots(figsize=(10, 6))

            # Create Gantt chart
            for res, chunks in inverse_dict_keys(res_chunks).items():
                for chunk in chunks:
                    ax.barh(str(res), chunk[1] - chunk[0], left=chunk[0], **kwargs)

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
                    "durations_gantt_lipidID_{}.pdf".format(lipid_id),
                )
            if self.title is None:
                self.title = "Lipid Contact Durations - Lipid ID: {}".format(lipid_id)
            if self.xlabel is None:
                self.xlabel = "Trajectory Frames"
            if self.ylabel is None:
                self.ylabel = "Residue ID"

            # Set labels and title
            ax.tick_params(axis="both", which="major", labelsize=12)
            plt.xticks(fontname="Arial Rounded MT Bold")
            plt.yticks(fontname="Arial Rounded MT Bold")
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            ax.set_title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=20,
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class LogoResidues(Plotter):
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
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def create_plot(self, lipid_type=None, metric_name=None, color_logo="silver", palette="Blues", **kwargs):
        metric_names = []
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names:
                        metric_names.append(key)
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        elif metric_name is None or metric_name not in metric_names:
            raise ValueError(
                "Please specify a valid metric name: {}.".format(metric_names)
            )
        else:
            df = create_logo_df(
                self.universe, self.metric, lipid=lipid_type, metric_name=metric_name
            )
            mat_df = lm.sequence_to_matrix("".join(df["Resname"].to_list()))

            def ceildiv(number):
                if number % 75 == 0:
                    return number // 75
                else:
                    return number // 75 + 1

            n_rows = ceildiv(len(df["Resname"]))

            fig, axs = plt.subplots(n_rows, figsize=[10, (n_rows * 0.75) + 0.4])

            # create colormap based on Metric
            cmap = mpl.cm.get_cmap(palette)
            norm = mpl.colors.Normalize(
                vmin=df["Metric"].min(), vmax=df["Metric"].max()
            )
            colors = [cmap(norm(value)) for value in df["Metric"]]
            # colors
            if n_rows == 1:
                ww_logo = lm.Logo(
                    mat_df,
                    ax=axs,
                    color_scheme="silver",
                    vpad=0.4,
                    font_name="Arial Rounded MT Bold",
                )
                for residx, metric in enumerate(df["Metric"]):
                    ww_logo.highlight_position(
                        p=residx, color=list(colors[residx][:3]), alpha=1
                    )
                ww_logo.style_spines(visible=False)

                plt.xticks(
                    range(0, len(df["Resname"]), 5), df["ResID"][::5], rotation=0
                )
                axs.tick_params(axis="x", which="major", labelsize=12)
                axs.set_yticks([])
                axs.grid(False)
            else:
                magic_number = math.ceil(len(df["Resname"]) / n_rows)
                for i in range(n_rows):
                    if i == n_rows - 1:
                        ww_logo = lm.Logo(
                            mat_df[i * magic_number :],
                            ax=axs[i],
                            color_scheme="silver",
                            vpad=0.4,
                            font_name="Arial Rounded MT Bold",
                            **kwargs
                        )
                        ww_logo.ax.set_xticks(
                            range(i * magic_number, len(df["Resname"]), 5),
                            df["ResID"][i * magic_number :: 5],
                            rotation=0,
                        )
                    else:
                        ww_logo = lm.Logo(
                            mat_df[i * magic_number : (i + 1) * magic_number],
                            ax=axs[i],
                            color_scheme=color_logo,
                            vpad=0.4,
                            font_name="Arial Rounded MT Bold",
                            **kwargs
                        )
                        ww_logo.ax.set_xticks(
                            range(i * magic_number, (i + 1) * magic_number, 5),
                            df["ResID"][i * magic_number : (i + 1) * magic_number : 5],
                            rotation=0,
                        )

                    for residx, metric in enumerate(df["Metric"]):
                        ww_logo.highlight_position(
                            p=residx, color=list(colors[residx][:3]), alpha=1
                        )
                    ww_logo.style_spines(visible=False)

                for ax in axs.flat:
                    ax.tick_params(axis="x", which="major", labelsize=12)
                    ax.set_yticks([])
                    ax.grid(False)

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "logo_{}_{}.pdf".format(lipid_type, metric_name),
                )
            if self.title is None:
                self.title = "Logo based on {}".format(metric_name)
            if self.xlabel is None:
                self.xlabel = "Residue ID"

            # add unique xlabel
            fig.text(
                0.5,
                -0.1 / (n_rows - n_rows / 2),
                self.xlabel,
                ha="center",
                fontsize=12,
                fontfamily="Arial Rounded MT Bold",
            )

            # # add color bar to the bottom
            cax = plt.axes([0.1, -0.3 / (n_rows - n_rows / 2), 0.8, 0.12 / n_rows])
            cbar = plt.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cax,
                orientation="horizontal",
                ticks=[
                    df["Metric"].min(),
                    df["Metric"].min() + (df["Metric"].max() - df["Metric"].min()) / 2,
                    df["Metric"].max(),
                ],
            )
            cbar.set_label(
                label=metric_name, size=12, fontfamily="Arial Rounded MT Bold"
            )
            cbar.ax.tick_params(labelsize=12)

            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class InteractionHeatMap(Plotter):
    def __init__(
        self,
        universe,
        contacts,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts

    def create_plot(self, residue_id=None, lipid_id=None, palette='Reds', **kwargs):
        if residue_id == None:
            raise ValueError("Please specify a residue_id.")
        elif lipid_id == None:
            raise ValueError("Please specify a lipid_id.")
        else:
            ri = SerialDistances(
                self.universe,
                self.universe.query,
                self.universe.database,
                lipid_id,
                residue_id,
                self.contacts.contact_frames[residue_id][lipid_id],
            )
            ri.run(verbose=False)

            hm_data = []
            for rx, ra in enumerate(ri.resid_atomnames):
                hm_data.append([])
                for lx, la in enumerate(ri.lipid_atomnames):
                    v = ri.distance_array[lx, rx]
                    hm_data[rx].append(float(v))
            min_value = ri.distance_array.min()
            max_value = ri.distance_array.max()

            fig, ax = plt.subplots(figsize=(8, 5))

            # Define color map
            cmap = mpl.cm.get_cmap(palette)
            cmap_r = reverse_colourmap(cmap)

            heatmap = ax.pcolor(hm_data, cmap=cmap_r, **kwargs)

            # Set axis labels
            ax.set_yticks(np.arange(len(ri.resid_atomnames)) + 0.5, [])
            ax.set_yticklabels(ri.resid_atomnames)
            ax.set_xticks(np.arange(len(ri.lipid_atomnames)) + 0.5, [])
            ax.set_xticklabels(ri.lipid_atomnames)

            # Create colorbar
            cbar = plt.colorbar(heatmap)
            cbar.set_label(
                r"Distance ($\AA$)", size=12, fontfamily="Arial Rounded MT Bold"
            )
            cbar.ax.tick_params(labelsize=12)
            cbar.set_ticks(
                [min_value, min_value + (max_value - min_value) / 2, max_value]
            )
            cbar.ax.invert_yaxis()

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "int_matrix_{}_{}_{}_{}.pdf".format(
                        self.universe.residues.resnames[
                            self.universe.residues.resids.tolist().index(residue_id)
                        ],
                        residue_id,
                        self.universe.residues.resnames[
                            self.universe.residues.resids.tolist().index(lipid_id)
                        ],
                        lipid_id,
                    ),
                )
            if self.title is None:
                self.title = "Interactions matrix between {} {} and {} {}".format(
                    self.universe.residues.resnames[
                        self.universe.residues.resids.tolist().index(residue_id)
                    ],
                    residue_id,
                    self.universe.residues.resnames[
                        self.universe.residues.resids.tolist().index(lipid_id)
                    ],
                    lipid_id,
                )
            if self.xlabel is None:
                self.xlabel = "Lipid atoms"
            if self.ylabel is None:
                self.ylabel = "Residue atoms"

            ax.tick_params(axis="both", which="major", labelsize=12)
            plt.xticks(fontname="Arial Unicode MS")
            plt.yticks(fontname="Arial Unicode MS")
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")

            ax.set_title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=20,
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class RadarMetrics(Plotter):
    def __init__(
        self,
        universe,
        metric,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(10, 10),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def create_plot(self, lipid=None, metric_name=None, palette='Reds', **kwargs):
        if lipid is None:
            raise ValueError("Please specify a lipid.")
        elif metric_name is None:
            raise ValueError("Please specify a metric_name.")
        else:
            metrics = get_metric_list_by_residues(
                self.universe, self.metric, lipid, metric_name
            )

            num_vars = len(metrics)
            theta = np.linspace(0, 2 * np.pi - np.pi / 6, num_vars, endpoint=False)
            theta += np.pi / 12

            fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw={"polar": True})
            ax.set_theta_zero_location("S")
            ax.set_theta_direction(-1)
            magic_number = len(metrics) // 32
            cmap = mpl.cm.get_cmap(palette)
            rescale = lambda metrics: (metrics - np.min(metrics)) / (
                np.max(metrics) - np.min(metrics)
            )

            # Plot the radar chart
            radar_bar = ax.bar(
                theta, metrics, width=0.01, alpha=0.5, color=cmap(rescale(metrics))
            )

            # customize x-axis
            ax.set_xticks(theta[::magic_number])
            ax.set_xticklabels(
                [
                    "{} {}".format(
                        self.universe.query.residues.resnames[i],
                        self.universe.query.residues.resids[i],
                    )
                    for i in range(
                        0, len(self.universe.query.residues.resids), magic_number
                    )
                ],
                fontfamily="Arial Rounded MT Bold",
            )
            ax.tick_params(axis="x", pad=20, labelsize=12)
            ax.xaxis.grid(True, linestyle="dashed", alpha=0.5)
            plt.xlim(0, 2 * np.pi)

            # customize y-axis
            ax.set_yticks([0, max(metrics) / 2, max(metrics)])
            for loc in [0, max(metrics) / 2, max(metrics)]:
                ax.text(
                    0.0,
                    loc,
                    loc,
                    ha="center",
                    va="top",
                    fontsize=12,
                    fontfamily="Arial Rounded MT Bold",
                )
            ax.set_yticklabels([])
            ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)
            ax.set_rorigin(-0.4 * max(metrics))
            ax.set_rlabel_position(0)  # Move radial labels away from plotted line
            ax.set_thetamin(15)
            ax.set_thetamax(345)

            # customize colorbar
            sm = ScalarMappable(
                cmap=cmap, norm=plt.Normalize(min(metrics), max(metrics))
            )
            sm.set_array([])
            cax = plt.axes([0.1, -0.04, 0.8, 0.02])
            cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
            cbar.set_label(metric_name, size=12, fontfamily="Arial Rounded MT Bold")
            cbar.ax.tick_params(labelsize=12)
            cbar.set_ticks(
                [
                    min(metrics),
                    min(metrics) + (max(metrics) - min(metrics)) / 2,
                    max(metrics),
                ]
            )

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "radar_metrics_{}_{}.pdf".format(lipid, metric_name),
                )
            if self.title is None:
                self.title = "Radar Chart | {}".format(lipid)
            if self.ylabel is None:
                self.ylabel = metric_name

            ax.text(
                0.0,
                -0.4 * max(metrics),
                self.ylabel,
                ha="center",
                va="top",
                fontsize=12,
                fontfamily="Arial Rounded MT Bold",
            )
            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )

            plt.tight_layout()


class SharedContacts(Plotter):
    def __init__(
        self,
        universe,
        contacts,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts

    def create_plot(self, lipid_type=None, label_size=6, palette='Reds', **kwargs):
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        else:
            top_lipids, residue_contact_freq = get_lipid_contact_frequencies(
                self.universe, self.contacts, lipid_type
            )
            chord_elements, hidden_node_indices, per_lipid_nodes = contact_chord(
                self.universe, self.contacts, top_lipids, residue_contact_freq
            )
            df = pd.DataFrame(chord_elements)
            df = df[df["value"] != 0]
            uniques_labels = np.unique(
                df["to"].unique().tolist() + df["from"].unique().tolist()
            )
            n_nodes = len(uniques_labels) + len(hidden_node_indices)
            edges = df.to_dict("records")

            # create plot
            fig, ax = plt.subplots(figsize=self.fig_size)

            # Create a directed graph
            G = nx.DiGraph()

            # create list with ordered unique labels
            labels = uniques_labels.tolist()
            labels.sort(key=lambda x: (int(x.split(" ")[0])))

            # create list of color values for nodes
            color_values = []
            for lab in labels:
                value = 0
                for j in edges:
                    if lab == j["from"] or lab == j["to"]:
                        value += j["value"]
                color_values.append(value)

            # Add nodes with labels
            for node in range(n_nodes):
                G.add_node(node)

            pos = nx.circular_layout(G)

            # Calculate an offset angle of 90 degrees in radians
            offset_angle = math.radians(-90)

            # Apply the offset angle to each node's angle
            for node, (x, y) in pos.items():
                angle = math.atan2(-y, x)
                new_angle = angle + offset_angle
                radius = math.sqrt(x**2 + y**2)
                pos[node] = (radius * math.cos(new_angle), radius * math.sin(new_angle))

            # Draw nodes and labels
            for i in hidden_node_indices:
                G.remove_node(i)

            # Create a layout
            nx.draw_networkx_nodes(
                G, pos, node_size=10, node_color=color_values, cmap=mpl.cm.get_cmap(palette)
            )

            # adding labels to nodes
            nx.set_node_attributes(
                G, {list(G.nodes)[i]: labels[i] for i in range(len(labels))}, "label"
            )

            node_labels = nx.get_node_attributes(G, "label")
            temp_labels = nx.draw_networkx_labels(
                G, pos, labels=node_labels, font_size=label_size
            )
            theta = {k: np.arctan2(v[1], v[0]) * 180 / np.pi for k, v in pos.items()}

            for key, t in temp_labels.items():
                if 90 < theta[key] or theta[key] < -90:
                    angle = 180 + theta[key]
                    t.set_ha("right")
                    sep_pad = 0.015
                    position = (
                        t.get_position()[0]
                        - sep_pad * np.cos(angle / (360.0 / (2.0 * np.pi))),
                        t.get_position()[1]
                        - sep_pad * np.sin(angle / (360.0 / (2.0 * np.pi))),
                    )
                else:
                    angle = theta[key]
                    t.set_ha("left")
                    position = (
                        t.get_position()[0]
                        + sep_pad * np.cos(angle / (360.0 / (2.0 * np.pi))),
                        t.get_position()[1]
                        + sep_pad * np.sin(angle / (360.0 / (2.0 * np.pi))),
                    )
                t.set_position(position)
                t.set_va("center_baseline")
                t.set_rotation(angle)
                t.set_rotation_mode("anchor")
                t.set_clip_on(False)

            # Draw edges
            for node_i, label_i in dict(G.nodes(data=True)).items():
                for node_j, label_j in dict(G.nodes(data=True)).items():
                    if (label_i["label"], label_j["label"]) in [
                        (x["from"], x["to"]) for x in edges
                    ]:
                        # print('OMG')
                        if node_i > node_j:
                            if node_i - node_j > n_nodes / 2:
                                G.add_edge(node_i, node_j)
                            else:
                                G.add_edge(node_j, node_i)
                        else:
                            if node_j - node_i > n_nodes / 2:
                                G.add_edge(node_j, node_i)
                            else:
                                G.add_edge(node_i, node_j)

            # Draw edges
            for edge in G.edges():
                rad = fig.get_figheight() * 2.5 / (edge[1] - edge[0])
                if rad < 0:
                    rad = -rad
                nx.draw_networkx_edges(
                    G,
                    pos,
                    edgelist=[edge],
                    width=0.5,
                    arrowstyle="-",
                    alpha=0.5,
                    edge_color="gray",
                    connectionstyle="arc3,rad={}".format(rad),
                )

            # Add color bar
            sm = plt.cm.ScalarMappable(
                cmap=mpl.cm.get_cmap(palette),
                norm=plt.Normalize(vmin=min(color_values), vmax=max(color_values)),
            )
            sm.set_array([])
            cax = plt.axes([0.1, -0.01, 0.8, 0.02])
            cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
            cbar.set_label(
                "Total number of shared contacts",
                size=12,
                fontfamily="Arial Rounded MT Bold",
            )
            cbar.ax.tick_params(labelsize=10)
            cbar.set_ticks(
                [
                    min(color_values),
                    min(color_values) + (max(color_values) - min(color_values)) / 2,
                    max(color_values),
                ]
            )

            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "shared_contacts_{}.pdf".format(lipid_type),
                )
            if self.title is None:
                self.title = "Shared Contacts | {}".format(lipid_type)

            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )
            plt.grid(False)
            ax.axis("off")
            plt.tight_layout()
