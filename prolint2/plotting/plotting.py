r""":mod:`prolint2.plotting.plotting`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

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
from prolint2.computers.distances import TwoPointDistances
from matplotlib.cm import ScalarMappable
import inspect
from scipy import interpolate
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
    "MetricsComparison",
    "SharedContacts",
    "TwoPointDistanceEvolution",
    "MosaicsGridData",
]


class Plotter:
    """
    Initialize the Plotter instance with provided parameters.

    Parameters
    ----------
    xlabel : str, optional
        X-axis label for the plot.
    ylabel : str, optional
        Y-axis label for the plot.
    fn : str, optional
        File name to save the plot as an image.
    title : str, optional
        Title of the plot.
    fig_size : tuple, optional
        Figure size (width, height) for the plot.

    """

    def __init__(
        self,
        xlabel: str = None,
        ylabel: str = None,
        fn: str = None,
        title: str = None,
        fig_size: tuple = (8, 8),
    ):
        self.xlabel = xlabel  # X-axis label for the plot
        self.ylabel = ylabel  # Y-axis label for the plot
        self.fn = fn  # File name to save the plot as an image
        self.title = title  # Title of the plot
        self.fig_size = fig_size  # Figure size (width, height) for the plot

    def save_plot(self, **kwargs):
        """
        Generates the plot using create_plot method and saves it to the specified file.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments to pass to create_plot.

        """
        # Generates the plot using create_plot method and saves it to the specified file
        self.create_plot(**kwargs)
        plt.savefig(self.fn, dpi=300, bbox_inches="tight")

    def generate_script(self, class_code, script_filename):
        """
        Generates a Python script that generates the plot using provided class code.

        Parameters
        ----------
        class_code : object
            The class code to generate the script for.
        script_filename : str
            The name of the script file to be generated.

        """

        # Extract the source code of the class_code passed as an argument
        plotting_function_source = inspect.getsource(class_code)

        # Construct different tail sections based on the class name for specific use cases
        if self.__class__.__name__ in ["PointDistribution"]:
            tail = point_distribution_tail(
                self.__class__.__name__
            )  # Tail script for PointDistribution case

        elif self.__class__.__name__ in ["MetricsComparison"]:
            tail = radar_tail(self.__class__.__name__)  # Tail script for Radar case

        elif self.__class__.__name__ in ["DensityMap"]:
            tail = density_map_tail(
                self.__class__.__name__
            )  # Tail script for DensityMap case

        elif self.__class__.__name__ in ["DurationGantt"]:
            tail = duration_gantt_tail(
                self.__class__.__name__,
            )  # Tail script for DurationGantt case

        elif self.__class__.__name__ in ["LogoResidues"]:
            tail = logo_tail(
                self.__class__.__name__
            )  # Tail script for LogoResidues case

        elif self.__class__.__name__ in ["InteractionHeatMap"]:
            tail = interaction_map_tail(
                self.__class__.__name__
            )  # Tail script for InteractionHeatMap case

        elif self.__class__.__name__ in ["Radar"]:
            tail = radar_metrics_tail(
                self.__class__.__name__
            )  # Tail script for RadarMetrics case

        elif self.__class__.__name__ in ["SharedContacts"]:
            tail = shared_contacts_tail(
                self.__class__.__name__
            )  # Tail script for SharedContacts case

        elif self.__class__.__name__ in ["TwoPointDistanceEvolution"]:
            tail = two_point_distance_evolution_tail(
                self.__class__.__name__
            )  # Tail script for TwoPointDistanceEvolution case
        elif self.__class__.__name__ in ["MosaicsGridData"]:
            tail = mosaics_grid_data_tail(
                self.__class__.__name__
            )  # Tail script for MosaicsGridData case

        # Combine class code and the corresponding tail script
        if self.__class__.__name__ in ["MosaicsGridData"]:
            script_code = use_mosaics_script_template(plotting_function_source, tail)
        else:
            script_code = use_1d_script_template(plotting_function_source, tail)

        # Write the generated script code to the specified file
        with open(script_filename, "w") as f:
            f.write(script_code)


class PointDistribution(Plotter):
    """
    Initialize the PointDistribution instance.

    Parameters
    ----------
    universe : Universe
        The simulation universe containing the data.
    metric : dict
        A dictionary containing metric data for residues and lipids.
    xlabel : str, optional
        The label for the x-axis of the plot.
    ylabel : str, optional
        The label for the y-axis of the plot.
    fn : str, optional
        The filename to save the plot. If None, a default filename is generated.
    title : str, optional
        The title of the plot.
    fig_size : tuple, optional
        A tuple specifying the figure size (width, height).

    Returns
    -------
    PointDistribution
        An instance of the PointDistribution class.

    Examples
    --------
    >>> pd = PointDistribution(universe, metric)
    >>> pd = PointDistribution(universe, metric, xlabel="X Label", ylabel="Y Label")
    """

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

    def create_plot(self, res_ids=None, lipid_type=None, metric_name=None, y_lim=None, **kwargs):
        """
        Plot the distribution of a metric for each residue.

        Parameters
        ----------
        res_ids : list of int, optional
            The list of residue IDs to include in the plot. If None, include all residues.
        lipid_type : str
            The lipid type to focus on in the plot.
        metric_name : str
            The name of the metric to plot.
        **kwargs
            Additional keyword arguments to customize the plot.

        Raises
        ------
        ValueError
            If lipid_type is None or metric_name is not found in the metric data.

        Examples
        --------
        >>> pd.create_plot(lipid_type="DOPC", metric_name="metric1")
        >>> pd.create_plot(lipid_type="POPC", metric_name="metric2", color="blue")
        """
        metric_names = []
        # Collect all unique metric names from the metric dictionary
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
            # Get the list of metric values for specified lipid_type and metric_name
            res_list, metric_list = get_metric_list_by_residues(
                self.universe, self.metric, lipid_type, metric_name, res_list=res_ids
            )

            # Create a scatter plot
            fig, ax = plt.subplots(figsize=self.fig_size)

            if "palette" in kwargs:
                # If palette is specified, use it for coloring scatter points
                ax = sns.scatterplot(
                    x=res_list, y=metric_list, hue=metric_list, **kwargs
                )
                norm = plt.Normalize(metric_list.min(), metric_list.max())
                sm = plt.cm.ScalarMappable(cmap=kwargs["palette"], norm=norm)
                sm.set_array([])

                # Remove legend and add color bar
                ax.get_legend().remove()
                cbar = ax.figure.colorbar(sm)
                cbar.ax.tick_params(labelsize=12)
            else:
                # Create scatter plot without palette
                ax = sns.scatterplot(x=res_list, y=metric_list, **kwargs)

            ax.tick_params(axis="both", which="major", labelsize=12)

            if y_lim is not None:
                ax.set_ylim(y_lim)

            # Set default labels, filename, and title if not provided
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

            # Add labels and title to the plot
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.title(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )

            plt.tight_layout()


class MetricsComparison(Plotter):
    """
    Initialize the MetricsComparison object with specified plot attributes.

    Parameters
    ----------
    metric : dict
        A dictionary containing metric data.
    xlabel : str, optional
        The label for the x-axis.
    ylabel : str, optional
        The label for the y-axis.
    fn : str, optional
        The filename for saving the plot.
    title : str, optional
        The title of the plot.
    fig_size : tuple, optional
        The size of the figure (width, height).

    Attributes
    ----------
    metric : dict
        A dictionary containing metric data.

    Example
    -------
    >>> metric_data = {...}
    >>> plot = MetricsComparison(metric_data, xlabel="X Label", ylabel="Y Label")
    >>> plot.create_plot()
    """

    def __init__(
        self,
        metric,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        # Initialize the MetricsComparison object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.metric = metric  # Store the input metric data

    def create_plot(self, resIDs=None, lipid_type=None, metric_names=None, **kwargs):
        """
        Create a radar plot comparing metrics for different residues.

        Parameters
        ----------
        resIDs : list, optional
            List of residue IDs to include in the plot.
        lipid_type : str, optional
            The type of lipid to focus on.
        metric_names : list, optional
            List of metric names to display on the plot.

        Raises
        ------
        ValueError
            If required input parameters are not provided or metric_names is invalid.

        Example
        -------
        >>> plot = MetricsComparison(metric_data)
        >>> plot.create_plot(resIDs=[1, 2, 3], lipid_type="DOPC", metric_names=["metric1", "metric2"])
        """
        metric_names_aux = []
        # Extract all unique metric names from the metric data
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names_aux:
                        metric_names_aux.append(key)

        # Check if the necessary inputs are provided
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
            # Retrieve relevant metric data for creating the radar plot
            metric_dict, metric_names = get_metrics_for_radar(
                self.metric, metric_names, resIDs=resIDs, lipid=lipid_type
            )

            # Calculate the number of variables (metrics) for the radar plot
            num_vars = len(metric_names)
            theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)
            theta += np.pi / 2

            # Create a polar plot
            fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw={"polar": True})
            ax.set_xticks(theta)
            ax.set_xticklabels(metric_names)
            ax.tick_params(axis="x", pad=20)
            ax.yaxis.grid(True, linestyle="dashed", alpha=0.5)

            # Plot radar lines and fill the areas
            for resi in metric_dict.keys():
                values = metric_dict[resi]
                xs = np.concatenate((theta, [theta[0]]))
                values = np.concatenate((values, [values[0]]))
                ax.plot(xs, values, label=resi, **kwargs)
                ax.fill(xs, values, alpha=0.1)

            ax.set_theta_offset(np.pi / 2)
            ax.set_theta_direction(-1)

            # Customize tick labels and legend
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

            # Set default values for filename and title if not provided
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(), "metrics_comparison_{}.pdf".format(lipid_type)
                )
            if self.title is None:
                self.title = "Metrics comparison - {}".format(lipid_type)

            # Add labels and title to the plot
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
    """
    Initialize the DensityMap object with specified plot attributes.

    Parameters
    ----------
    universe : object
        The universe object containing the data for the density map.
    xlabel : str, optional
        The label for the x-axis. Defaults to None.
    ylabel : str, optional
        The label for the y-axis. Defaults to None.
    fn : str, optional
        The filename for the saved plot. Defaults to None.
    title : str, optional
        The title of the plot. Defaults to None.
    fig_size : tuple, optional
        The size of the figure (width, height) in inches. Defaults to (8, 8).
    """

    def __init__(
        self,
        universe,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 8),
    ):
        # Initialize the DensityMap object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe

    def create_plot(
        self, lipid_type=None, bins=150, size_in_mb=50000, start=0, stop=-1, step=1, **kwargs
    ):
        """
        Plot the preferential localization of lipids using 2D density maps.

        Parameters
        ----------
        lipid_type : str
            The type of lipid to create the density map for.
        bins : int, optional
            The number of bins for the histogram. Defaults to 150.
        size_in_mb : int, optional
            The size of the data (in megabytes) used for computing density. Defaults to 50000.
        frame : int, optional
            The number of pixels to exclude from the edges of the density map. Defaults to None.
        **kwargs
            Additional keyword arguments for customizing the plot.

        Raises
        ------
        ValueError
            If `lipid_type` is not specified.

        Returns
        -------
        None
        """
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        else:
            # Compute the lipid coordinates
            computed_coords = compute_density(self.universe, lipid_type, size_in_mb, start, stop, step)

            # Compute the lipid density using a 2D histogram
            H, xe, ye = np.histogram2d(
                computed_coords[:, 0], computed_coords[:, 1], bins=bins, density=True
            )

            # Generate a figure and axes for the plot
            fig, ax = plt.subplots(figsize=self.fig_size)

            # Plot the density map using imshow
            im = ax.imshow(
                    H, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]], **kwargs
                )
            ax.grid(False)

            # Create a colorbar of the same size as the plot
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            # Set colorbar label and font size
            cbar.set_label(
                label="Density distribution",
                size=12,
                fontfamily="Arial Rounded MT Bold",
            )
            cbar.ax.tick_params(labelsize=12)

            # Remove labels and ticks from the axes
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            # Set default filename and title if not provided
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(), "density_map_{}.pdf".format(lipid_type)
                )
            if self.title is None:
                self.title = "Density map - {}".format(lipid_type)

            # Add title to the plot
            ax.set_title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=15,
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class DurationGantt(Plotter):
    """
    Initialize the DurationGantt object with specified plot attributes.

    Parameters
    ----------
    universe : Universe
        The molecular dynamics simulation universe.
    contacts : Contacts
        The contacts data.
    xlabel : str, optional
        The label for the x-axis. Default is None.
    ylabel : str, optional
        The label for the y-axis. Default is None.
    fn : str, optional
        The filename for saving the plot. Default is None.
    title : str, optional
        The title of the plot. Default is None.
    fig_size : tuple, optional
        The figure size (width, height). Default is (8, 8).
    """

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
        # Initialize the DurationGantt object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts

    def get_contact_durations(self, lipid_type, **kwargs):
        """
        Call the external function to get lipid contact durations.

        Parameters
        ----------
        lipid_type : str
            The type of lipid to get contact durations for.
        **kwargs
            Additional keyword arguments to pass to the external function.

        Returns
        -------
        dict
            A dictionary of contact durations for the specified lipid type.
        """
        return get_lipid_contact_durations(
            self.universe, self.contacts, lipid_type, **kwargs
        )

    def create_plot(
        self, lipid_id=None, top_filter=10, continuity_filter=10, tolerance=6, **kwargs
    ):
        """
        Create a Gantt chart to visualize lipid contact durations.

        Parameters
        ----------
        lipid_id : str, optional
            The ID of the lipid to create the chart for.
        top_filter : int, optional
            The number of top residues to display. Default is 10.
        continuity_filter : int, optional
            The continuity filter for contact intervals. Default is 10.
        tolerance : int, optional
            The tolerance for contact intervals. Default is 6.
        **kwargs
            Additional keyword arguments to customize the appearance of the Gantt chart.

        Raises
        ------
        ValueError
            If lipid_id is not specified.

        """
        if lipid_id is None:
            raise ValueError("Please specify a lipid_id.")
        else:
            # Get residues associated with the given lipid_id
            residues = self.contacts.get_residues_by_lipid_id(lipid_id=lipid_id)
            per_res_freq = {}

            # Count the number of contacts for each residue
            for res in residues:
                per_res_freq[res] = len(self.contacts.contact_frames[res][lipid_id])

            # Sort residues by contact frequency
            temp = sorted(per_res_freq.items(), key=lambda x: x[1], reverse=True)
            top_residues = temp[:top_filter]

            res_chunks = {}
            # Calculate contact intervals for top residues
            for res in top_residues:
                res_chunks[res[0]] = get_frame_contact_intervals(
                    self.contacts.contact_frames[res[0]][lipid_id],
                    continuity_filter,
                    tolerance,
                )

            # Create the Gantt chart figure and axis
            fig, ax = plt.subplots(figsize=(10, 6))

            # Plot Gantt chart for each residue's contact intervals
            for res, chunks in inverse_dict_keys(res_chunks).items():
                for chunk in chunks:
                    ax.barh(str(res), chunk[1] - chunk[0], left=chunk[0], **kwargs)

            new_patches = []
            # Customize the appearance of the bars in the Gantt chart
            for patch in reversed(ax.patches):
                bb = patch.get_bbox()
                color = patch.get_facecolor()
                p_bbox = FancyBboxPatch(
                    (bb.xmin, bb.ymin),
                    abs(bb.width),
                    abs(bb.height),
                    boxstyle="round,pad=-0.003,rounding_size=4",
                    ec="none",
                    fc=color,
                    mutation_aspect=0.2,
                )
                patch.remove()
                new_patches.append(p_bbox)
            for patch in new_patches:
                ax.add_patch(patch)

            # Set default values for plot properties if not specified
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
    """
    Initialize the LogoResidues object with specified plot attributes.

    Parameters
    ----------
    universe : Universe
        The molecular dynamics universe or system to analyze.
    metric : dict
        A dictionary containing metric data for the universe.
    xlabel : str, optional
        The label for the x-axis of the plot. Default is None.
    ylabel : str, optional
        The label for the y-axis of the plot. Default is None.
    fn : str, optional
        The filename for saving the plot. Default is None.
    title : str, optional
        The title for the plot. Default is None.
    fig_size : tuple, optional
        The size of the plot figure (width, height). Default is (8, 8).
    """

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
        # Initialize the LogoResidues object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def create_plot(
        self,
        lipid_type=None,
        metric_name=None,
        color_logo="silver",
        palette="Blues",
        res_ids=None,
        limits=None,
        **kwargs
    ):
        """
        Create a Logo plot for a specific lipid and metric.

        Parameters
        ----------
        lipid_type : str, optional
            The type of lipid to analyze.
        metric_name : str, optional
            The name of the metric to visualize.
        color_logo : str, optional
            The color scheme for the Logo plot. Default is 'silver'.
        palette : str, optional
            The color palette for the plot. Default is 'Blues'.
        res_ids : list, optional
            A list of residue IDs to include in the plot. Default is None.
        **kwargs
            Additional keyword arguments for customizing the Logo plot.

        Raises
        ------
        ValueError
            If lipid_type is not specified or if an invalid metric_name is provided.

        Returns
        -------
        None
        """
        # Extract unique metric names
        metric_names = []
        for res in self.metric:
            for lip in self.metric[res]:
                for key in self.metric[res][lip]:
                    if key not in metric_names:
                        metric_names.append(key)

        # Check input values and raise exceptions if necessary
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        elif metric_name is None or metric_name not in metric_names:
            raise ValueError(
                "Please specify a valid metric name: {}.".format(metric_names)
            )
        else:
            # Create DataFrame and matrix from the universe and metric data
            df = create_logo_df(
                self.universe,
                self.metric,
                lipid=lipid_type,
                metric_name=metric_name,
                res_list=res_ids,
            )
            mat_df = lm.sequence_to_matrix("".join(df["Resname"].to_list()))

            # Function to calculate number of rows for the subplot grid
            def ceildiv(number):
                if number % 75 == 0:
                    return number // 75
                else:
                    return number // 75 + 1

            n_rows = ceildiv(len(df["Resname"]))

            # Create a subplot grid
            fig, axs = plt.subplots(n_rows, figsize=[10, (n_rows * 0.75) + 0.4])

            # Create colormap based on Metric and normalize values
            cmap = mpl.cm.get_cmap(palette)
            if limits is not None:
                norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
            else:
                norm = mpl.colors.Normalize(
                    vmin=df["Metric"].min(), vmax=df["Metric"].max()
                )
            colors = [cmap(norm(value)) for value in df["Metric"]]

            # Process the subplot grid based on the number of rows
            if n_rows == 1:
                # Create a Logo plot
                ww_logo = lm.Logo(
                    mat_df,
                    ax=axs,
                    color_scheme=color_logo,
                    vpad=0.4,
                    font_name="Arial Rounded MT Bold",
                )
                # Highlight positions based on the color map
                for residx, metric in enumerate(df["Metric"]):
                    ww_logo.highlight_position(
                        p=residx, color=list(colors[residx][:3]), alpha=1
                    )
                ww_logo.style_spines(visible=False)
                # Set x-axis ticks and labels
                plt.xticks(
                    range(0, len(df["Resname"]), 5), df["ResID"][::5], rotation=0
                )
                axs.tick_params(axis="x", which="major", labelsize=12)
                axs.set_yticks([])
                axs.grid(False)
            else:
                # Determine the number of rows and create Logo plots accordingly
                magic_number = math.ceil(len(df["Resname"]) / n_rows)
                for i in range(n_rows):
                    if i == n_rows - 1:
                        # Create a Logo plot for the last row
                        ww_logo = lm.Logo(
                            mat_df[i * magic_number :],
                            ax=axs[i],
                            color_scheme=color_logo,
                            vpad=0.4,
                            font_name="Arial Rounded MT Bold",
                            **kwargs,
                        )
                        ww_logo.ax.set_xticks(
                            range(i * magic_number, len(df["Resname"]), 5),
                            df["ResID"][i * magic_number :: 5],
                            rotation=0,
                        )
                    else:
                        # Create a Logo plot for non-last rows
                        ww_logo = lm.Logo(
                            mat_df[i * magic_number : (i + 1) * magic_number],
                            ax=axs[i],
                            color_scheme=color_logo,
                            vpad=0.4,
                            font_name="Arial Rounded MT Bold",
                            **kwargs,
                        )
                        ww_logo.ax.set_xticks(
                            range(i * magic_number, (i + 1) * magic_number, 5),
                            df["ResID"][i * magic_number : (i + 1) * magic_number : 5],
                            rotation=0,
                        )

                    # Highlight positions based on the color map
                    for residx, metric in enumerate(df["Metric"]):
                        ww_logo.highlight_position(
                            p=residx, color=list(colors[residx][:3]), alpha=1
                        )
                    ww_logo.style_spines(visible=False)

                for ax in axs.flat:
                    ax.tick_params(axis="x", which="major", labelsize=12)
                    ax.set_yticks([])
                    ax.grid(False)

            # Set default values for file name, title, and x-axis label if not provided
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "logo_{}_{}.pdf".format(lipid_type, metric_name),
                )
            if self.title is None:
                self.title = "Logo based on {}".format(metric_name)
            if self.xlabel is None:
                self.xlabel = "Residue ID"

            # Add x-axis label
            fig.text(
                0.5,
                -0.1 / (n_rows - n_rows / 2),
                self.xlabel,
                ha="center",
                fontsize=12,
                fontfamily="Arial Rounded MT Bold",
            )

            # Add color bar to the bottom
            cax = plt.axes([0.1, -0.3 / (n_rows - n_rows / 2), 0.8, 0.12 / n_rows])
            if limits is not None:
                cb_ticks=[limits[0], limits[0] + (limits[1] - limits[0]) / 2, limits[1]]
            else:
                cb_ticks=[
                    df["Metric"].min(),
                    df["Metric"].min() + (df["Metric"].max() - df["Metric"].min()) / 2,
                    df["Metric"].max(),
                ]
            cbar = plt.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cax,
                orientation="horizontal",
                ticks=cb_ticks,
            )
            cbar.set_label(
                label=metric_name, size=12, fontfamily="Arial Rounded MT Bold"
            )
            cbar.ax.tick_params(labelsize=12)

            # Set figure title and adjust layout
            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )
            plt.tight_layout()


class InteractionHeatMap(Plotter):
    """
    Initialize the InteractionHeatMap object with specified plot attributes.

    Parameters
    ----------
    universe : Universe
        The universe containing molecular dynamics data.
    contacts : Contacts
        The contacts information for the interactions.
    xlabel : str, optional
        The label for the x-axis.
    ylabel : str, optional
        The label for the y-axis.
    fn : str, optional
        The filename to save the plot to.
    title : str, optional
        The title of the plot.
    fig_size : tuple, optional
        The figure size (width, height).

    Returns
    -------
    None
    """

    def __init__(
        self,
        universe,
        contacts,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(8, 5),
    ):
        # Initialize the InteractionHeatMap object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        # Store universe and contacts information
        self.universe = universe
        self.contacts = contacts

    def create_plot(self, residue_id=None, lipid_id=None, palette="Reds", **kwargs):
        """
        Create an interaction heatmap for a specific residue and lipid pair.

        Parameters
        ----------
        residue_id : int
            The ID of the residue.
        lipid_id : int
            The ID of the lipid.
        palette : str, optional
            The color palette for the heatmap.
        **kwargs
            Additional keyword arguments for matplotlib.

        Returns
        -------
        None
        """
        # Check if both residue_id and lipid_id are specified
        if residue_id is None:
            raise ValueError("Please specify a residue_id.")
        elif lipid_id is None:
            raise ValueError("Please specify a lipid_id.")
        else:
            # Calculate distances using SerialDistances
            ri = SerialDistances(
                self.universe,
                self.universe.query,
                self.universe.database,
                lipid_id,
                residue_id,
                self.contacts.contact_frames[residue_id][lipid_id],
            )
            ri.run(verbose=False)

            # Initialize heatmap data
            hm_data = []
            for rx, ra in enumerate(ri.resid_atomnames):
                hm_data.append([])
                for lx, la in enumerate(ri.lipid_atomnames):
                    v = ri.distance_array[lx, rx]
                    hm_data[rx].append(float(v))
            min_value = ri.distance_array.min()
            max_value = ri.distance_array.max()

            # Create a new figure and axis
            fig, ax = plt.subplots(figsize=self.fig_size)

            # Define color map
            cmap = mpl.cm.get_cmap(palette)
            cmap_r = reverse_colourmap(cmap)

            # Create the heatmap
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

            # Handle default values for fn, title, xlabel, and ylabel
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

            # Customize plot appearance
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

            # Ensure plot layout
            plt.tight_layout()


class TwoPointDistanceEvolution(Plotter):
    """
    Initialize the TwoPointDistanceEvolution object.

    Parameters
    ----------
    universe : Universe
        The universe object containing simulation data.
    xlabel : str, optional
        The label for the x-axis of the plot.
    ylabel : str, optional
        The label for the y-axis of the plot.
    fn : str, optional
        The filename for saving the plot.
    title : str, optional
        The title for the plot.
    fig_size : tuple, optional
        The size of the figure (width, height).

    Notes
    -----
    This class inherits from Plotter and sets various plot attributes for creating distance evolution plots.

    """

    def __init__(
        self,
        universe,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(7, 5),
    ):
        # Initialize the TwoPointDistances object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        # Store universe and contacts information
        self.universe = universe

    def create_plot(
        self,
        lipid_id=None,
        residue_id=None,
        lipid_atomname=None,
        residue_atomname=None,
        unit="frame",
        smooth_line=False,
        n_points=250,
        useOffset=True,
        **kwargs
    ):
        """
        Create a distance evolution plot.

        Parameters
        ----------
        lipid_id : int
            The ID of the lipid of interest.
        residue_id : int
            The ID of the residue of interest.
        lipid_atomname : str, optional
            The name of the lipid atom to consider.
        residue_atomname : str, optional
            The name of the residue atom to consider.
        unit : str, optional
            The unit of the x-axis, either 'frame' or 'time'.
        smooth_line : bool, optional
            Whether to smooth the line in the plot.
        n_points : int, optional
            The number of points in the smoothed line.
        useOffset : bool, optional
            Whether to use offset notation for axis values.
        **kwargs
            Additional keyword arguments for the line plot.

        Raises
        ------
        ValueError
            If both residue_id and lipid_id are not specified.

        Notes
        -----
        This method calculates distances using the TwoPointDistances class and creates a distance evolution plot.

        """
        # check if both residue_id and lipid_id are specified
        if residue_id is None:
            raise ValueError("Please specify a residue_id.")
        elif lipid_id is None:
            raise ValueError("Please specify a lipid_id.")
        else:
            # calculate distances using TwoPointDistances
            tpd = TwoPointDistances(
                self.universe,
                self.universe.query,
                self.universe.database,
                lipid_id,
                residue_id,
                lipid_sel=lipid_atomname,
                residue_sel=residue_atomname,
                unit=unit,
            )
            tpd.run(verbose=False)

            # create a new figure and axis
            fig, ax = plt.subplots(figsize=self.fig_size)

            # create the line plot
            if smooth_line:
                # x_new, bspline, y_new
                x_new = np.linspace(
                    tpd.time_array.min(), tpd.time_array.max(), n_points
                )
                bspline = interpolate.make_interp_spline(
                    tpd.time_array, tpd.result_array
                )
                y_new = bspline(x_new)
                sns.lineplot(x=x_new, y=y_new, ax=ax, **kwargs)
            else:
                sns.lineplot(x=tpd.time_array, y=tpd.result_array, ax=ax, **kwargs)

            # Handle default values for fn, title, xlabel, and ylabel
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "distances_lipidID_{}_residueID_{}.pdf".format(
                        lipid_id,
                        residue_id,
                    ),
                )
            if self.title is None:
                self.title = "Distances between Lipid ID: {} and Residue ID: {}".format(
                    lipid_id,
                    residue_id,
                )
            if self.xlabel is None:
                if unit == "frame":
                    self.xlabel = "Trajectory Frames"
                elif unit == "time":
                    self.xlabel = "Trajectory Time (ps)"
            if self.ylabel is None:
                self.ylabel = r"Distance ($\AA$)"

            # Customize plot appearance
            ax.tick_params(axis="both", which="major", labelsize=12)
            plt.xticks(fontname="Arial Unicode MS")
            plt.yticks(fontname="Arial Unicode MS")
            plt.xlabel(self.xlabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ylabel(self.ylabel, fontsize=12, fontfamily="Arial Rounded MT Bold")
            plt.ticklabel_format(axis="both", useOffset=useOffset)
            ax.set_title(
                self.title,
                fontsize=14,
                weight="bold",
                pad=15,
                fontfamily="Arial Rounded MT Bold",
            )

            # Ensure plot layout
            plt.tight_layout()


class Radar(Plotter):
    """
    Initialize the Radar object with specified plot attributes.

    Parameters
    ----------
    universe : Universe
        The universe containing the data to be plotted.
    metric : str
        The metric to be displayed on the radar chart.
    xlabel : str, optional
        The label for the x-axis. Default is None.
    ylabel : str, optional
        The label for the y-axis. Default is None.
    fn : str, optional
        The filename for saving the radar chart. Default is None.
    title : str, optional
        The title for the radar chart. Default is None.
    fig_size : tuple, optional
        The size of the radar chart figure in inches. Default is (10, 10).
    """

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
        # Initialize the Radar object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.metric = metric

    def create_plot(
        self, res_ids=None, lipid=None, metric_name=None, palette="Reds", **kwargs
    ):
        """
        Create a radar chart plot.

        Parameters
        ----------
        res_ids : list of int, optional
            The list of residue IDs to include in the radar chart. Default is None.
        lipid : str
            The name of the lipid for which the metric is calculated.
        metric_name : str
            The name of the metric to be displayed on the radar chart.
        palette : str, optional
            The color palette to be used for the radar chart. Default is "Reds".
        **kwargs : dict, optional
            Additional keyword arguments.

        Raises
        ------
        ValueError
            If `lipid` or `metric_name` is not specified.

        Returns
        -------
        None
        """
        # Check for required arguments
        if lipid is None:
            raise ValueError("Please specify a lipid.")
        elif metric_name is None:
            raise ValueError("Please specify a metric_name.")
        else:
            # Get a list of metrics for the given lipid and metric_name
            resids, metrics = get_metric_list_by_residues(
                self.universe, self.metric, lipid, metric_name, res_list=res_ids
            )

            # Calculate the number of variables and the corresponding angles for the radar chart
            num_vars = len(metrics)
            theta = np.linspace(0, 2 * np.pi - np.pi / 6, num_vars, endpoint=False)
            theta += np.pi / 12

            # Create the radar chart figure and axis
            fig, ax = plt.subplots(figsize=self.fig_size, subplot_kw={"polar": True})
            ax.set_theta_zero_location("S")
            ax.set_theta_direction(-1)
            magic_number = len(metrics) // 32 + 1
            cmap = mpl.cm.get_cmap(palette)
            rescale = lambda metrics: (metrics - np.min(metrics)) / (
                np.max(metrics) - np.min(metrics)
            )

            # Plot the radar chart bars
            radar_bar = ax.bar(
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
            ax.set_xticks(theta[::magic_number])
            ax.set_xticklabels(
                [
                    "{} {}".format(
                        resnames[i],
                        resids[i],
                    )
                    for i in range(0, len(resids), magic_number)
                ],
                fontfamily="Arial Rounded MT Bold",
            )
            ax.tick_params(axis="x", pad=20, labelsize=12)
            ax.xaxis.grid(True, linestyle="dashed", alpha=0.5)
            plt.xlim(0, 2 * np.pi)

            # Customize y-axis
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

            # Customize colorbar
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

            # Set default values if not provided
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "radar_metrics_{}_{}.pdf".format(lipid, metric_name),
                )
            if self.title is None:
                self.title = "Radar Chart | {}".format(lipid)
            if self.ylabel is None:
                self.ylabel = metric_name

            # Add y-axis label
            ax.text(
                0.0,
                -0.4 * max(metrics),
                self.ylabel,
                ha="center",
                va="top",
                fontsize=12,
                fontfamily="Arial Rounded MT Bold",
            )
            # Add main title
            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )

            plt.tight_layout()  # Adjust layout to prevent overlapping


class SharedContacts(Plotter):
    """
    Initialize the SharedContacts object with specified plot attributes.

    Parameters
    ----------
    universe : str
        The universe object representing the simulation data.
    contacts : list
        A list of contact data.
    xlabel : str, optional
        The label for the x-axis of the plot. Default is None.
    ylabel : str, optional
        The label for the y-axis of the plot. Default is None.
    fn : str, optional
        The filename for saving the plot. Default is None.
    title : str, optional
        The title of the plot. Default is None.
    fig_size : tuple, optional
        The size of the figure (width, height) in inches. Default is (8, 8).
    """

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
        # Initialize the SharedContacts object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        self.universe = universe
        self.contacts = contacts

    def create_plot(self, lipid_type=None, label_size=6, palette="Reds", **kwargs):
        """
        Create a chord diagram plot for shared contacts.

        Parameters
        ----------
        lipid_type : str, optional
            The type of lipid for which the shared contacts will be visualized.
        label_size : int, optional
            The font size for node labels. Default is 6.
        palette : str, optional
            The color palette to use for node coloring. Default is "Reds".
        **kwargs
            Additional keyword arguments for customization.

        Raises
        ------
        ValueError
            If `lipid_type` is not specified.

        Returns
        -------
        None
        """
        # Ensure a lipid_type is provided
        if lipid_type is None:
            raise ValueError("Please specify a lipid_type.")
        else:
            # Get top lipid contact frequencies and related data
            top_lipids, residue_contact_freq = get_lipid_contact_frequencies(
                self.universe, self.contacts, lipid_type
            )

            # Calculate chord diagram elements and related data
            chord_elements, hidden_node_indices, per_lipid_nodes = contact_chord(
                self.universe, self.contacts, top_lipids, residue_contact_freq
            )

            # Create DataFrame from chord elements and filter non-zero values
            df = pd.DataFrame(chord_elements)
            df = df[df["value"] != 0]

            # Extract unique labels from data
            unique_labels = np.unique(
                df["to"].unique().tolist() + df["from"].unique().tolist()
            )

            # Calculate total number of nodes
            n_nodes = len(unique_labels) + len(hidden_node_indices)

            # Extract edges data from DataFrame
            edges = df.to_dict("records")

            # Create plot figure and axis
            fig, ax = plt.subplots(figsize=self.fig_size)

            # Create a directed graph
            G = nx.DiGraph()

            # Sort unique labels
            labels = unique_labels.tolist()
            labels.sort(key=lambda x: (int(x.split(" ")[0])))

            # Calculate color values for nodes
            color_values = []
            for lab in labels:
                value = 0
                for j in edges:
                    if lab == j["from"] or lab == j["to"]:
                        value += j["value"]
                color_values.append(value)

            # Add nodes with labels to the graph
            for node in range(n_nodes):
                G.add_node(node)

            # Create a circular layout for the graph
            pos = nx.circular_layout(G)

            # Calculate an offset angle of 90 degrees in radians
            offset_angle = math.radians(-90)

            # Apply offset angle to each node's position
            for node, (x, y) in pos.items():
                angle = math.atan2(-y, x)
                new_angle = angle + offset_angle
                radius = math.sqrt(x**2 + y**2)
                pos[node] = (radius * math.cos(new_angle), radius * math.sin(new_angle))

            # Draw nodes and labels while considering angles for alignment
            for i in hidden_node_indices:
                G.remove_node(i)

            # Create a layout
            nx.draw_networkx_nodes(
                G,
                pos,
                node_size=10,
                node_color=color_values,
                cmap=mpl.cm.get_cmap(palette),
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

            # Draw edges and apply radial arcs
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

            # Add color bar to the plot
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

            # Set default filename and title if not provided
            if self.fn is None:
                self.fn = os.path.join(
                    os.getcwd(),
                    "shared_contacts_{}.pdf".format(lipid_type),
                )

            # Add title, adjust layout, and finalize the plot
            if self.title is None:
                self.title = "Shared Contacts - {}".format(lipid_type)

            fig.suptitle(
                self.title,
                fontsize=14,
                weight="bold",
                fontfamily="Arial Rounded MT Bold",
            )

            plt.grid(False)
            ax.axis("off")
            plt.tight_layout()


class MosaicsGridData(Plotter):
    """
    Initialize the MosaicsGridData object with specified plot attributes.

    Parameters
    ----------
    grid_data : str
        The path to the grid data file.
    xlabel : str, optional
        The label for the x-axis. Default is None.
    ylabel : str, optional
        The label for the y-axis. Default is None.
    fn : str, optional
        The filename for the plot. Default is None.
    title : str, optional
        The title of the plot. Default is None.
    fig_size : tuple, optional
        The size of the figure (width, height). Default is (10, 10).

    Returns
    -------
    None

    Examples
    --------
    >>> mosaics_grid = MosaicsGridData("grid_data.txt", "X Label", "Y Label")
    >>> mosaics_grid.create_plot()
    >>> mosaics_grid.create_plot(frame=2, prop_label="Property", nan_to_zero=True)
    """

    def __init__(
        self,
        grid_data,
        xlabel=None,
        ylabel=None,
        fn=None,
        title=None,
        fig_size=(10, 10),
    ):
        # Initialize the MosaicsGridData object with specified plot attributes
        # Inherits from Plotter and sets plot labels, title, and figure size.
        super().__init__(xlabel, ylabel, fn, title, fig_size)
        # load grid data from file
        self.grid_data = np.loadtxt(grid_data)

    def create_plot(self, frame=None, prop_label=None, nan_to_zero=False, **kwargs):
        """
        Create a heatmap plot using the grid data.

        Parameters
        ----------
        frame : int, optional
            The frame to display in the plot. Default is None.
        prop_label : str, optional
            The label for the colorbar. Default is "Property".
        nan_to_zero : bool, optional
            Convert NaN values to zero if True. Default is False.
        **kwargs
            Additional keyword arguments to pass to imshow.

        Returns
        -------
        None

        Examples
        --------
        >>> mosaics_grid.create_plot(frame=1, prop_label="Density", nan_to_zero=True)
        >>> mosaics_grid.create_plot(prop_label="Temperature", cmap="hot")
        """
        if prop_label is None:
            prop_label = "Property"
        if nan_to_zero:
            self.grid_data = np.nan_to_num(self.grid_data)

        # Generate a figure and axis for the plot
        fig, ax = plt.subplots(figsize=self.fig_size)

        # Create heatmap using the grid data
        if frame is None:
            im = ax.imshow(self.grid_data, origin="lower", **kwargs)
        else:
            im = ax.imshow(
                self.grid_data[frame:-frame, frame:-frame], origin="lower", **kwargs
            )

        ax.grid(False)

        # Create a colorbar of the same size as the plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        # Set colorbar label and font size
        cbar.set_label(
            label=prop_label,
            size=12,
            fontfamily="Arial Rounded MT Bold",
        )
        cbar.ax.tick_params(labelsize=12)

        # Remove labels and ticks from the axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        # Set default filename and title if not provided
        if self.fn is None:
            self.fn = os.path.join(
                os.getcwd(), "mosaics_grid_data_{}.pdf".format(prop_label)
            )
        if self.title is None:
            self.title = "Mosaics Grid Data - {}".format(prop_label)

        # Add title to the plot
        ax.set_title(
            self.title,
            fontsize=14,
            weight="bold",
            pad=15,
            fontfamily="Arial Rounded MT Bold",
        )
        plt.tight_layout()
