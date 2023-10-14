r""":mod:`prolint2.plotting.projection`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os
import matplotlib as mpl
from matplotlib.pyplot import cm
from .utils import *
import nglview as nv


def show_contact_projection(
    universe,
    metric,
    metric_name=None,
    lipid=None,
    res_ids=None,
    query_repr="surface",
    database_repr="spacefill",
    cmap="Reds",
    frame_idx=None,
):
    """
    Generate a contact projection visualization based on specified parameters.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe containing the trajectory data.
    metric : str
        The metric to use for coloring. Must be one of 'max', 'sum', 'mean'.
    metric_name : str, optional
        Name of the specific metric to use (if available).
    lipid : str, optional
        The name of the target lipid for visualization.
    res_ids : list, optional
        List of residue IDs to include in the visualization.
    query_repr : str, optional
        Representation style for query structure ('surface' or other).
    database_repr : str, optional
        Representation style for the database structure ('spacefill' or other).
    cmap : str, optional
        Colormap to use for coloring.
    frame_idx : int, optional
        Index of a specific frame to visualize.

    Returns
    -------
    nglview.NGLWidget
        An NGLWidget displaying the contact projection visualization.
    nv.color._ColorScheme
        The color scheme used in the visualization.

    Notes
    -----
    This function creates a contact projection visualization that colors residues
    based on the specified metric values.

    Examples
    --------
    >>> universe = MDAnalysis.Universe("trajectory.dcd", "topology.pdb")
    >>> view, scheme = show_contact_projection(universe, 'max', lipid="DOPC")
    >>> view
    >>> scheme
    """
    # Obtain a list of metrics calculated for each residue
    res_list, metric_list = get_metric_list_by_residues(
        universe, metric, lipid=lipid, metric_name=metric_name, res_list=res_ids
    )

    all_resids = universe.query.residues.resids.tolist()

    # Create a dictionary mapping residue indices to their corresponding metrics
    metric_dict = {
        res_id: metric_list[res_list.index(all_resids[idx])]
        if all_resids[idx] in res_list
        else 0
        for idx, res_id in enumerate(universe.query.residues.resindices)
    }

    # Calculate the maximum metric value
    max_metric = max(metric_list)

    # Obtain the list of residue names in the database
    dat_resnames = universe.database.residues.resnames

    # Create a dictionary that assigns the maximum metric value to lipids and 0 to other residues
    metric_dist2 = {
        res_id: max_metric if dat_resnames[idx] == lipid else 0
        for idx, res_id in enumerate(universe.database.residues.resindices)
    }

    # Update the metric dictionary with lipid-specific metrics
    metric_dict.update(metric_dist2)

    # Generate a list of atomic B-factors using the calculated metric values
    atomic_bfactors = [metric_dict[x] for x in universe.atoms.resindices]

    # Get a colormap for B-factor coloring
    bf_cmap = cm.get_cmap(cmap)

    # Convert B-factor values to color hex codes using the colormap
    colors = [mpl.colors.to_hex(x) for x in bf_cmap(shift_range(atomic_bfactors))]

    # Create a color scheme with B-factor-based coloring
    cs = [[y, str(universe.atoms.resindices[x])] for x, y in enumerate(colors)]
    scheme = nv.color._ColorScheme(cs, "bf")

    # Generate visualization using NGLView
    if frame_idx is not None:
        # Visualize a specific frame
        universe.trajectory[frame_idx]
        universe.query.atoms.write("temp.pdb")
        universe.database.atoms.write("temp2.pdb")

        # Load structures and create the visualization
        with open("temp.pdb", "r") as myfile:
            with open("temp2.pdb") as myfile2:
                query_structure = nv.TextStructure(myfile.read())
                database_structure = nv.TextStructure(myfile2.read())
                view = nv.NGLWidget()
                view.add_component(query_structure)
                view.add_component(database_structure)
                view.clear_representations(component=0)
                view.clear_representations(component=1)

                # Add representations based on chosen styles (surface or other)
                if query_repr == "surface":
                    view.add_representation(
                        "surface",
                        surfaceType="av",
                        probeRadius=2.1,
                        component=0,
                        color=scheme,
                    )
                else:
                    view.add_representation(query_repr, component=0, color=scheme)
                view.add_representation(database_repr, component=1, color=scheme)
                view.center(component=0)

        # Remove temporary files
        os.remove("temp.pdb")
        os.remove("temp2.pdb")
    else:
        # Visualize the entire trajectory
        view = nv.NGLWidget()
        view.add_trajectory(nv.MDAnalysisTrajectory(universe.query))
        view.add_trajectory(nv.MDAnalysisTrajectory(universe.database))
        view.clear_representations(component=0)
        view.clear_representations(component=1)

        # Add representations based on chosen styles (surface or other)
        if query_repr == "surface":
            view.add_representation(
                "surface", surfaceType="av", probeRadius=2.1, component=0, color=scheme
            )
        else:
            view.add_representation(query_repr, component=0, color=scheme)
        view.add_representation(database_repr, component=1, color=scheme)
        view.center(component=0)

    # Return the visualization and color scheme
    return view, scheme
