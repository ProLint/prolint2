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
    query_repr="surface",
    database_repr="spacefill",
    cmap="Reds",
    frame_idx=None,
):
    metric_list = get_metric_list_by_residues(
        universe, metric, lipid=lipid, metric_name=metric_name
    )

    metric_dict = {
        res_id: metric_list[idx]
        for idx, res_id in enumerate(universe.query.residues.resindices)
    }

    max_metric = max(metric_list)
    dat_resnames = universe.database.residues.resnames
    metric_dist2 = {
        res_id: max_metric if dat_resnames[idx] == lipid else 0
        for idx, res_id in enumerate(universe.database.residues.resindices)
    }

    metric_dict.update(metric_dist2)

    atomic_bfactors = [metric_dict[x] for x in universe.atoms.resindices]

    bf_cmap = cm.get_cmap(cmap)

    colors = [mpl.colors.to_hex(x) for x in bf_cmap(shift_range(atomic_bfactors))]
    cs = [[y, str(universe.atoms.resindices[x])] for x, y in enumerate(colors)]

    scheme = nv.color._ColorScheme(cs, "bf")

    if frame_idx is not None:
        universe.trajectory[frame_idx]

        universe.query.atoms.write("temp.pdb")
        universe.database.atoms.write("temp2.pdb")

        with open("temp.pdb", "r") as myfile:
            with open("temp2.pdb") as myfile2:
                query_structure = nv.TextStructure(myfile.read())
                database_structure = nv.TextStructure(myfile2.read())
                view = nv.NGLWidget()
                view.add_component(query_structure)
                view.add_component(database_structure)
                view.clear_representations(component=0)
                view.clear_representations(component=1)
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

        os.remove("temp.pdb")
        os.remove("temp2.pdb")

    else:
        view = nv.NGLWidget()
        view.add_trajectory(nv.MDAnalysisTrajectory(universe.query))
        view.add_trajectory(nv.MDAnalysisTrajectory(universe.database))
        view.clear_representations(component=0)
        view.clear_representations(component=1)
        if query_repr == "surface":
            view.add_representation(
                "surface", surfaceType="av", probeRadius=2.1, component=0, color=scheme
            )
        else:
            view.add_representation(query_repr, component=0, color=scheme)
        view.add_representation(database_repr, component=1, color=scheme)
        view.center(component=0)

    return view, scheme
