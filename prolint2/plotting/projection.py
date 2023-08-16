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
    ngl_repr="surface",
    cmap="Reds",
    frame_idx=None,
    show_database=False,
):
    metric_list = get_metric_list_by_residues(
        universe, metric, lipid=lipid, metric_name=metric_name
    )

    metric_dict = {
        res_id: metric_list[idx]
        for idx, res_id in enumerate(universe.query.residues.resids)
    }

    atomic_bfactors = [metric_dict[x] for x in universe.query.atoms.resids]

    bf_cmap = cm.get_cmap(cmap)

    colors = [mpl.colors.to_hex(x) for x in bf_cmap(shift_range(atomic_bfactors))]
    cs = [[y, str(universe.query.atoms.resids[x])] for x, y in enumerate(colors)]

    scheme = nv.color._ColorScheme(cs, "bf")

    reprs = [
        "point",
        "line",
        "rope",
        "tube",
        "trace",
        "label",
        "cartoon",
        "licorice",
        "ribbon",
        "backbone",
        "spacefill",
        "ball+stick",
    ]

    params = {
        "surface": [
            {
                "type": "surface",
                "params": {"surfaceType": "av", "probeRadius": 2.1, "color": scheme},
            }
        ],
    }

    for rep in reprs:
        params[rep] = [{"type": rep, "params": {"color": scheme}}]

    if frame_idx != None:
        universe.trajectory[frame_idx]
        universe.query.atoms.write("temp.pdb")

        with open("temp.pdb", "r") as myfile:
            query_structure = nv.TextStructure(myfile.read())
            view = nv.NGLWidget()
            view.add_component(query_structure)
            view.set_representations(params[ngl_repr])

        # close file and remove it
        myfile.close()
        os.remove("temp.pdb")
    else:
        view = nv.show_mdanalysis(universe.query)
        view.set_representations(params[ngl_repr])

    return view
