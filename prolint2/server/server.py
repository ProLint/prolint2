from collections import Counter
from prolint2.interactive_sel import interactive_selection
import os
import ast
import json

from bottle import Bottle, redirect, route, run, template, debug, static_file, request
from prolint2.contacts import SerialDistances

import MDAnalysis as mda
from prolint2.prolint2 import PL2
from io import StringIO

from .chord_utils import contact_chord

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))

@route("/prolint2")
def prolint2():
    import sys

    print(request.body.getvalue().decode("utf-8"), file=sys.stdout)
    return request.body

class ProLintDashboard:
    """
    A dashboard for ProLint2. It is a fully functional web application that can be used to
    visualize the results of the ProLint2 analysis.
    """
    def __init__(self, port=8351, debug_bool=False, reloader=False):
        self.backend_data = None
        self.ts = None
        self.args = None
        self.data = None
        self.data_loaded = False
        self.port = port
        self.debug_bool = debug_bool
        self.reloader = reloader
        self.app = Bottle()
        self.setup_routes()

    def setup_routes(self):
        self.app.route("/", method="GET", callback=self.redirect_to_app)
        self.app.route("/app", method="GET", callback=self.serve_app)
        self.app.route("/static/<filepath:path>", method="GET", callback=self.server_static)
        self.app.route("/data/<metadata>", method="GET", callback=self.serve_data)
        self.app.route("/pdb/<metadata>", method="GET", callback=self.serve_pdb)
        self.app.route("/network/<metadata>", method="GET", callback=self.serve_network)
        self.app.route("/tabledata/<metadata>", method="GET", callback=self.serve_table_data)
        self.app.route("/toplipids/<metadata>", method="GET", callback=self.serve_top_lipids)
        self.app.route("/distance/<metadata>", method="GET", callback=self.serve_distance_array)

    def serve_app(self):
        return static_file("index.html", root=SERVER_PATH)

    def server_static(self, filepath):
        print ('server path', SERVER_PATH)
        return static_file(filepath, root=os.path.join(SERVER_PATH, "static"))

    def redirect_to_app(self):
        redirect("/app")

    def serve_pdb(self, metadata):
        u = mda.Universe(self.args.structure, self.args.trajectory)
        protein = u.select_atoms("protein")
        pstream = mda.lib.util.NamedStream(StringIO(), "dummy.pdb")
        with mda.Writer(pstream, format="PDB") as w:
            w.write(protein)

        return pstream.read()

    def get_gantt_app_data(self, g, lipid_id, residues_to_show=15, intervals_to_filter_out=10):
        gantt_data = []
        for res, _ in g[lipid_id][:residues_to_show]:
            frame_numbers = self.ts.contacts.contact_frames[f"{res},{lipid_id}"]
            frame_intervals = self.get_frame_contact_intervals(frame_numbers)
            for start, end in frame_intervals:
                if end - start < intervals_to_filter_out:
                    continue
                gantt_data.append(
                    {
                        # "category": f'{res}',
                        "category": res,
                        "startFrame": start,
                        "endFrame": end,
                        "lipid_id": lipid_id,
                    }
                )

        # TODO:
        # `categories` is now just the `gantt_data` keys.
        # replace with: `list(gantt_data.keys())` or remove entirely
        categories = []
        for y in [x["category"] for x in gantt_data]:
            if y not in categories:
                categories.append(y)
        return gantt_data, categories

    @staticmethod
    def sort_lipids(ts):
        """
        Sort lipid contacts according to their contact frequency, all the while keeping track
        of residue IDs, and number of contacts with each residue.

        Returns:
            t (dict): Stores lipid IDs and their contact frequency, sorted in descending order
            g (dict): For each lipid ID, stores the residues in contact and the corresponding
                    frequency.
        """

        def sort_tuple(tup):
            tup.sort(key=lambda x: x[1], reverse=True)
            return tup

        # TODO:
        # top lipid number should be put in the config.
        contact_threshold = ts.n_frames * 0.05

        # initialize dictionary to store values:
        t = {k: {} for k in ts.database_unique}
        g = {}
        for ix, (residue, lipid_contacts) in enumerate(ts.contacts.contacts.items()):
            for lipid, contact_counter in lipid_contacts.items():
                top10_counter = contact_counter.most_common()
                for (lipid_id, lipid_counter) in top10_counter:
                    # Exclude short-lived contacts
                    if lipid_counter <= contact_threshold:
                        continue
                    if lipid_id in t[lipid]:
                        t[lipid][lipid_id] += lipid_counter
                        g[lipid_id].append((residue, lipid_counter))
                    else:
                        t[lipid][lipid_id] = lipid_counter
                        g[lipid_id] = [(residue, lipid_counter)]

        for lipid, values in t.items():
            t[lipid] = Counter(values).most_common()

        # for lipid, values in g.items():
        for lipid_id, vals in g.items():
            g[lipid_id] = sort_tuple(vals)

        return t, g

    @staticmethod
    def get_frame_contact_intervals(frames, tolerance=6):
        """
        Get frame ranges
        """
        ranges_collect = []
        range_start = 0
        for ix, el in enumerate(frames):
            if ix == 0:
                range_start = el
                continue

            prev_el = frames[ix - 1]
            if not el - tolerance <= prev_el:
                ranges_collect.append((range_start, prev_el))
                range_start = el
            if ix == len(frames) - 1:
                ranges_collect.append((range_start, el))
        return ranges_collect

    def serve_data(self, metadata):
        metadata = ast.literal_eval(metadata)

        lipid = metadata["lipid"]
        protein = metadata["protein"]
        metric = metadata.get('metric', '')

        print(lipid, protein, metric)
        if lipid == "" and protein == "":
            # Starting setup:
            lipid = self.backend_data["lipids"][0]
            protein = self.backend_data["proteins"][0]

        table_data = []
        for ix, (lipid_id, freq) in enumerate(self.backend_data["top_lipids"][lipid]):
            table_data.append({"id": ix, "lipidID": lipid_id, "contactFrequency": freq})

        # Initiate ganttApp with the top lipid data
        lipid_id = self.backend_data["top_lipids"][lipid][0][0]
        gantt_data, categories = self.get_gantt_app_data(
            self.backend_data["lipid_contact_frames"], lipid_id
        )

        # Initiate heatmapApp with the top residue
        residue_id = self.backend_data["lipid_contact_frames"][lipid_id][0][0]
        ri = SerialDistances(
            self.ts.query.selected.universe,
            self.ts.query.selected,
            self.ts.database.selected,
            lipid_id,
            residue_id,
            self.ts.contacts.contact_frames[f"{residue_id},{lipid_id}"],
        )
        ri.run(verbose=False)

        hm_data, la_data = [], []
        for lx, la in enumerate(ri.lipid_atomnames):
            la_data.append({"LipidAtoms": la})
            for rx, ra in enumerate(ri.resid_atomnames):
                v = ri.distance_array[lx, rx]
                hm_data.append({"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)})
        ra_data = [{"ResidueAtoms": x} for x in ri.resid_atomnames]

        # TODO:
        # Possibly, avoid single point of failure on these dictionary lookups?
        response = {
            "data": self.backend_data["data"][protein][lipid],
            "proteins": self.backend_data["proteins"],
            "lipids": self.backend_data["lipids"],
            "pieData": self.backend_data["pie_data"],
            "ganttData": gantt_data,
            "topLipids": categories,
            "globalTopLipids": self.backend_data["top_lipids"],
            "lipidContactFrames": self.backend_data["lipid_contact_frames"],
            "tableData": table_data,
            "heatmapData": hm_data,
            "lipidAtomsData": la_data,
            "residueAtomsData": ra_data,
            "frameNumber": self.ts.n_frames,
        }
        return response

    def serve_network(self, metadata):
        metadata = ast.literal_eval(metadata)
        lipid = metadata['lipid']

        top_lipid_ids = [x[0] for x in self.backend_data['top_lipids'][lipid]]
        chord_elements, hidden_node_indices, per_lipid_nodes = contact_chord(
            self.ts,
            top_lipid_ids,
            self.backend_data['lipid_contact_frames'],
            cutoff=100
            )

        return {
            "chordElements": chord_elements,
            "positionResidues": hidden_node_indices,
            "lipidNodes": per_lipid_nodes
        }

    def serve_table_data(self, metadata):
        metadata = ast.literal_eval(metadata)
        lipid = metadata["lipid"]

        table_data = []
        for ix, (lipid_id, freq) in enumerate(self.backend_data["top_lipids"][lipid]):
            table_data.append({"id": ix, "lipidID": lipid_id, "contactFrequency": freq})

        return {
            "tableData": table_data,
        }

    def serve_top_lipids(self, metadata):
        metadata = ast.literal_eval(metadata)
        lipid_id = metadata["lipidID"]

        gantt_data, categories = self.get_gantt_app_data(
            self.backend_data["lipid_contact_frames"], lipid_id
        )

        return {
            "ganttData": gantt_data,
            "topLipids": categories,
        }

    def serve_distance_array(self, metadata):
        metadata = ast.literal_eval(metadata)
        lipid_id = metadata["lipidID"]
        residue_id = int(metadata["residueID"])

        ri = SerialDistances(
            self.ts.query.selected.universe,
            self.ts.query.selected,
            self.ts.database.selected,
            lipid_id,
            residue_id,
            self.ts.contacts.contact_frames[f"{residue_id},{lipid_id}"],
        )
        ri.run(verbose=False)

        hm_data, la_data = [], []
        for lx, la in enumerate(ri.lipid_atomnames):
            la_data.append({"LipidAtoms": la})
            for rx, ra in enumerate(ri.resid_atomnames):
                v = ri.distance_array[lx, rx]
                hm_data.append({"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)})
        ra_data = [{"ResidueAtoms": x} for x in ri.resid_atomnames]

        return {
            "heatmapData": hm_data,
            "lipidAtomsData": la_data,
            "residueAtomsData": ra_data,
        }


    def start_server(self, payload=None):
        self.args = payload
        self.ts = PL2(self.args.structure, self.args.trajectory, add_lipid_types=self.args.other_lipids)
        self.ts.contacts.compute(cutoff=self.args.cutoff)

        payload = self.ts.contacts.server_payload()

        t, g = self.sort_lipids(self.ts)
        payload["top_lipids"] = t
        payload["lipid_contact_frames"] = g

        self.backend_data = payload

        self.app.run(reloader=self.reloader, host="localhost", port=self.port, debug=self.debug_bool)

if __name__ == "__main__":
    # start_server(debug_bool=True)
    app = ProLintDashboard(debug_bool=True)
    app.start_server()
