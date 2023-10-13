import os
import ast
from io import StringIO
from collections import Counter

import MDAnalysis as mda
from bottle import Bottle, redirect, static_file

from prolint2.core.universe import Universe
from prolint2.server.chord_utils import contact_chord
from prolint2.computers.payload import ServerPayload
from prolint2.computers.distances import SerialDistances

SERVER_PATH = os.path.abspath(os.path.dirname(__file__))

import os
import configparser

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini"))
parameters_config = config["Parameters"]


class ProLintDashboard:
    """
    A dashboard for ProLint2. It is a fully functional web application that can be used to
    visualize the results of the ProLint2 analysis.
    """

    def __init__(self, port=8351, debug_bool=False, reloader=False):
        self.backend_data = None
        self.ts = None
        self.contacts = None
        self.payload = None
        self.args = None
        self.data = None
        self.data_loaded = False
        self.port = port
        self.debug_bool = debug_bool
        self.reloader = reloader
        self.response = None
        self.app = Bottle()
        self.setup_routes()

    def setup_routes(self):
        """
        Setup the routes for the dashboard.
        """
        self.app.route("/", method="GET", callback=self.redirect_to_app)
        self.app.route("/app", method="GET", callback=self.serve_app)
        self.app.route(
            "/static/<filepath:path>", method="GET", callback=self.server_static
        )
        self.app.route("/data/<metadata>", method="GET", callback=self.serve_data)
        self.app.route("/pdb/<metadata>", method="GET", callback=self.serve_pdb)
        self.app.route("/network/<metadata>", method="GET", callback=self.serve_network)
        self.app.route(
            "/tabledata/<metadata>", method="GET", callback=self.serve_table_data
        )
        self.app.route(
            "/toplipids/<metadata>", method="GET", callback=self.serve_top_lipids
        )
        self.app.route(
            "/distance/<metadata>", method="GET", callback=self.serve_distance_array
        )
        self.app.route("/metric/<metadata>", method="GET", callback=self.update_metric)

    def serve_app(self):
        """
        Serve the main application.
        """
        return static_file("index.html", root=SERVER_PATH)

    def server_static(self, filepath):
        """
        Serve static files.
        """
        return static_file(filepath, root=os.path.join(SERVER_PATH, "static"))

    def redirect_to_app(self):
        """
        Redirect to the main application.
        """
        redirect("/app")

    def serve_pdb(self, metadata):
        """
        Serve the PDB file to the client for use with Mol*.
        """
        u = mda.Universe(self.args["structure"], self.args["trajectory"])
        protein = u.select_atoms("protein")
        pstream = mda.lib.util.NamedStream(StringIO(), "dummy.pdb")
        with mda.Writer(pstream, format="PDB") as w:
            w.write(protein)

        return pstream.read()

    def get_gantt_app_data(
        self,
        g,
        lipid_id,
        residues_to_show=int(parameters_config["residues_to_show"]),
        intervals_to_filter_out=int(parameters_config["intervals_to_filter_out"]),
    ):
        """
        Get the data for the Gantt chart in the application.

        Args:
            g (dict): For each lipid ID, stores the residues in contact and the corresponding
                    frequency.
            lipid_id (str): The lipid ID to use.
            residues_to_show (int): The number of residues to show in the Gantt chart.
            intervals_to_filter_out (int): The number of frames to filter out.

        Returns:
            gantt_data (list): A list of dictionaries containing the data for the Gantt chart.
            categories (list): A list of residue IDs.
        """

        gantt_data = []
        for res, _ in g[lipid_id][:residues_to_show]:
            frame_numbers = self.contacts.contact_frames[res][lipid_id]
            frame_intervals = self.get_frame_contact_intervals(frame_numbers)
            for start, end in frame_intervals:
                if end - start < intervals_to_filter_out:
                    continue
                gantt_data.append(
                    {
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

    def sort_lipids(self):
        """
        Sort lipid contacts according to their contact frequency, all the while keeping track
        of residue IDs, and number of contacts with each residue.

        Args:
            ts (PL2): The ProLint2 object.

        Returns:
            lipid_frequency (dict): Stores lipid IDs and their contact frequency, sorted in descending order
            residue_contact_freq (dict): For each lipid ID, stores the residues in contact and the corresponding
                    frequency.

        """

        def sort_by_frequency(contact_list):
            contact_list.sort(key=lambda x: x[1], reverse=True)
            return contact_list

        # TODO:
        # top lipid number should be put in the config.
        # contact_threshold = self.ts.trajectory.n_frames * 0.05
        contact_threshold = 0

        # Initialize dictionaries to store values:
        lipid_frequency = {lipid: {} for lipid in self.ts.database.unique_resnames}
        residue_contact_freq = {}

        for residue, lipid_contacts in self.contacts.compute_metric("sum").items():
            for lipid, contact_counter in lipid_contacts.items():
                # Sort the contact_counter dictionary by its values
                sorted_contacts = sorted(
                    contact_counter.items(), key=lambda x: x[1], reverse=True
                )

                for lipid_id, freq in sorted_contacts:
                    # Exclude short-lived contacts
                    if freq <= contact_threshold:
                        continue

                    # Update lipid_frequency
                    if lipid_id in lipid_frequency[lipid]:
                        lipid_frequency[lipid][int(lipid_id)] += freq
                    else:
                        lipid_frequency[lipid][int(lipid_id)] = freq

                    # Update residue_contact_freq
                    if int(lipid_id) in residue_contact_freq:
                        residue_contact_freq[int(lipid_id)].append((int(residue), freq))
                    else:
                        residue_contact_freq[int(lipid_id)] = [(int(residue), freq)]

        for lipid, values in lipid_frequency.items():
            lipid_frequency[lipid] = Counter(values).most_common()

        for lipid_id, vals in residue_contact_freq.items():
            residue_contact_freq[lipid_id] = sort_by_frequency(vals)

        return lipid_frequency, residue_contact_freq

    @staticmethod
    def get_frame_contact_intervals(
        frames, tolerance=int(parameters_config["tolerance"])
    ):
        """
        Get the intervals of frames in which a contact is present.

        Args:
            frames (list): A list of frames in which a contact is present.
            tolerance (int): The number of frames to tolerate before considering a new interval.

        Returns:
            ranges_collect (list): A list of tuples containing the start and end frames of each
                interval.

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
        """
        Serve the data to the client for use with the application.
        """

        metadata = ast.literal_eval(metadata)

        lipid = metadata["lipid"]
        protein = metadata["protein"]
        metric = metadata.get("metric", "")

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
            self.ts.query.universe,
            self.ts.query,
            self.ts.database,
            lipid_id,
            residue_id,
            self.contacts.contact_frames[residue_id][lipid_id],
        )
        ri.run(verbose=False)

        hm_data, la_data = [], []
        for lx, la in enumerate(ri.lipid_atomnames):
            la_data.append({"LipidAtoms": la})
            for rx, ra in enumerate(ri.resid_atomnames):
                v = ri.distance_array[lx, rx]
                hm_data.append(
                    {"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)}
                )
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
            "frameNumber": self.ts.trajectory.n_frames,
        }
        self.response = response
        return response

    def update_metric(self, metadata):
        """
        Update the metric used for the backend data.
        """
        metadata = ast.literal_eval(metadata)

        lipid = metadata["lipid"]
        metric = metadata["metric"]

        residue_contacts = self.payload.residue_contacts(
            lipid_type=lipid, metric=metric
        )

        self.response["data"] = residue_contacts[lipid]
        return self.response

    def serve_network(self, metadata):
        """
        Serve the data to the client for use with the network application.
        """
        metadata = ast.literal_eval(metadata)
        lipid = metadata["lipid"]

        top_lipid_ids = [x[0] for x in self.backend_data["top_lipids"][lipid]]
        chord_elements, hidden_node_indices, per_lipid_nodes = contact_chord(
            self.ts,
            self.contacts,
            top_lipid_ids,
            self.backend_data["lipid_contact_frames"],
            cutoff=100,
        )

        return {
            "chordElements": chord_elements,
            "positionResidues": hidden_node_indices,
            "lipidNodes": per_lipid_nodes,
        }

    def serve_table_data(self, metadata):
        """
        Serve the data to the client for use with the table application.
        """
        metadata = ast.literal_eval(metadata)
        lipid = metadata["lipid"]

        table_data = []
        for ix, (lipid_id, freq) in enumerate(self.backend_data["top_lipids"][lipid]):
            table_data.append({"id": ix, "lipidID": lipid_id, "contactFrequency": freq})

        return {
            "tableData": table_data,
        }

    def serve_top_lipids(self, metadata):
        """
        Serve the data to the client for use with the gantt application.
        """
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
        """
        Serve the data to the client for use with the heatmap application.
        """
        metadata = ast.literal_eval(metadata)
        lipid_id = metadata["lipidID"]
        residue_id = int(metadata["residueID"])

        ri = SerialDistances(
            self.ts.query.universe,
            self.ts.query,
            self.ts.database,
            lipid_id,
            residue_id,
            self.contacts.contact_frames[residue_id][lipid_id],
        )
        ri.run(verbose=False)

        hm_data, la_data = [], []
        for lx, la in enumerate(ri.lipid_atomnames):
            la_data.append({"LipidAtoms": la})
            for rx, ra in enumerate(ri.resid_atomnames):
                v = ri.distance_array[lx, rx]
                hm_data.append(
                    {"LipidAtoms": la, "ResidueAtoms": ra, "value": float(v)}
                )
        ra_data = [{"ResidueAtoms": x} for x in ri.resid_atomnames]

        return {
            "heatmapData": hm_data,
            "lipidAtomsData": la_data,
            "residueAtomsData": ra_data,
        }

    def start_server(self, payload=None):
        """
        Start the server.
        """
        if payload is None:
            import sys  # pylint: disable=import-outside-toplevel

            print("Please provide a payload")
            sys.exit(1)

        self.args = payload
        self.port = self.args["port"]
        self.ts = Universe(self.args["structure"], self.args["trajectory"])
        # load contacts from file if provided
        if "contacts" in self.args.keys():
            self.contacts = self.ts.load_contacts_from_file(self.args["contacts"])
        else:
            self.contacts = self.ts.compute_contacts(
                cutoff=float(parameters_config["cutoff"])
            )

        self.payload = ServerPayload(self.contacts, self.ts)
        payload = self.payload.payload

        lipid_frequency, residue_contact_freq = self.sort_lipids()
        payload["top_lipids"] = lipid_frequency
        payload["lipid_contact_frames"] = residue_contact_freq

        self.backend_data = payload

        self.app.run(
            reloader=self.reloader,
            host=self.args["host"],
            port=self.port,
            debug=self.debug_bool,
        )


if __name__ == "__main__":
    app = ProLintDashboard(debug_bool=True)
    app.start_server()
