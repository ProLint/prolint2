from prolint2.metrics.metrics import create_metric


class ServerPayload:
    """Class that provides the data for the dashboard.

    Parameters
    ----------
    contacts : :class:`ContactsProvider`
        The contacts provider object.
    ts : :class:`Universe`
        The universe object.

    """

    def __init__(self, contacts, ts):
        self.contacts = contacts

        self.registry = ts.registry

        self.database_resnames = ts.database.unique_resnames.tolist()
        self.database_resname_counts = ts.database.resname_counts
        self.residue_names = ts.query.residues.resnames
        self.residue_ids = ts.query.residues.resids

        self.dt = ts.trajectory.dt
        self.totaltime = ts.trajectory.totaltime

        self.ordered_lipid_names = list(self.database_resname_counts.keys())

        self._compute()

    def residue_contacts(self, lipid_type: str = None, metric="sum", dt=1, totaltime=1):
        """Compute residue contacts."""
        metric_instance = create_metric(
            self.contacts,
            metrics=[metric],
            metric_registry=self.registry,
            output_format="dashboard",
            lipid_type=self.ordered_lipid_names[0]
            if lipid_type is None
            else lipid_type,
            residue_names=self.residue_names,
            residue_ids=self.residue_ids,
        )
        return metric_instance.compute(dt=dt, totaltime=totaltime)
        # return metric_instance.compute(dt=self.dt, totaltime=self.totaltime)

    def _compute(self, lipid_type: str = None, metric="sum"):
        # protein name is hardcoded -> read protein name(s) dynamically
        # update code to handle multiple identical proteins
        # update code to handle multiple copies of different proteins
        protein_name = (
            "Protein"  # TODO: we'll need to update this into a list and iterate over it
        )
        proteins = [protein_name]
        protein_counts = {protein_name: 1}

        residue_contacts = self.residue_contacts(lipid_type=lipid_type, metric=metric)
        # print ('residue_contacts', residue_contacts)

        lipid_counts = self.database_resname_counts
        total_lipid_sum = sum(lipid_counts.values())
        sub_data = []
        for lipid, count in lipid_counts.items():
            sub_data.append(
                {"category": lipid, "value": "{:.2f}".format(count / total_lipid_sum)}
            )

        pie_data = []
        for protein in proteins:
            value = protein_counts[protein] / sum(protein_counts.values())

            protein_pdata = {
                "category": protein_name,
                "value": "{:.2f}".format(value),
                "subData": sub_data,
            }
            pie_data.append(protein_pdata)

        self._payload = {
            "data": {protein_name: residue_contacts},
            "proteins": [protein_name],
            "lipids": self.ordered_lipid_names,
            "pie_data": pie_data,  # TODO: include protein info
        }

    @property
    def payload(self):
        """The payload."""
        return self._payload

    def get_payload(self):
        """Return the payload."""
        return self.payload
