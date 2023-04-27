from prolint2.metrics.metrics import create_metric

class ServerPayload:
    def __init__(self, contacts, ts):
        self.contacts = contacts

        self.registry = ts.registry
        self.database = ts.database
        self.dt = ts.dt
        self.totaltime = ts.totaltime
        self.residue_names = ts.query.selected.residues.resnames
        self.residue_ids = ts.query.selected.residues.resids

        self._compute()

    def residue_contacts(self, lipid_type: str = None, metric="max"):
        metric_instance = create_metric(
            self.contacts,
            metrics=[metric],
            metric_registry=self.registry,
            output_format="dashboard",
            lipid_type=self.database.lipid_types().tolist()[0] if lipid_type is None else lipid_type,
            residue_names=self.residue_names,
            residue_ids=self.residue_ids,
        )
        return metric_instance.compute(dt=self.dt, totaltime=self.totaltime)

    def _compute(self, lipid_type: str = None, metric="max"):
        # protein name is hardcoded -> read protein name(s) dynamically
        # update code to handle multiple identical proteins
        # update code to handle multiple copies of different proteins
        protein_name = "Protein"  # TODO: we'll need to update this into a list and iterate over it
        proteins = [protein_name]
        protein_counts = {protein_name: 1}

        residue_contacts = self.residue_contacts(lipid_type=lipid_type, metric=metric)

        lipid_counts = self.database.lipid_count()
        total_lipid_sum = sum(lipid_counts.values())
        sub_data = []
        for lipid, count in lipid_counts.items():
            sub_data.append({"category": lipid, "value": "{:.2f}".format(count / total_lipid_sum)})

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
            "lipids": self.database.lipid_types().tolist(),
            "pie_data": pie_data,  # TODO: include protein info
        }

    @property
    def payload(self):
        return self._payload
    
    def get_payload(self):
        return self.payload
