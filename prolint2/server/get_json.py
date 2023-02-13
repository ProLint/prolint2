#!/usr/bin/env python

import csv
import json

# First data setup:
# For 1 protein systems
# protein_data = {Lipid: [
#   Residue: [FrameValues]
# ]}

# TODO:
# For multiple protein systems:
# {Protein: protein_data}

# TODO:
# Generalize to any system composition
# with system-agnostic terminology:
# {Reference: {InteractionObject: ReferenceUnit: [FrameValue]}}

csv_in = open("out_girk.csv")

js = {}
for row in csv.DictReader(csv_in):
    lipid = row["Lipids"]
    protein = row["Protein"]
    # if lipid != "CHOL": continue
    residue_id = row["ResName"] + " " + row["ResID"]
    lipid_number_value = float(row["Lipid_Number"])
    value = [residue_id, float("{:.2f}".format(lipid_number_value))]

    if js.get(protein):
        if js.get(protein).get(lipid):
            js[protein][lipid].append(value)
        else:
            js[protein][lipid] = [value]
    else:
        js[protein] = {lipid: [value]}

with open("girk.json", "w") as fp:
    json.dump(js, fp)
