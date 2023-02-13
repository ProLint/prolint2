r"""Interactive selection of Database and Query
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import MDAnalysis as mda
from collections import OrderedDict

# Function to select the query and the database groups interactively starting from a PL2 object. It returns a new PL2 object with the selected groups for the query and the database.
def interactive_selection(target_system):

    # default selections for query and database
    db_qr = {
        "Database": target_system.database.selected,
        "Query": target_system.query.selected,
    }

    # dictionary to store the groups for selections
    groups_dict = {
        0: ["System", target_system.query.whole.universe.atoms],
        1: ["Protein", target_system.query.whole],
        2: ["Lipids", target_system.database.whole],
    }

    # adding groups for each lipid type in the database
    groups_to_add = target_system.database.whole.universe.select_atoms(
        "not protein"
    ).groupby("resnames")
    groups_to_add = OrderedDict(sorted(groups_to_add.items()))
    for label, ag in groups_to_add.items():
        groups_dict[max(groups_dict.keys()) + 1] = ["{}".format(label), ag]

    # function to print action labels
    def print_action_keys():
        print(
            """\nActions keys:\ndb\t:  update database (i.e. >> db Group_ID).\nqr\t:  update query (i.e. >> qr Group_ID).\ngb\t:  group by topological attribute (i.e. >> gp Group_ID resnames).\nsl\t:  split group (i.e. >> sl Group_ID residue)\nadd\t:  merge two or more groups (i.e. >> add Group_ID1 Group_ID2 ... Group_IDn).\ndel\t:  delete a group (i.e. >> del Group_ID).\nlg\t:  print the groups.\nh\t:  print the list of action keys.\ne\t:  exit interactive selection mode.\nFor a more detailed description of the actions, please refer to the README.md file.
            """
        )

    # function to print the groups for selections
    def print_list_groups():
        print("\nSelection groups:\n")
        for key, value in groups_dict.items():
            print("({}) : {} --> {} atoms".format(key, value[0], value[1].n_atoms))

    # function to print the query and database selections
    def print_db_qr():
        print("\nDatabase and Query groups:\n")
        for key, value in db_qr.items():
            print("{} --> {} atoms".format(key, value.n_atoms))
        print("-----------------------------------------------------")

    # function to combine two or more AtomGroups
    def add_atomgroups(ag_list):
        combined = ag_list[0]
        for ag in ag_list[1:]:
            combined = mda.Merge(combined.atoms, ag)
        return combined

    # calling print functions
    print_action_keys()
    print_db_qr()
    print_list_groups()

    # stariting the interactive selection mode
    input_key = ""
    # loop to keep the interactive selection mode open until the user decides to exit ("e")
    while input_key != "e":
        input_key = input(">> ")
        if input_key != "e" and len(input_key) > 0:
            # splitting the input key into a list of strings
            input_list = input_key.split()

            # updating database action
            if input_list[0] == "db":
                if len(input_list) > 2:
                    print("Error: Too many arguments.")
                elif len(input_list) < 2:
                    print("Error: Too few arguments.")
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print("Group_ID must be an integer.")
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print("Error: Group_ID must be a valid group ID.")
                    else:
                        target_system.database.select(
                            groups_dict[int(input_list[1])][1]
                        )
                        db_qr["Database"] = groups_dict[int(input_list[1])][1]
                        print(
                            "Database group updated to {} atoms.".format(
                                groups_dict[int(input_list[1])][1].n_atoms
                            )
                        )

            # updating query action
            elif input_list[0] == "qr":
                if len(input_list) > 2:
                    print("Error: Too many arguments.")
                elif len(input_list) < 2:
                    print("Error: Too few arguments.")
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print("Group_ID must be an integer.")
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print("Error: Group_ID must be a valid group ID.")
                    else:
                        target_system.query.select(groups_dict[int(input_list[1])][1])
                        db_qr["Query"] = groups_dict[int(input_list[1])][1]
                        print(
                            "Query group updated to {} atoms.".format(
                                groups_dict[int(input_list[1])][1].n_atoms
                            )
                        )

            # grouping by topological attribute action
            elif input_list[0] == "gb":
                if len(input_list) > 3:
                    print("Error: Too many arguments.")
                elif len(input_list) < 3:
                    print("Error: Too few arguments.")
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print("Group_ID must be an integer.")
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print("Error: Group_ID must be a valid group ID.")
                    elif not hasattr(groups_dict[int(input_list[1])][1], input_list[2]):
                        print(
                            "The attribute {} is not included in the topological attributes of this system.".format(
                                input_list[2]
                            )
                        )
                    else:
                        groups_to_add = groups_dict[int(input_list[1])][1].groupby(
                            input_list[2]
                        )
                        groups_to_add = OrderedDict(sorted(groups_to_add.items()))
                        for label, ag in groups_to_add.items():
                            groups_dict[max(groups_dict.keys()) + 1] = [
                                "{} grouped by {} from {}".format(
                                    label,
                                    input_list[2],
                                    groups_dict[int(input_list[1])][0],
                                ),
                                ag,
                            ]

            # splitting group action
            elif input_list[0] == "sl":
                if len(input_list) > 3:
                    print("Error: Too many arguments.")
                elif len(input_list) < 3:
                    print("Error: Too few arguments.")
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print("Group_ID must be an integer.")
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print("Error: Group_ID must be a valid group ID.")
                    elif input_list[2] not in [
                        "segment",
                        "residue",
                        "atom",
                    ]:
                        print(
                            "Error: Level must be a valid level: segment, residue, atom."
                        )
                    else:
                        groups_to_add = groups_dict[int(input_list[1])][1].split(
                            input_list[2]
                        )
                        if input_list[2] == "atom":
                            groups_to_add = sorted(
                                groups_to_add, key=lambda x: x.atoms.ids[0]
                            )
                            for ag in groups_to_add:
                                groups_dict[max(groups_dict.keys()) + 1] = [
                                    "{} {} splitted from {}".format(
                                        input_list[2],
                                        ag.atoms.ids[0],
                                        groups_dict[int(input_list[1])][0],
                                    ),
                                    ag,
                                ]
                        elif input_list[2] == "residue":
                            groups_to_add = sorted(
                                groups_to_add, key=lambda x: x.resids[0]
                            )
                            for ag in groups_to_add:
                                groups_dict[max(groups_dict.keys()) + 1] = [
                                    "{} {} splitted from {}".format(
                                        input_list[2],
                                        ag.resids[0],
                                        groups_dict[int(input_list[1])][0],
                                    ),
                                    ag,
                                ]
                        elif input_list[2] == "segment":
                            groups_to_add = sorted(
                                groups_to_add, key=lambda x: x.segids[0]
                            )
                            for ag in groups_to_add:
                                groups_dict[max(groups_dict.keys()) + 1] = [
                                    "{} {} splitted from {}".format(
                                        input_list[2],
                                        ag.segids[0],
                                        groups_dict[int(input_list[1])][0],
                                    ),
                                    ag,
                                ]

            # combining groups action
            elif input_list[0] == "add":
                if len(input_list) < 3:
                    print("Error: Too few arguments.")
                else:
                    try:
                        for i in range(1, len(input_list)):
                            int(input_list[i])
                    except ValueError:
                        print("All Group_IDs must be integers.")
                        continue
                    if not all(
                        [int(item) in groups_dict.keys() for item in input_list[1:]]
                    ):
                        print("Error: All Group_IDs must be valid group IDs.")
                    else:
                        ag_list = [groups_dict[int(x)][1] for x in input_list[1:]]
                        combined = add_atomgroups(ag_list).atoms
                        combined.dimensions = groups_dict[int(input_list[2])][
                            1
                        ].dimensions
                        groups_dict[max(groups_dict.keys()) + 1] = [
                            "Group combined from {}".format(input_list[1:]),
                            combined,
                        ]

            # removing group action
            elif input_list[0] == "del":
                if len(input_list) > 2:
                    print("Error: Too many arguments.")
                elif len(input_list) < 2:
                    print("Error: Too few arguments.")
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print("Group_ID must be an integer.")
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print("Error: Group_ID must be a valid group ID.")
                    else:
                        groups_dict.pop(int(input_list[1]))

            # help actions
            elif input_list[0] == "lg":
                if len(input_list) > 1:
                    print("Error: Too many arguments.")
                else:
                    print_db_qr()
                    print_list_groups()
            elif input_list[0] == "h":
                if len(input_list) > 1:
                    print("Error: Too many arguments.")
                else:
                    print_action_keys()

            # unknown actions
            else:
                print(
                    "Unknown action key. Please type (h) to see the available action keys."
                )

        # if no action key is given keep the loop going
        elif input_key != "e" and len(input_key) == 0:
            continue

        # exiting the interactive selection mode
        else:
            print("Exiting interactive selection.")
            break
    return target_system
