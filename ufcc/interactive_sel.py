r"""Interactive selection of Database and Query
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import MDAnalysis as mda


def interactive_selection(target_system):
    db_qr = {
       'Database' : target_system.database.selected,
        'Query' : target_system.query.selected
    }

    groups_dict = {
        0: ['System', target_system.query.whole.universe.atoms],
        1: ['Protein', target_system.query.whole],
        2: ['Lipids', target_system.database.whole]
    }

    groups_to_add = target_system.database.whole.universe.select_atoms('not protein').groupby('resnames')
    for label, ag in groups_to_add.items():
        groups_dict[max(groups_dict.keys()) + 1] = ['{}'.format(label), ag]

    def print_action_keys():
        print('''\nActions keys:
db  :  update database group for contacts calculation. (i.e. >> d Group_ID) where Group_ID is the identifier of the group on the list of groups.
qr  :  update query group for contacts calculation. (i.e. >> d Group_ID) where Group_ID is the identifier of the group on the list of groups.
gb :  create subgroups grouped by any topology attribute ("names", "resnames", "masses", etc). (i.e. >> gp Group_ID resnames) creates new groups for every different resname in the Group_ID selected.
sl  :  create new groups by splitting a previous group at the desired level ("segment", "residue", "atom", "molecule"). (i.e. >> sl Group_ID residue) splits the Group_ID selected into all the residues that make it up, and create a new group for each of them.
add:  merge two or more groups and create a new group with the combination. (i.e. >> add Group_ID1 Group_ID2 ... Group_IDn) creates a new group for the combination of the groups Group_ID1 Group_ID2 ... Group_IDn.
del:  delete a group from the list of groups. (i.e. >> del Group_ID) remove the Group_ID selected.
lg :  print the list of groups.
h  :  print the help for all the available action keys.
e  :  exit interactive selection and calculate de contacts between the Query and Database groups. 
''')

    def print_list_groups():
        print('\nSelection groups:\n')
        for key, value in groups_dict.items():
            print('({}) : {} --> {} atoms'.format(key, value[0], value[1].n_atoms))

    def print_db_qr():
        print('\nDatabase and Query groups:\n')
        for key, value in db_qr.items():
            print('{} --> {} atoms'.format(key, value.n_atoms))
        print('-----------------------------------------------------')
        
    def print_help_message():
        print('You can type (lg) to list the groups at any time, or (h) to review the available action keys with examples.')

    def add_atomgroups(ag_list):
        combined = ag_list[0]
        for ag in ag_list[1:]:
            combined = mda.Merge(combined, ag)
        return combined

    print_db_qr()
    print_list_groups()
    print_action_keys()
    
    input_key = ''
    while input_key != 'e':
        print_help_message()
        input_key = input('>> ')
        if input_key != 'e' and len(input_key) > 0:
            input_list = input_key.split()

            # updating database and query groups
            if input_list[0] == 'db':
                if len(input_list) > 2:
                    print('Error: Too many arguments.')
                elif len(input_list) < 2:
                    print('Error: Too few arguments.')
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print('Group_ID must be an integer.')
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print('Error: Group_ID must be a valid group ID.')
                    else:
                        target_system.database.select(groups_dict[int(input_list[1])][1])
                        db_qr['Database'] = groups_dict[int(input_list[1])][1]
            elif input_list[0] == 'qr':
                if len(input_list) > 2:
                    print('Error: Too many arguments.')
                elif len(input_list) < 2:
                    print('Error: Too few arguments.')
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print('Group_ID must be an integer.')
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print('Error: Group_ID must be a valid group ID.')
                    else:
                        target_system.query.select(groups_dict[int(input_list[1])][1])
                        db_qr['Query'] = groups_dict[int(input_list[1])][1]

            # creating subgroups
            elif input_list[0] == 'gb':
                if len(input_list) > 3:
                    print('Error: Too many arguments.')
                elif len(input_list) < 3:
                    print('Error: Too few arguments.')
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print('Group_ID must be an integer.')
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print('Error: Group_ID must be a valid group ID.')
                    elif not hasattr(groups_dict[int(input_list[1])][1], input_list[2]):
                        print('The attribute {} is not included in the topological attributes of this system.'.format(input_list[2]))
                    else:
                        groups_to_add = groups_dict[int(input_list[1])][1].groupby(input_list[2])
                        for label, ag in groups_to_add.items():
                            groups_dict[max(groups_dict.keys()) + 1] = ['{} grouped by {} from {}'.format(label, input_list[2], groups_dict[int(input_list[1])][0]), ag]
            elif input_list[0] == 'sl':
                if len(input_list) > 3:
                    print('Error: Too many arguments.')
                elif len(input_list) < 3:
                    print('Error: Too few arguments.')
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print('Group_ID must be an integer.')
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print('Error: Group_ID must be a valid group ID.')
                    elif input_list[2] not in ['segment', 'residue', 'atom', 'molecule']:
                        print('Error: Level must be a valid level: segment, residue, atom, molecule.')
                    else:
                        groups_to_add = groups_dict[int(input_list[1])][1].split(input_list[2])
                        for ag in groups_to_add:
                            groups_dict[max(groups_dict.keys()) + 1] = ['{} splitted from {}'.format(input_list[2],groups_dict[int(input_list[1])][0]), ag]

            # logical operations with groups
            elif input_list[0] == 'add':
                if len(input_list) < 3:
                    print('Error: Too few arguments.')
                else:
                    try:
                        for i in range(1, len(input_list)):
                            int(input_list[i])
                    except ValueError:
                        print('All Group_IDs must be integers.')
                        continue
                    if not all([int(item) in groups_dict.keys() for item in input_list[1:]]):  
                        print('Error: All Group_IDs must be valid group IDs.')
                    else:
                        ag_list = [groups_dict[int(x)][1] for x in input_list[1:]]
                        combined = add_atomgroups(ag_list).atoms
                        groups_dict[max(groups_dict.keys()) + 1] = ['Group combined from {}'.format(input_list[1:]), combined]
            elif input_list[0] == 'del':
                if len(input_list) > 2:
                    print('Error: Too many arguments.')
                elif len(input_list) < 2:
                    print('Error: Too few arguments.')
                else:
                    try:
                        int(input_list[1])
                    except ValueError:
                        print('Group_ID must be an integer.')
                        continue
                    if int(input_list[1]) not in groups_dict.keys():
                        print('Error: Group_ID must be a valid group ID.')
                    else:
                        groups_dict.pop(int(input_list[1]))

            # help actions
            elif input_list[0] == 'lg':
                if len(input_list) > 1:
                    print('Error: Too many arguments.')
                else:
                    print_db_qr()
                    print_list_groups()
            elif input_list[0] == 'h':
                if len(input_list) > 1:
                    print('Error: Too many arguments.')
                else:
                    print_action_keys()

            # exit actions
            else:
                print('Unknown action key. Please type (h) to see the available action keys.')

        elif input_key != 'e' and len(input_key) == 0:
            continue
        else:
            print('Exiting interactive selection.')
            break
    return target_system
