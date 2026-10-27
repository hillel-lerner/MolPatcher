def validate_names(pdb_list, topol_list, config_names, label):
    pdb_list = set(pdb_list)
    topol_list = set(topol_list)
    config_names = set(config_names)

    outputs = []
    for n in sorted(config_names):
        if n in pdb_list and n not in topol_list:
            message = f"[{label}] {n} missing from topology"
            outputs.append(message)
        if n in topol_list and n not in pdb_list:
            message = f"[{label}] {n} missing from pdb"
            outputs.append(message)
        if n not in topol_list and n not in pdb_list:
            message = f"[{label}] {n} missing from topology and pdb"
            outputs.append(message)

    return outputs
