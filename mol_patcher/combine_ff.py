import os

def build_master_forcefield_complete(base_path, pfp_path, pi_path, output_path):
    mods = {}
    
    def add_mod(sec, line, tag):
        # Prevent proper (9) and improper (2) dihedrals from mixing
        if sec == "[ dihedraltypes ]":
            parts = line.split(';')[0].split()
            if len(parts) >= 5:
                func = parts[4]
                if func in ['1', '9']:
                    sec = "[ dihedraltypes_proper ]"
                elif func in ['2', '4']:
                    sec = "[ dihedraltypes_improper ]"
        
        if sec not in mods:
            mods[sec] = []
        
        modified_line = line.split(';')[0].rstrip() + f" ; {tag}\n"
        mods[sec].append(modified_line)

    # 1. Load PFP base parameters
    current_sec = None
    with open(pfp_path, "r") as f_pfp:
        for line in f_pfp:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                current_sec = stripped
            elif stripped and not stripped.startswith(";"):
                if current_sec and current_sec != "[ defaults ]":
                    add_mod(current_sec, line, "added from pfp base")

    # 2. Load PI's modifications
    current_sec = None
    with open(pi_path, "r") as f_pi:
        for line in f_pi:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                current_sec = stripped
            elif "MM" in line and current_sec:
                add_mod(current_sec, line, "added by HL (link)")

    # 3. Inject explicit atom types
    extra_atomtypes = [
        " CLA       17   35.450000   -1.00    A   0.404468018036      0.6276       ; explicit addition\n",
        " HT         1   1.0080      0.417    A   4.00013524445e-02   1.924640e-01 ; explicit addition\n",
        " NG311      7   14.0070    -0.092    A   3.56359487256e-01   1.882800e-01 ; explicit addition\n",
        " OT         8   15.9994    -0.834    A   3.15057422683e-01   6.363864e-01 ; explicit addition\n",
        " POT       19   39.102000   1.00     A   0.314264522824      0.364008     ; explicit addition\n",
        " SOD       11   22.989770   1.000    A   0.251367073323      0.19623      ; explicit addition\n"
    ]
    if "[ atomtypes ]" not in mods:
        mods["[ atomtypes ]"] = []
    mods["[ atomtypes ]"].extend(extra_atomtypes)

    # 4. Merge, deduplicate, sort, and write
    with open(base_path, "r") as f_base, open(output_path, "w") as f_out:
        current_sec = None
        active_mod_key = None
        section_lines = []
        dihedral_count = 0
        
        def write_sorted_section():
            if current_sec:
                combined_lines = section_lines + mods.get(active_mod_key, [])
                
                # CRITICAL FIX: No deduplication and no sorting for matrices/defaults!
                if current_sec in ["[ defaults ]", "[ cmaptypes ]"]:
                    unique_lines = combined_lines
                else:
                    seen_standard = set()
                    unique_lines = []
                    for l in combined_lines:
                        data_only = l.split(';')[0].strip()
                        standardized = " ".join(data_only.split())
                        if not standardized or standardized not in seen_standard:
                            if standardized:
                                seen_standard.add(standardized)
                            unique_lines.append(l)
                    unique_lines.sort()

                for l in unique_lines:
                    f_out.write(l)

        for line in f_base:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                write_sorted_section()
                current_sec = stripped
                section_lines = []
                
                # Track which dihedral section we are currently in
                if current_sec == "[ dihedraltypes ]":
                    dihedral_count += 1
                    if dihedral_count == 1:
                        active_mod_key = "[ dihedraltypes_proper ]"
                    else:
                        active_mod_key = "[ dihedraltypes_improper ]"
                else:
                    active_mod_key = current_sec
                    
                f_out.write(line)
            elif stripped.startswith(";") or not stripped:
                f_out.write(line)
            else:
                section_lines.append(line)
                
        write_sorted_section()

# Execute
build_master_forcefield_complete(
    base_path="forcefields/6oge_all_forcefield.itp", 
    pfp_path="forcefields/ff_pfp_donor.itp",
    pi_path="forcefields/ff_PRODE_E136.itp", 
    output_path="forcefields/forcefield_master.itp"
)