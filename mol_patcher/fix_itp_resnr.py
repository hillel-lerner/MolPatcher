import sys
import os

def fix_itp_residues(input_file, output_file):
    print(f"Processing {input_file}...")
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        in_atoms = False
        res_counter = 1
        last_seen_raw = None
        modified_count = 0

        for line in infile:
            # Toggle state flag
            if line.strip().startswith("[ atoms ]"):
                in_atoms = True
                outfile.write(line)
                continue
            elif line.strip().startswith("[") and in_atoms:
                in_atoms = False

            # Process data lines inside [ atoms ]
            if in_atoms and line.strip() and not line.strip().startswith(';'):
                parts = line.split()
                raw_res = parts[2]
                
                # Extract the base integer to detect chain resets (e.g., '1A' -> 1)
                base_int = int("".join(filter(str.isdigit, raw_res)))

                # Mimic the PDB's state machine logic:
                if last_seen_raw is None:
                    # Initialize the very first atom
                    last_seen_raw = raw_res
                    res_counter = 1
                elif base_int == 1 and raw_res != last_seen_raw:
                    # Chain reset detected, Restart the counter to 1
                    res_counter = 1
                    last_seen_raw = raw_res
                    print(f"Chain reset detected. Resetting counter to 1.")
                elif raw_res != last_seen_raw:
                    # Residue changed (e.g., '52' -> '52A' or '52' -> '53'). Increment counter.
                    res_counter += 1
                    last_seen_raw = raw_res

                # Overwrite the column with our sequential counter
                parts[2] = str(res_counter)
                
                # Format and save
                formatted_line = f"{parts[0]:>6} {parts[1]:>10} {parts[2]:>6} {parts[3]:>6} {parts[4]:>6} {parts[5]:>6} {parts[6]:>12} {parts[7]:>12}"
                if len(parts) > 8:
                    formatted_line += "   " + " ".join(parts[8:])
                
                outfile.write(formatted_line + "\n")
                modified_count += 1
            else:
                # Copy all other lines directly
                outfile.write(line)

    print(f"Done. Processed {modified_count} atoms. Saved clean topology to {output_file}")

if __name__ == "__main__":
    input_itp = "topology_files/PROBC.itp"
    output_itp = "topology_files/PROBC_fixed.itp"
    
    if not os.path.exists(input_itp):
        print(f"Error: {input_itp} not found.")
        sys.exit(1)
        
    fix_itp_residues(input_itp, output_itp)