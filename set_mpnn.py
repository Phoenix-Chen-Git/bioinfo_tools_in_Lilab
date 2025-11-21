import os
import glob
import shutil

# ==========================================
# CONFIGURATION: SLURM HEADER TEMPLATE
# ==========================================
SLURM_HEADER = """#!/bin/bash
#SBATCH -J MPNN
#SBATCH -p gpu_l40
#SBATCH -N 1
#SBATCH -o RFdiff_%j.out
#SBATCH -e RFdiff_%j.err
#SBATCH --no-requeue
#SBATCH -A yulongli_g1
#SBATCH --qos=yulonglil40
#SBATCH --gres=gpu:1
#SBATCH --overcommit
#SBATCH --mincpus=4

source ~/.bashrc
source activate mlfold
"""

# ==========================================
# CONFIGURATION: BASH BLOCK TEMPLATE
# ==========================================
# Note: We use {{variables}} for Python formatting.
# The $variables are literal bash variables.
BLOCK_TEMPLATE = """
# --------------------------------------------------------------
# Group: Length {g_len}, Start Residue {g_start}
# --------------------------------------------------------------
folder_with_pdbs="{group_folder_path}"

output_dir="{output_base_dir}/output_len_{g_len}_start_{g_start}"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="A"
# The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
design_only_positions="{positions_str}" # design only these residues; use flag --specify_non_fixed

python /home/yulongli_pkuhpc/lustre3/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python /home/yulongli_pkuhpc/lustre3/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python /home/yulongli_pkuhpc/lustre3/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python /lustre3/yulongli_pkuhpc/ProteinMPNN/protein_mpnn_run.py \\
        --jsonl_path $path_for_parsed_chains \\
        --chain_id_jsonl $path_for_assigned_chains \\
        --fixed_positions_jsonl $path_for_fixed_positions \\
        --out_folder $output_dir \\
        --num_seq_per_target 6 \\
        --sampling_temp "0.1" \\
        --seed 37 \\
        --batch_size 1
"""


def get_longest_gly_pattern(pdb_path):
    """
    Parses PDB to find the longest continuous GLY pattern.
    Returns: (length, list_of_residue_ids)
    """
    residues = []
    with open(pdb_path, 'r') as f:
        seen = set()
        for line in f:
            if line.startswith("ENDMDL"): break
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                chain = line[21]
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                i_code = line[26].strip()
                full_id = res_seq + i_code

                key = (chain, full_id)
                if key in seen: continue
                seen.add(key)

                residues.append({'name': res_name, 'id': full_id, 'chain': chain})

    if not residues: return 0, []

    max_len = 0
    best_ids = []
    curr_len = 0
    curr_ids = []
    last_chain = None

    for res in residues:
        if res['chain'] != last_chain:
            curr_len = 0
            curr_ids = []
            last_chain = res['chain']

        if res['name'] == 'GLY':
            curr_len += 1
            curr_ids.append(res['id'])
            if curr_len > max_len:
                max_len = curr_len
                best_ids = list(curr_ids)
        else:
            curr_len = 0
            curr_ids = []

    return max_len, best_ids


def main():
    # 1. Setup Directories
    base_input_dir = os.path.abspath(os.getcwd())  # Current dir assumed to have PDBs
    pdb_files = glob.glob("*.pdb")

    if not pdb_files:
        print("No .pdb files found here.")
        return

    # 2. Parse and Group Files
    # We group by (Length, Tuple(Positions)) to ensure MPNN runs correctly.
    # If we only grouped by length, different sequence positions would conflict in one batch.
    groups = {}

    print(f"Processing {len(pdb_files)} files...")

    for pdb in pdb_files:
        length, positions = get_longest_gly_pattern(pdb)

        if length == 0:
            print(f"Skipping {pdb} (No Glycines found)")
            continue

        # Create a unique key for grouping
        # positions is a list like ['93', '94', '95']
        pos_tuple = tuple(positions)
        group_key = (length, pos_tuple)

        if group_key not in groups:
            groups[group_key] = []
        groups[group_key].append(pdb)

    # 3. Move Files and Generate Bash Script
    bash_script_content = SLURM_HEADER

    # Create a base output directory for the results referenced in the bash script
    # The user specified "../outputs/..." in the bash text, we will use that structure.
    bash_output_base = "/home/yulongli_pkuhpc/lustre3/HACE/rigid_linker/design1/MPNNresults"

    for (length, pos_tuple), files in groups.items():
        start_pos = pos_tuple[0]

        # Create specific folder name: G_len_X_start_Y
        # This ensures that all PDBs in this folder have the exact same design positions
        dir_name = f"G_len_{length}_start_{start_pos}"
        full_dir_path = os.path.join(base_input_dir, dir_name)

        if not os.path.exists(full_dir_path):
            os.makedirs(full_dir_path)

        # Generate the positions string (e.g., "93 94 95")
        pos_str = " ".join(pos_tuple)

        # Move files and create txt report
        txt_path = os.path.join(full_dir_path, "positions.txt")
        with open(txt_path, "w") as txt:
            for f in files:
                shutil.move(f, os.path.join(full_dir_path, f))
                txt.write(f"{f}: {pos_str}\n")

        print(f"Created group {dir_name} with {len(files)} files.")

        # Append to Bash Script
        # Note: The python script in the bash block assumes it is running from a specific location.
        # We set 'folder_with_pdbs' to the absolute path we just created to be safe.
        bash_block = BLOCK_TEMPLATE.format(
            g_len=length,
            g_start=start_pos,
            group_folder_path=full_dir_path,
            output_base_dir=bash_output_base,
            positions_str=pos_str
        )
        bash_script_content += bash_block

    # 4. Write the final run_mpnn.sh
    with open("run_mpnn.sh", "w") as f:
        f.write(bash_script_content)

    print("\nDone!")
    print("1. PDBs have been grouped into subdirectories.")
    print("2. 'run_mpnn.sh' has been generated.")
    print("3. Check 'positions.txt' in each subfolder.")


if __name__ == "__main__":
    main()