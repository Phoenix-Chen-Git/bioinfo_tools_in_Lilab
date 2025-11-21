import os
import glob

# Configuration
input_pattern = "*.fa"
output_dir = "split_fastas"

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created directory: {output_dir}")

# Find files
files = glob.glob(input_pattern)
print(f"Processing {len(files)} files...")

count = 0

for file_path in files:
    # Get filename without extension (e.g., rigid_HACE_11)
    base_name = os.path.splitext(os.path.basename(file_path))[0]

    with open(file_path, 'r') as f:
        content = f.read().strip()

    # Split by '>' to separate sequences
    entries = content.split('>')

    for entry in entries:
        if not entry.strip():
            continue

        lines = entry.strip().split('\n')
        header = lines[0]
        sequence = ''.join(lines[1:])  # Join lines in case sequence is wrapped

        # FILTER: Only process if it is a designed sample
        # MPNN original sequences do not have 'sample=' in the header
        if "sample=" in header:

            # Extract sample number (e.g., from "T=0.1, sample=1, score=...")
            parts = header.split(',')
            sample_id = "unknown"
            for part in parts:
                if "sample=" in part:
                    sample_id = part.split('=')[1].strip()
                    break

            # Create new filename: rigid_HACE_11_sample_1.fa
            output_filename = f"{base_name}_sample_{sample_id}.fa"
            output_path = os.path.join(output_dir, output_filename)

            with open(output_path, 'w') as out_f:
                out_f.write(f">{header}\n{sequence}\n")

            count += 1

print(f"Done! Extracted {count} designed sequences to '{output_dir}/'.")