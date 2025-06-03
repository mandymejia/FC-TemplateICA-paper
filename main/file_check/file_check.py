from pathlib import Path
from datetime import datetime

# Constants
ROOT = Path(__file__).resolve().parents[1]  # points to /main
GEN_PATH = Path(__file__).parent / "generated_files"
GEN_PATH.mkdir(exist_ok=True)

# --- Load all previously recorded file paths ---
previous_files = set()
for txt_file in GEN_PATH.glob("*.txt"):
    with txt_file.open("r") as f:
        for line in f:
            if line.strip():
                previous_files.add(Path(line.strip()).resolve())

# --- Take a fresh snapshot ---
current_files = set()
for ext in ('*.pdf', '*.png'):
    current_files.update(p.resolve() for p in ROOT.rglob(ext))

# --- Find new files not previously recorded ---
new_files = sorted(current_files - previous_files)

# --- Save snapshot with timestamp ---

outfile = GEN_PATH / f"1_make_templates.txt"
with outfile.open("w") as f:
    for p in sorted(new_files):
        f.write(str(p) + "\n")

# --- Print summary ---
print(f"Found {len(current_files)} total .pdf/.png files.")
print(f"Detected {len(new_files)} new file(s) not in any prior snapshot.")

if new_files:
    print("\nNew files:")
    for f in new_files:
        print(f"  + {f}")
    print(f"\nSaved snapshot of new files to: {outfile}")
else:
    print("No new files detected.")