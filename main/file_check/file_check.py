from pathlib import Path
from expected_paths_static import expected_paths

# ---- CONFIGURATION ----
NEW_SNAPSHOT_FILENAME = "1_make_templates.txt"  # Change this name before each run
# -----------------------

ROOT = Path(__file__).resolve().parent.parent
GEN_PATH = Path(__file__).parent / "generated_files"
GEN_PATH.mkdir(exist_ok=True)

# 1. Load previously recorded file paths
previous_files = set()
for txt_file in GEN_PATH.glob("*.txt"):
    with txt_file.open("r") as f:
        for line in f:
            if line.strip():
                previous_files.add(Path(line.strip()).resolve())

# 2. Prepare current expected paths, removing already handled ones
expected_paths_set = set((ROOT / p).resolve() for p in expected_paths)
expected_paths_set -= previous_files

# 3. Search for all .pdf and .png files under ROOT
all_found_files = set()
for ext in ('*.pdf', '*.png'):
    all_found_files.update(p.resolve() for p in ROOT.rglob(ext))

# 4. Compute status categories
found = sorted(p for p in expected_paths_set if p in all_found_files)
missing = sorted(p for p in expected_paths_set if p not in all_found_files)
extra = sorted(p for p in all_found_files if p not in expected_paths_set and p not in previous_files)

# 5. Save newly matched files to NEW_SNAPSHOT_FILENAME
outfile = GEN_PATH / NEW_SNAPSHOT_FILENAME
with outfile.open("w") as f:
    for p in found:
        f.write(str(p) + "\n")

# 6. Print summaries
print(f"\nSaved {len(found)} new matched file(s) to: {outfile}")

if found:
    print("\nNewly matched files:")
    for fpath in found:
        print(f"  + {fpath}")

print(f"\nExtra: {len(extra)} unexpected file(s)")
for f in extra:
    print(f"  - {f}")