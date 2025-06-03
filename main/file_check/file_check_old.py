from pathlib import Path
from expected_paths_static import expected_paths

# Define the root (parent directory of the script)
ROOT = Path(__file__).resolve().parent.parent

# Normalize expected paths relative to ROOT
expected_paths_set = set((ROOT / p).resolve() for p in expected_paths)

# Search for all .pdf and .png files under the parent directory
all_found_files = set()
for ext in ('*.pdf', '*.png'):
    for file in ROOT.rglob(ext):
        all_found_files.add(file.resolve())

# Found and missing from expected list
found = [str(p) for p in expected_paths_set if p in all_found_files]
missing = [str(p) for p in expected_paths_set if p not in all_found_files]
extra = [str(p) for p in all_found_files if p not in expected_paths_set]

# Summary
print(f"\n Found: {len(found)} files")
for f in found:
    print(f"  ✔ {f}")

# print(f"\n Missing: {len(missing)} files")
# for f in missing:
#     print(f"  ✘ {f}")

print(f"\n🧹 Extra: {len(extra)} unexpected files")
for f in extra:
    print(f"  ⚠ {f}")