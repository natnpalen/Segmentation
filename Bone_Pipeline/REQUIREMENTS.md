# Bone Segmentation Pipeline

Automated segmentation of multiple excised-in-air bones from CT scans. Designed for cadaveric hand/wrist specimens (metacarpals, carpals) with embedded lead identification markers. The pipeline takes a CT scan containing multiple bones, finds each bone automatically, segments it into cortical (hard outer shell) and cancellous (spongy interior) regions, and optionally determines how many mechanical test specimens can be cut from each bone. All results are exported as 3D masks and meshes.

## MATLAB Version

- MATLAB R2020b or later (tested on R2024a+)

## Required Toolboxes

| Toolbox | Required? | Purpose |
|---------|-----------|---------|
| **Image Processing Toolbox** | Yes | Core image processing: morphological operations, distance transforms, connected components, fast marching, mesh operations, DICOM/NIfTI I/O |
| **Parallel Computing Toolbox** | No (optional) | Processes multiple bones simultaneously using `parfor`. Without it, bones are processed one at a time — same results, just slower. |

## How to Run

1. Open MATLAB and navigate to `Bone_Pipeline/src/`
2. Edit `run_scan.m` — set `dicomFolder` to your DICOM series folder and `stlFolder` to your mechanical specimen STL folder
3. Run:
   ```matlab
   run_scan
   ```

The pipeline is self-contained — no dependencies on other folders in the repository.

---

## Pipeline Overview

The pipeline runs 6 stages in sequence. A typical scan with 4 bones takes ~2 minutes without packing, or ~30 minutes with packing enabled.

### Stage 1: DICOM Loading (~25-50s)

Reads all DICOM files from the scan folder and assembles them into a 3D volume of Hounsfield Unit (HU) values — a standardized density scale where air is around -1000, water is 0, cancellous bone is 100-400, cortical bone is 400-1500, and lead markers are 4000-7000.

The loader reads each file's header individually and sorts slices by their physical position. This approach is more robust than MATLAB's built-in DICOM functions, which fail on scanner-exported files that use non-standard naming (e.g. hex filenames without `.dcm` extensions). If the folder contains multiple series (e.g. the main CT volume plus a smaller scout/localizer image), it automatically keeps only the series with the most slices.

**Output:** A 3D HU volume with voxel spacing (typically 0.25 x 0.25 x 0.50 mm for these scans) and coordinate transforms.

### Stage 2: Bone Separation (~30-50s)

Finds and isolates each individual bone in the scan. The bones are excised (cut out) and scanned in air, so they appear as bright objects (high HU) against a dark air background (low HU). Each bone has a small lead letter marker attached to it for identification.

**How it works:**

1. **Find the markers** — Lead letters show up at very high HU (>3000). The pipeline finds these, then "grows" the marker mask outward to capture the attached metal flag tabs. This marker mask is used to prevent the bone segmentation from including marker material.

2. **Find seed points** — The volume is thresholded to separate bone-like material from air. Connected regions are found and scored by shape (roundness, elongation) and size. Fragments split by marker exclusion are merged back together if they're within 5mm. One seed point is placed at the deepest interior point of each region.

3. **Grow each bone** — Starting from each seed, the bone region is expanded outward using a Fast Marching Method (FMM) — essentially a "smart flood fill" that follows bone-like densities and avoids markers and air. The growth speed is weighted by how bone-like each voxel is (based on HU) and how far it is from marker artifacts. Multiple growth thresholds are tested and the one producing the best-shaped result (scored by boundary sharpness and interior density) is kept.

4. **Clean up** — The raw bone masks are refined: the outer shell is sealed (small gaps closed), marker material is carved out, low-density surface tissue is scrubbed off, and small disconnected blobs are removed. Non-bone objects (mean HU < 50) are rejected. Each lead marker is associated with its nearest bone.

**Output:** One binary mask per bone, plus marker information and tag associations.

### Stage 3: Cortical / Cancellous Segmentation (~20-40s)

Divides each bone into its two tissue types:
- **Cortical bone** — the dense, hard outer shell (HU typically 400-1500)
- **Cancellous bone** — the spongy, porous interior (HU typically 100-400)

**How it works:**

The bone is divided along its long axis into slabs (4mm wide). Within each slab, the pipeline builds a depth-vs-density profile: starting from the bone surface and moving inward, it measures the average HU at each depth. Cortical bone shows up as a high-density layer near the surface that drops off sharply into lower-density cancellous bone. The boundary is placed at the depth where this density drop is steepest (the maximum negative gradient).

Bones are classified by shape — "elongated" bones like metacarpals get a thicker cortical allowance (up to 2.5mm) while "compact" bones like carpals get a thinner one (up to 1.2mm). The transition depth is smoothed across slabs so the cortical shell varies gradually along the bone's length.

**Output:** Cortical mask, cancellous mask, and metrics (cortical thickness, cortical fraction, bone shape classification).

### Stage 4: Specimen Packing (~20-35 min, optional)

Determines how many mechanical test specimens (Bend, Compression, Punch, Shear) can be physically cut from each bone. This stage is the slowest and can be disabled with `PackSpecimens = false`.

Two packing modes are available:
- **Cortical/cancellous mode** (default) — packs specimens into cortical and cancellous regions separately, so you know how many of each type come from each tissue
- **Whole bone mode** (`PackWholeBone = true`) — ignores the cortical/cancellous boundary and packs into the entire bone volume

**How it works:**

1. **Build templates** — Each specimen STL mesh is loaded, rotated to several orientations aligned with the bone's long axis, and converted to a 3D voxel grid (voxelized) at the scan's resolution.

2. **Find valid positions** — For each template, a 3D convolution slides it across the bone region and measures what fraction of the specimen overlaps with available bone at every position. Positions where at least 95% of the specimen fits inside the bone are considered valid.

3. **Place specimens** — The best position is selected (highest overlap + greatest depth from the bone surface). That space is marked as used, and the search repeats. Priority phase places one of each type first, then a greedy phase fills remaining space.

The pipeline reports which shapes fit and which don't, along with the best overlap percentage achieved for shapes that couldn't be placed. A shape that reports "best overlap 72%" means at best only 72% of the specimen fits inside the bone — the specimen is too large for that bone in every orientation.

**Output:** List of placed specimens with positions, orientations, and tissue classification.

### Stage 5: Visualization (~10-15s)

Generates interactive 3D figures in MATLAB:

1. **Bone Separation** — Each bone shown as a colored 3D surface, labeled with volume, mean HU, and tag ID
2. **Cortical / Cancellous** — Translucent cortical shells over opaque cancellous interiors, labeled with cortical fraction and thickness
3. **Specimen Packing** — Transparent bone outlines with colored specimens placed inside, labeled by type

Figures are also saved as PNG images when output saving is enabled.

### Stage 6: Output Saving (~25-35s)

Writes all results to `bone_pipeline_outputs/<series_name>/<timestamp>/` next to the DICOM folder.

**Per-scan:**
- `pipeline_results.mat` — full MATLAB results struct (for further analysis)
- `pipeline_summary.txt` — human-readable text summary

**Per-bone:**
| File | Description |
|------|-------------|
| `bone_XX_mask.nii.gz` | Binary bone mask — 1 inside bone, 0 outside. NIfTI format, compressed. |
| `bone_XX_cortical.nii.gz` | Cortical (dense shell) region mask |
| `bone_XX_cancellous.nii.gz` | Cancellous (spongy interior) region mask |
| `bone_XX_hu.nii.gz` | HU density values within the bone (non-bone voxels set to -3000) |
| `bone_XX_voxelized.stl` | 3D bone mesh — voxel-accurate surface, minimal smoothing. Useful for measurements. |
| `bone_XX_smooth.stl` | 3D bone mesh — smoothed and decimated for visualization and CAD import. |

The NIfTI files can be opened in 3D Slicer, ITK-SNAP, or similar medical imaging software. The STL files can be opened in SolidWorks, MeshLab, Blender, or any CAD/mesh viewer.

---

## Pipeline Options

Set these as name-value pairs in the `run_bone_pipeline()` call inside `run_scan.m`:

| Option | Default | Description |
|--------|---------|-------------|
| `PackSpecimens` | `true` | Run specimen packing stage. Set to `false` to skip (saves ~30 min). |
| `PackWholeBone` | `false` | Pack into the full bone volume as one region, ignoring cortical/cancellous boundaries. |
| `SaveOutputs` | `true` | Export MAT, NIfTI, and STL files. |
| `ShowViewer` | `true` | Show interactive 3D visualization figures. |
| `PackingOrientations` | `6` | Number of rotations to try per specimen shape. More orientations = better packing but slower. |
| `TagHUMin` | `1200` | HU threshold for metal tag detection. |
| `MinBoneVolMM3` | `500` | Minimum bone volume (mm^3) to keep. Objects smaller than this are discarded. |
| `ClosingRadiusMM` | `3.0` | Morphological closing radius for sealing small gaps in bone masks. |
| `ArtifactSigmaMM` | `3.0` | Controls how far the marker artifact suppression extends from each marker. |
| `TargetIsoMM` | `[]` (off) | Resample to isotropic voxels at this spacing (mm). Leave empty to keep original spacing. |
| `Smoothing` | `false` | Apply edge-preserving smoothing to the HU volume before processing. |
| `OutputDir` | `''` (auto) | Output directory. If empty, auto-creates a timestamped folder next to the DICOM folder. |

---

## File Structure

```
Bone_Pipeline/
  REQUIREMENTS.md          <- this file
  src/
    run_scan.m             <- entry point: set your paths here and run
    run_bone_pipeline.m    <- main pipeline orchestrator (6 stages)
    +dicom/
      series_load.m        <- DICOM CT series loader
    +bone/
      separate_bones.m     <- multi-bone separation (FMM-based)
      cortical_cancellous.m <- cortical/cancellous segmentation
      pack_specimens.m     <- mechanical specimen packing
      visualize_results.m  <- 3D visualization
    +meshing/
      write_stl_binary.m   <- binary STL file writer
    +utils/
      parse_opts.m         <- name-value option parser
```

---

## Input Requirements

### DICOM Folder
- A folder containing CT scan DICOM files (one series, or multiple — the loader picks the dominant one)
- Scanner-exported files without `.dcm` extensions are supported (e.g. hex-named files like `0000004F`)
- Typical scans: micro-CT or clinical CT of excised bone specimens in air

### STL Folder (for specimen packing)
- Must contain one or more of: `Bend.STL`, `Compression.STL`, `Punch.STL`, `Shear.STL`
- These are the mechanical test specimen shapes that will be virtually "cut" from each bone
- Files should be in mm units
- Both `.STL` and `.stl` extensions are accepted
