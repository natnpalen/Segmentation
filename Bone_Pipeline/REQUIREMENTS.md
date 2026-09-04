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

There are two entry points: **single scan** (full pipeline) and **batch** (many scans, segmentation only).

### Single scan — `run_scan.m`

1. Open MATLAB and navigate to `Bone_Pipeline/src/`
2. Edit `run_scan.m` — set `dicomFolder` to your DICOM series folder and `stlFolder` to your mechanical specimen STL folder
3. Run:
   ```matlab
   run_scan
   ```

### Batch — `run_batch.m`

1. Open MATLAB and navigate to `Bone_Pipeline/src/`
2. Edit `run_batch.m` — set `rootFolder` to the folder that holds one subfolder per scan
3. Run:
   ```matlab
   run_batch
   ```

See [Batch Mode](#batch-mode) below for what it does and how to configure it.

The pipeline is self-contained — no dependencies on other folders in the repository.

---

## Pipeline Overview

The pipeline runs 6 stages in sequence. A typical scan with 4 bones takes ~2 minutes without packing, or ~30 minutes with packing enabled.

Stages 3, 4 and 5 are optional. Turning off `CorticalCancellous` reduces the run to loading, bone separation and saving — which is what [batch mode](#batch-mode) does, at roughly 1 minute per scan.

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

5. **Judge the result** — Each mask is scored against the material immediately around it (see [Difficult scans](#difficult-scans-the-fallback-pass)). If it fails, a second pass runs under different assumptions and the better of the two is kept.

**Output:** One binary mask per bone, plus marker information, tag associations, and a quality assessment per bone.

### Difficult scans: the fallback pass

The thresholds in stage 2 are tuned for healthy cortical bone. A scan that is low-dose, very osteoporotic, or still carrying soft tissue can defeat them in one of two directions: the HU floors carve away real bone that never reaches them, or the growth swallows tissue that should have been left behind.

Rather than loosen the tuning for everyone — which would cost accuracy on the scans that already work — the pipeline runs the standard pass first, judges the result, and only retries when there is a reason to.

**How a mask is judged.** Absolute HU numbers say little when the whole bone is demineralized, so the mask is measured against its own surroundings:

| Measure | Meaning |
|---------|---------|
| `core` | median HU of the mask interior, more than 1 mm from the surface |
| `rind` | median HU of the non-air material just outside the mask — the tissue we chose not to grade |
| `contrast` | `core - rind`. A mask boundary that isn't backed by a density step isn't a real boundary. |
| `low_frac` | fraction of the mask below `TissueCeilingHU` |
| `fill_ratio` | mask volume as a fraction of the blob it grew from |

The voxel shell touching the mask is skipped when measuring the rind, because partial-volume voxels there are a blend of bone and air and would read as tissue on any specimen.

**Flags** (reported, never silently acted on):

| Flag | Trigger |
|------|---------|
| `low_density` | `core` below `LowDensityHU` — osteoporotic, not necessarily wrong |
| `tissue_suspect` | contrast below `MinContrastHU`, or `low_frac` above `MaxLowFrac` |
| `under_segmented` | `fill_ratio` below `MinFillRatio` — the mask is a fragment of its blob |

**What triggers a retry, and with what:**

| Situation | Second pass | What changes |
|-----------|-------------|--------------|
| No bone found at all | `lowdensity` | HU floors scaled to 0.45, core percentile 94 → 85, growth sweep extended, minimum volume halved |
| `tissue_suspect` | `tissue` | HU floors scaled to 1.25, core percentile → 96, surface tissue scrub 1.7x harder |
| `under_segmented` | `lowdensity` | as above |
| `low_density` **and** `fill_ratio` < 0.75 | `lowdensity` | as above |

`low_density` on its own does **not** trigger a retry — plenty of osteoporotic bone segments cleanly, and the flag alone is not evidence of failure.

**The second pass has to earn its place.** The standard result is kept unless the fallback clearly beats it:

- A `lowdensity` pass is accepted only if it recovered **more than 5% more volume**, its score did not drop, and it did not newly become `tissue_suspect`. A mask that grew by swallowing soft tissue fails all three.
- A `tissue` pass is accepted only if its score improved by more than 0.05 **and** it kept over half the volume. Stripping tissue should trim a mask, not gut it.
- Anything ambiguous keeps the standard result.

**Some scans will still be imperfect.** A severely osteoporotic bone with tissue still attached may have no density step to find, and no threshold setting recovers one. In that case the pipeline keeps the best mask it has and flags it rather than pretending. The flags land in `pipeline_summary.txt`, in the console table, and in `batch_summary.csv` as `quality_flags` and a `review` column — so the scans worth checking by eye can be pulled out directly instead of being trusted silently.

To turn the whole mechanism off, pass `'Fallback', false`. To force one preset for every scan, pass `'FallbackPreset', 'lowdensity'` (or `'tissue'`).

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
| `bone_XX_cortical.nii.gz` | Cortical (dense shell) region mask — only when `CorticalCancellous` is on |
| `bone_XX_cancellous.nii.gz` | Cancellous (spongy interior) region mask — only when `CorticalCancellous` is on |
| `bone_XX_hu.nii.gz` | HU density values within the bone (non-bone voxels set to -3000) |
| `bone_XX_voxelized.stl` | 3D bone mesh — voxel-accurate surface, minimal smoothing. Useful for measurements. |
| `bone_XX_smooth.stl` | 3D bone mesh — smoothed and decimated for visualization and CAD import. |

#### Mesh smoothing

The smooth STL uses **Taubin smoothing**, which alternates a shrinking pass with a slightly larger inflating pass so the surface loses its voxel staircase without pulling in off the bone.

This replaced 15 iterations of plain Laplacian smoothing, which shrank the mesh on every pass. Measured on a 10 mm test sphere carrying a 0.5 mm anatomical ridge at 4 mm wavelength, meshed at CT-like resolution (0.38 mm edges):

| Pass | Ridge detail kept | Surface noise | Volume |
|------|------------------:|--------------:|-------:|
| Unsmoothed input | 100% | 0.119 mm | — |
| Old: Laplacian x15, λ=0.5 | **45%** | 0.147 mm | −1.8% |
| New: Taubin x8, λ=0.5, μ=−0.53 | **99%** | 0.047 mm | ±0.0% |

The old pass was erasing over half of real surface detail at that scale — and because it pulled the surface away from the true shape, it did not even measure as cleaner. The new pass removes 61% of the voxel noise while keeping the detail and the volume.

Bulk volume understates the old behavior: shrinkage scales with curvature, so ridges, the scaphoid waist and other high-curvature features lost far more than the whole-bone figure suggests.

The pipeline prints the smooth STL's enclosed volume against its mask volume after saving, so the effect of any smoothing setting is visible per bone rather than assumed.

Tuning: `MeshSmoothIterations` (default 8; more is smoother, and unlike the old pass it does not shrink), `MeshSmoothLambda`, `MeshPreSmoothSigma` (Gaussian applied to the mask before isosurfacing, 1.0), `MeshDecimate` (fraction of faces kept, 0.5 — was 0.3). Set `MeshSmoothIterations` to 0 for no smoothing at all.

The NIfTI files can be opened in 3D Slicer, ITK-SNAP, or similar medical imaging software. The STL files can be opened in SolidWorks, MeshLab, Blender, or any CAD/mesh viewer.

---

## Pipeline Options

Set these as name-value pairs in the `run_bone_pipeline()` call inside `run_scan.m`:

| Option | Default | Description |
|--------|---------|-------------|
| `CorticalCancellous` | `true` | Run stage 3. Set to `false` for bone masks only — this also disables packing and visualization, which both need the cortical/cancellous split. |
| `MaxBones` | `[]` (all) | Keep only the N largest bones found in the scan. Set to `1` for single-bone scans so stray objects are discarded. |
| `SaveMat` | `true` | Write `pipeline_results.mat`. Set to `false` to skip it — it is by far the largest output file. |
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

### Fallback options

| Option | Default | Description |
|--------|---------|-------------|
| `Fallback` | `true` | Retry difficult scans with a different preset. See [Difficult scans](#difficult-scans-the-fallback-pass). |
| `FallbackPreset` | `''` (auto) | Force a preset instead of choosing one: `'lowdensity'`, `'tissue'`, or `'none'`. |
| `LowDensityHU` | `250` | Interior median below this flags the bone as osteoporotic. |
| `MinContrastHU` | `150` | Minimum density step between bone interior and its surroundings. |
| `MaxLowFrac` | `0.35` | Fraction of the mask allowed below `TissueCeilingHU`. |
| `TissueCeilingHU` | `150` | HU below which a voxel is not clearly bone. |
| `MinFillRatio` | `0.50` | Mask volume / source blob volume, below which the bone is under-segmented. |
| `LowDensityFillRatio` | `0.75` | A low-density bone below this fill ratio is retried. |

### Segmentation knobs

Set these to override the tuning directly. Leave them empty to use the preset values — see `bone.segment_options`.

| Option | Default | Description |
|--------|---------|-------------|
| `DensityScale` | `1.0` | Scales every HU floor in core selection, the growth sweep and boundary refinement. Below 1 keeps low-density bone the defaults would carve away. |
| `CorePrctile` | `94` | Percentile used for the dense-core seed. Lower it when the bone has no dense core. |
| `FMMThreshMax` | `0.42` | Upper end of the growth sweep. Higher grows further before scoring. |
| `TissueScrub` | `1.0` | Scales the surface tissue-removal threshold. Above 1 strips more clinging tissue. |
| `MinBoneHU` | `50` | A candidate whose mean HU is below this is not bone. |

### Mesh options

| Option | Default | Description |
|--------|---------|-------------|
| `MeshSmoothIterations` | `8` | Taubin smoothing passes on the smooth STL. `0` disables smoothing. |
| `MeshSmoothLambda` | `0.5` | Taubin shrink weight. |
| `MeshPreSmoothSigma` | `1.0` | Gaussian sigma applied to the mask before isosurfacing. |
| `MeshDecimate` | `0.5` | Fraction of faces kept in the smooth STL. |

---

## Batch Mode

Batch mode segments every scan under one root folder without any per-scan setup. It is meant for large sets of single-bone scans: it runs **stages 1, 2 and 6 only** — DICOM loading, bone separation, and saving. No cortical/cancellous split, no specimen packing, no figures. Output is NIfTI masks and STL meshes.

### Folder layout

Point `rootFolder` at the folder holding one subfolder per scan. Nested image folders are found automatically, so both of these work:

```
New Bone Scans/                       New Bone Scans/
  156L-1/DICOMOBJ/0000004F ...          156L-1/0000004F ...
  156R-2/DICOMOBJ/...                   156R-2/...
```

A folder counts as a scan when it directly contains at least `MinFiles` DICOM images. Files are recognized by their `DICM` header bytes, so extensionless scanner exports are picked up. Generic wrapper folder names (`DICOMOBJ`, `DICOM`, `IMAGES`, ...) are dropped from the case name, so `156L-1/DICOMOBJ` becomes case `156L-1`.

### Output

Everything lands under `<rootFolder>/bone_pipeline_batch/` (override with `OutputRoot`). By default results are **filed by file type**, pooled across every scan, with HU volumes kept in their own folder apart from the masks. Each file is prefixed with its case name so it stays traceable:

```
bone_pipeline_batch/
  batch_summary.txt        <- per-case status, bone volumes, failures
  batch_summary.csv        <- one row per bone, for Excel/analysis
                              (n_bones_found flags scans where MaxBones
                               discarded extras; pass / quality_flags /
                               review flag scans to check by eye)
  nifti_mask/
    156L-1_bone_01_mask.nii.gz
    156R-2_bone_01_mask.nii.gz
  nifti_hu/
    156L-1_bone_01_hu.nii.gz
    156R-2_bone_01_hu.nii.gz
  stl_smooth/
    156L-1_bone_01_smooth.stl
    156R-2_bone_01_smooth.stl
  stl_voxelized/
    156L-1_bone_01_voxelized.stl
  summaries/
    156L-1_pipeline_summary.txt
```

Both summary files are rewritten after every scan, so a long run can be inspected while it is still going and survives an interrupted session.

### Output layout — the `Organize` option

| Value | Layout |
|-------|--------|
| `'type'` (default) | Pooled by file type across all scans: `<OutputRoot>/nifti_mask/<case>_bone_01_mask.nii.gz`. Best for large batches — every mask, or every STL, in one place. |
| `'case'` | One folder per scan, split by type inside it: `<OutputRoot>/156L-1/nifti_mask/bone_01_mask.nii.gz`. Best when you work one specimen at a time. |
| `'flat'` | One folder per scan with all its files together: `<OutputRoot>/156L-1/bone_01_mask.nii.gz`. Matches the single-scan pipeline's layout. |

Type folders are `nifti_mask`, `nifti_hu`, `stl_smooth`, `stl_voxelized`, `summaries`, plus `nifti_cortical` / `nifti_cancellous` / `mat` / `figures` if you re-enable those stages through `PipelineArgs`. Anything unrecognized goes to `other`.

Scans run into a temporary `_staging/<case>/` folder and are filed once they finish, so a case folder can never collide with a type folder. `_staging` is removed at the end of the run; if it survives, it holds the partial output of scans that failed.

### Behavior on long runs

- **Failures don't stop the batch.** A scan that errors is recorded with its message in the summary and the run continues.
- **Re-running resumes.** Cases that already have a `pipeline_summary.txt` and a bone mask are skipped. Pass `'Overwrite', true` to redo them.
- **Check before committing.** Pass `'DryRun', true` to list the scans that would be processed without running any.

### Batch options

Set these in the `run_batch_pipeline()` call inside `run_batch.m`:

| Option | Default | Description |
|--------|---------|-------------|
| `OutputRoot` | `''` (auto) | Where results are written. Empty = `<rootFolder>/bone_pipeline_batch`. |
| `Organize` | `'type'` | Output layout: `'type'`, `'case'` or `'flat'` — see the table above. |
| `MaxBones` | `1` | Bones to keep per scan. `1` suits single-bone scans; `[]` keeps everything found. |
| `MinFiles` | `5` | Minimum DICOM files for a folder to count as a scan. Lower it for very short series. |
| `Include` | `''` | Regular expression — only run cases whose name matches (e.g. `'^156'`). |
| `Exclude` | `''` | Regular expression — skip cases whose name matches. |
| `Overwrite` | `false` | Re-run cases that already have outputs. Files are replaced in place, not cleared first — if a re-run finds fewer bones than the previous one did, the leftover `bone_02_*` files from the earlier run stay behind. Delete that case's files first if that matters. |
| `DryRun` | `false` | List the discovered scans and stop. |
| `SaveMat` | `false` | Also write the large `pipeline_results.mat` for each case. |
| `PipelineArgs` | `{}` | Extra name-value pairs forwarded to `run_bone_pipeline`, e.g. `{'MinBoneVolMM3', 300}`. |

Example — dry run first, then process only the left-hand specimens with a lower volume floor:

```matlab
run_batch_pipeline(rootFolder, 'DryRun', true);

run_batch_pipeline(rootFolder, ...
    'Include',      'L', ...
    'PipelineArgs', {'MinBoneVolMM3', 300});
```

Batch mode runs scans one at a time; each scan holds a full CT volume in memory. To use more cores, split the root folder and run several MATLAB instances, or use `Include` to partition by name.

---

## File Structure

```
Bone_Pipeline/
  REQUIREMENTS.md          <- this file
  src/
    run_scan.m             <- entry point: single scan, full pipeline
    run_batch.m            <- entry point: batch, segmentation only
    run_bone_pipeline.m    <- main pipeline orchestrator (6 stages)
    run_batch_pipeline.m   <- batch orchestrator (finds and runs every scan)
    +dicom/
      series_load.m        <- DICOM CT series loader
      find_series_dirs.m   <- finds DICOM series folders under a root
    +bone/
      separate_bones.m     <- multi-bone separation (FMM-based)
      segment_options.m    <- density/tissue knobs and fallback presets
      mask_quality.m       <- judges a mask against its surroundings
      quality_defaults.m   <- quality thresholds
      cortical_cancellous.m <- cortical/cancellous segmentation
      pack_specimens.m     <- mechanical specimen packing
      visualize_results.m  <- 3D visualization
    +meshing/
      write_stl_binary.m   <- binary STL file writer
      smooth_mesh_taubin.m <- volume-preserving surface smoothing
    +utils/
      parse_opts.m         <- name-value option parser
      organize_outputs.m   <- files batch outputs into per-type folders
```

---

## Input Requirements

### DICOM Folder
- A folder containing CT scan DICOM files (one series, or multiple — the loader picks the dominant one)
- For batch mode, a root folder holding one such folder per scan (see [Batch Mode](#batch-mode))
- Scanner-exported files without `.dcm` extensions are supported (e.g. hex-named files like `0000004F`)
- Typical scans: micro-CT or clinical CT of excised bone specimens in air

### STL Folder (for specimen packing)
- Must contain one or more of: `Bend.STL`, `Compression.STL`, `Punch.STL`, `Shear.STL`
- These are the mechanical test specimen shapes that will be virtually "cut" from each bone
- Files should be in mm units
- Both `.STL` and `.stl` extensions are accepted
