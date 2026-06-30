# Bone Segmentation Pipeline — Requirements

## MATLAB Version
- MATLAB R2020b or later (tested on R2024a+)

## Required Toolboxes

| Toolbox | Used For |
|---------|----------|
| **Image Processing Toolbox** | `bwdist`, `bwconncomp`, `regionprops3`, `imclose`, `imopen`, `imerode`, `imdilate`, `imfill`, `imreconstruct`, `bwareaopen`, `bwperim`, `strel`, `imgradient3`, `imsegfmm`, `imresize3`, `isosurface`, `smooth3`, `mat2gray`, `gradientweight`, `imdiffusefilt`, `stlread`, `reducepatch`, `niftiwrite`, `niftiinfo`, `dicominfo`, `dicomread` |
| **Parallel Computing Toolbox** *(optional)* | `parfor` — used to process bones in parallel. Pipeline works without it (falls back to serial `for` loops). |

## How to Run

1. Open MATLAB and navigate to `Bone_Pipeline/src/`
2. Edit `run_scan.m` — set `dicomFolder` and `stlFolder` to your paths
3. Run:
   ```matlab
   run_scan
   ```

## Pipeline Options

Set these as name-value pairs in the `run_bone_pipeline()` call in `run_scan.m`:

| Option | Default | Description |
|--------|---------|-------------|
| `PackSpecimens` | `true` | Run specimen packing (slow) |
| `SaveOutputs` | `true` | Export MAT, NIfTI, and STL files |
| `ShowViewer` | `true` | Show 3D visualization figures |
| `PackingOrientations` | `6` | Number of orientations per specimen shape |
| `TagHUMin` | `1200` | HU threshold for metal tag detection |
| `MinBoneVolMM3` | `500` | Minimum bone component volume (mm^3) |

## Outputs

When `SaveOutputs` is enabled, the pipeline writes to `bone_pipeline_outputs/` next to the DICOM folder:

- `pipeline_results.mat` — full results struct
- `pipeline_summary.txt` — text summary
- `bone_XX_mask.nii.gz` — binary bone mask (NIfTI)
- `bone_XX_cortical.nii.gz` — cortical region mask
- `bone_XX_cancellous.nii.gz` — cancellous region mask
- `bone_XX_hu.nii.gz` — HU values within bone
- `bone_XX_voxelized.stl` — bone mesh (voxel-accurate)
- `bone_XX_smooth.stl` — bone mesh (smoothed, decimated)
