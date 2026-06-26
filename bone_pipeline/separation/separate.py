"""
Bone separation for excised specimens scanned in air.

For excised bones in air, the non-air non-tag region IS the bone.
The scaphoid pipeline uses FMM-based segmentation tuned for a single
bone per scan — here we need to separate multiple bones, and the
envelope approach (specimen minus tags = bone) is simpler and avoids
the reconstruction/FMM machinery that can under-detect less-dense
bone ends when thresholds are tuned for a single-bone scenario.

Pipeline:
  1. Specimen isolation: non-air material, physical-space closing to
     bridge internal porosity, 3D hole-fill
  2. Tag detection: lead markers (HU > 1200) with exclusion zone
  3. Bone detection: specimen minus tags, split into components
  4. Per-bone interior fill (per-slice 2D + 3D) to capture marrow cavities
  5. Validation: each component must have some dense bone tissue
  6. Tag-to-bone association by proximity
"""

import numpy as np
from pathlib import Path
from scipy import ndimage
from skimage.measure import regionprops
from skimage.morphology import ball

from .dicom_io import load_dicom_series


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def separate_bones(dicom_folder, tag_hu_min=1200, min_bone_volume_mm3=200.0,
                   closing_radius_mm=3.0):
    """Separate individual bones from a multi-bone DICOM CT scan.

    For excised bones scanned in air, uses envelope detection:
    the non-air non-tag region is the bone.

    Parameters
    ----------
    dicom_folder : str or Path
        Folder containing DICOM files for one scan.
    tag_hu_min : float
        HU threshold for lead tag detection (default 1200).
    min_bone_volume_mm3 : float
        Minimum volume for a component to count as a bone.
    closing_radius_mm : float
        Morphological closing radius (mm) for bridging internal porosity.

    Returns
    -------
    dict  with keys 'volume', 'spacing', 'bones' (list of bone dicts).
    """
    dicom_folder = Path(dicom_folder)
    volume, spacing = load_dicom_series(dicom_folder)

    if volume.ndim != 3:
        raise ValueError(f"Expected 3-D volume, got shape {volume.shape}.")

    vol = volume.astype(np.float32)
    voxel_vol = float(np.prod(spacing))

    # === Stage 1: Specimen isolation ===
    print("  Stage 1: Specimen isolation...")
    specimen = _isolate_specimen(vol, spacing, closing_radius_mm)
    spec_vol = np.sum(specimen) * voxel_vol
    print(f"    Specimen: {spec_vol:.0f} mm³ "
          f"({np.sum(specimen)} voxels)")

    # === Stage 2: Tag detection ===
    print("  Stage 2: Tag detection...")
    lead_mask = vol > tag_hu_min
    tag_components = _find_tags(lead_mask, spacing, voxel_vol)
    print(f"    Found {len(tag_components)} metal tags")

    # Build a mask that excludes tags and their immediate proximity.
    # Tags are near bones but not touching — a small exclusion zone
    # removes metal artifacts without cutting into bone.
    if np.any(lead_mask):
        tag_dist = ndimage.distance_transform_edt(~lead_mask,
                                                   sampling=spacing)
        tag_exclusion = tag_dist < 1.0  # 1 mm exclusion around tags
    else:
        tag_exclusion = np.zeros_like(vol, dtype=bool)

    # === Stage 3: Bone detection ===
    print("  Stage 3: Bone envelope detection...")

    # For excised bones in air, the bone is everything in the specimen
    # that isn't a tag or tag artifact.
    bone_region = specimen & ~tag_exclusion

    # Split into connected components (each should be one bone)
    labeled, n_comp = ndimage.label(bone_region)
    raw_vol = np.sum(bone_region) * voxel_vol
    print(f"    Raw bone region: {raw_vol:.0f} mm³, "
          f"{n_comp} components")

    # === Stage 4: Per-bone fill and validation ===
    print("  Stage 4: Per-bone fill & validation...")

    bones = []
    small_count = 0

    for i in range(1, n_comp + 1):
        comp = labeled == i
        comp_vol = np.sum(comp) * voxel_vol

        if comp_vol < min_bone_volume_mm3:
            small_count += 1
            continue

        # Validate: must have some dense bone tissue (not just noise)
        n_dense = int(np.sum(vol[comp] > 200))
        dense_frac = n_dense / max(1, np.sum(comp))
        if dense_frac < 0.02:
            print(f"    Component {i}: {comp_vol:.0f} mm³ — "
                  f"skipped (only {dense_frac:.1%} dense)")
            continue

        mean_hu = float(np.mean(vol[comp]))
        if mean_hu > tag_hu_min:
            print(f"    Component {i}: {comp_vol:.0f} mm³ — "
                  f"skipped (mean HU {mean_hu:.0f}, likely tag)")
            continue

        # Per-slice 2D fill to capture enclosed marrow cavities.
        # In each axial cross-section, the cortical shell forms a ring
        # that encloses the marrow cavity.
        filled = comp.copy()
        for z in range(filled.shape[0]):
            if np.any(filled[z]):
                filled[z] = ndimage.binary_fill_holes(filled[z])
        filled = ndimage.binary_fill_holes(filled)

        # Stay within the specimen bounds
        filled = filled & specimen

        filled_vol = np.sum(filled) * voxel_vol
        filled_hu = float(np.mean(vol[filled]))

        # Compute centroid in mm
        coords = np.argwhere(filled)
        centroid_vox = coords.mean(axis=0)
        centroid_mm = centroid_vox * np.array(spacing)

        # Bounding box
        z_min, y_min, x_min = coords.min(axis=0)
        z_max, y_max, x_max = coords.max(axis=0)
        bbox = (int(z_min), int(y_min), int(x_min),
                int(z_max), int(y_max), int(x_max))

        bones.append({
            "label": i,
            "centroid": centroid_mm,
            "mask": filled,
            "bbox": bbox,
            "volume_mm3": filled_vol,
            "mean_hu": filled_hu,
            "dense_fraction": dense_frac,
        })

        print(f"    Bone: {comp_vol:.0f} → {filled_vol:.0f} mm³ "
              f"(fill +{filled_vol - comp_vol:.0f}), "
              f"mean HU {filled_hu:.0f}, "
              f"dense {dense_frac:.0%}")

    if small_count > 0:
        print(f"    Filtered {small_count} small components "
              f"(< {min_bone_volume_mm3} mm³)")

    # === Stage 5: Tag association ===
    _associate_tags(bones, tag_components)
    bones.sort(key=lambda b: b["volume_mm3"], reverse=True)

    print(f"\nFound {len(bones)} bones and {len(tag_components)} tags in scan")
    for i, b in enumerate(bones):
        tag = f"tag {b['tag_id']}" if b.get("tag_id") else "no tag"
        print(f"  Bone {i + 1}: {b['volume_mm3']:.1f} mm³, "
              f"mean HU {b['mean_hu']:.0f}, {tag}")

    return {"volume": volume, "spacing": spacing, "bones": bones}


# ---------------------------------------------------------------------------
# Helpers — morphological operations
# ---------------------------------------------------------------------------

def _physical_close(mask, radius_mm, spacing):
    """Morphological closing with correct physical radius.

    Uses distance transforms instead of structuring elements so the
    closing sphere is truly round in physical space regardless of
    voxel anisotropy.
    """
    dt_bg = ndimage.distance_transform_edt(~mask, sampling=spacing)
    dilated = dt_bg <= radius_mm
    dt_fg = ndimage.distance_transform_edt(dilated, sampling=spacing)
    return dt_fg >= radius_mm


# ---------------------------------------------------------------------------
# Helpers — specimen & markers
# ---------------------------------------------------------------------------

def _isolate_specimen(vol, spacing, closing_radius_mm=3.0):
    """All non-air material with physical-space closing.

    Generous closing bridges internal trabecular porosity so each bone
    forms a connected region.  3D fill captures fully enclosed cavities.
    Per-slice fill is NOT done here to avoid merging adjacent bones in
    slices where they're close — that's done per-bone after splitting.
    All components above 50 mm³ are kept.
    """
    non_air = vol > -500
    non_air = _physical_close(non_air, closing_radius_mm, spacing)
    non_air = ndimage.binary_fill_holes(non_air)
    voxel_vol = float(np.prod(spacing))
    min_vox = max(100, int(50.0 / voxel_vol))
    non_air = _remove_tiny(non_air, min_vox)
    return non_air


def _find_tags(lead_mask, spacing, voxel_vol):
    if not np.any(lead_mask):
        return []
    labeled, n = ndimage.label(lead_mask)
    tags = []
    for prop in regionprops(labeled):
        tags.append({
            "label": prop.label,
            "centroid": np.array(prop.centroid) * np.array(spacing),
            "volume_mm3": prop.area * voxel_vol,
        })
    return tags


def _associate_tags(bones, tags):
    for b in bones:
        b["tag_id"] = None
        b["tag_dist"] = None
    if not tags or not bones:
        return
    centroids = np.array([b["centroid"] for b in bones])
    for tag in tags:
        dists = np.linalg.norm(centroids - tag["centroid"], axis=1)
        idx = int(np.argmin(dists))
        d = float(dists[idx])
        if bones[idx]["tag_dist"] is None or d < bones[idx]["tag_dist"]:
            bones[idx]["tag_id"] = tag["label"]
            bones[idx]["tag_dist"] = d


# ---------------------------------------------------------------------------
# Helpers — utilities
# ---------------------------------------------------------------------------

def _remove_tiny(mask, min_voxels):
    """Remove connected components smaller than *min_voxels*."""
    labeled, n = ndimage.label(mask)
    if n == 0:
        return mask
    sizes = ndimage.sum(mask, labeled, range(1, n + 1))
    keep = np.zeros(n + 1, dtype=bool)
    for i in range(n):
        if sizes[i] >= min_voxels:
            keep[i + 1] = True
    return keep[labeled]
