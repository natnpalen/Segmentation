"""
Bone separation module — adapted from the scaphoid pipeline's approach.

Pipeline (mirrors the proven scaphoid segmentation strategy):
  1. Air detection: HU < -500 → identify non-air specimen region
  2. Metal tag detection: HU > 1200 → mask tags with Gaussian falloff
  3. Specimen isolation: non-air, morphologically closed, largest component
  4. Bone seed finding: dense interior points within specimen
  5. Morphological reconstruction: grow from seeds through allow region
     where voxels are included if HU > low_threshold OR gradient is high
  6. Interior protection: deep interior (>1mm) never removed
  7. Connected components to separate individual bones
  8. Tag-to-bone association by proximity
"""

import numpy as np
from pathlib import Path
from scipy import ndimage
from skimage.filters import threshold_otsu
from skimage.measure import regionprops, label as sk_label
from skimage.morphology import ball, binary_closing, binary_opening
from skimage.segmentation import watershed

from .dicom_io import load_dicom_series


def separate_bones(dicom_folder, tag_hu_min=1200, min_bone_volume_mm3=200.0,
                   closing_radius_mm=2.0, bone_hu_floor=0):
    """Separate individual bones from a multi-bone DICOM CT scan.

    Uses a multi-stage approach adapted from the scaphoid pipeline:
    air detection → tag masking → specimen isolation → seed-based growth.

    Parameters
    ----------
    dicom_folder : str or Path
        Path to folder containing DICOM files for one scan.
    tag_hu_min : float
        HU threshold for lead tag detection (default 1200).
    min_bone_volume_mm3 : float
        Minimum volume for a component to be considered a bone.
    closing_radius_mm : float
        Morphological closing radius in mm for bridging trabecular gaps.
    bone_hu_floor : float
        Minimum HU for bone tissue. Voxels above this within the specimen
        (and not tags) are candidates for bone. Default 0 (above water).

    Returns
    -------
    dict with keys:
        'volume', 'spacing', 'bones' (list of bone dicts)
    """
    dicom_folder = Path(dicom_folder)
    volume, spacing = load_dicom_series(dicom_folder)

    if volume.ndim != 3:
        raise ValueError(
            f"Expected 3D volume, got shape {volume.shape}.")

    voxel_vol_mm3 = float(np.prod(spacing))
    mean_spacing = float(np.mean(spacing))

    # --- Stage 1: Air detection & specimen isolation ---
    print("  Stage 1: Air detection & specimen isolation...")
    air_mask = volume < -500
    non_air = ~air_mask

    close_r = max(1, int(round(1.0 / mean_spacing)))
    non_air = ndimage.binary_closing(non_air, structure=ball(close_r))

    labeled_spec, n_spec = ndimage.label(non_air)
    if n_spec > 0:
        sizes = ndimage.sum(non_air, labeled_spec, range(1, n_spec + 1))
        largest = int(np.argmax(sizes)) + 1
        specimen_mask = labeled_spec == largest
    else:
        specimen_mask = non_air

    specimen_vol = np.sum(specimen_mask) * voxel_vol_mm3
    print(f"    Specimen volume: {specimen_vol:.0f} mm³")

    # --- Stage 2: Metal tag detection ---
    print("  Stage 2: Metal tag detection...")
    tag_mask = specimen_mask & (volume > tag_hu_min)

    # Dilate tags slightly to capture immediate surroundings
    tag_dilated = ndimage.binary_dilation(tag_mask, structure=ball(2))

    # Build artifact weight map (Gaussian falloff from tags)
    if np.any(tag_mask):
        tag_dist = ndimage.distance_transform_edt(~tag_mask, sampling=spacing)
        artifact_sigma_mm = 3.0
        artifact_weight = np.exp(-(tag_dist / artifact_sigma_mm) ** 2)
    else:
        artifact_weight = np.zeros_like(volume, dtype=float)

    n_tags_found = 0
    tag_components = []
    if np.any(tag_mask):
        tag_labeled, n_tags_found = ndimage.label(tag_mask)
        tag_props = regionprops(tag_labeled)
        for prop in tag_props:
            tag_components.append({
                'label': prop.label,
                'centroid': np.array(prop.centroid) * np.array(spacing),
                'mask': tag_labeled == prop.label,
                'volume_mm3': prop.area * voxel_vol_mm3,
            })

    print(f"    Found {n_tags_found} metal tags")

    # --- Stage 3: Bone detection with gradient-assisted growth ---
    print("  Stage 3: Bone detection (seed-based growth)...")

    # Compute gradient magnitude for edge detection
    grad = _compute_gradient_magnitude(volume, spacing)

    # Allow region: within specimen, above bone floor OR high gradient
    # This is the key insight from the scaphoid pipeline — OR logic
    # includes trabecular bone that has good boundary gradients
    specimen_no_tags = specimen_mask & (~tag_dilated)

    grad_thresh = np.percentile(grad[specimen_no_tags], 85)
    hu_allow = specimen_no_tags & (
        (volume > bone_hu_floor) | (grad > grad_thresh)
    )

    allow_vol = np.sum(hu_allow) * voxel_vol_mm3
    print(f"    Allow region (HU>{bone_hu_floor} OR grad>{grad_thresh:.1f}): "
          f"{allow_vol:.0f} mm³")

    # Find dense core seeds (high HU within specimen, away from tags)
    core_hu_thresh = np.percentile(volume[specimen_no_tags], 90)
    core_hu_thresh = max(200, min(700, core_hu_thresh))
    core_seeds = specimen_no_tags & (volume > core_hu_thresh)

    # Remove small core fragments
    core_seeds = _remove_small_components(core_seeds, 10)

    print(f"    Core seed threshold: {core_hu_thresh:.0f} HU, "
          f"seed volume: {np.sum(core_seeds) * voxel_vol_mm3:.0f} mm³")

    # Morphological reconstruction: grow from core seeds through allow region
    # This is equivalent to the scaphoid pipeline's imreconstruct
    bone_mask = _morphological_reconstruct(core_seeds, hu_allow)

    # Morphological closing to bridge remaining trabecular gaps
    close_r_bone = max(1, int(round(closing_radius_mm / mean_spacing)))
    bone_mask = ndimage.binary_closing(bone_mask, structure=ball(close_r_bone))
    bone_mask = ndimage.binary_fill_holes(bone_mask)

    # Keep only within specimen
    bone_mask = bone_mask & specimen_mask

    bone_vol = np.sum(bone_mask) * voxel_vol_mm3
    print(f"    Bone mask volume after reconstruction: {bone_vol:.0f} mm³")

    # --- Stage 4: Interior protection & boundary refinement ---
    print("  Stage 4: Interior protection & boundary carving...")

    bone_mask = _refine_boundary(volume, bone_mask, grad, spacing,
                                 artifact_weight)

    final_vol = np.sum(bone_mask) * voxel_vol_mm3
    print(f"    Final bone volume: {final_vol:.0f} mm³")

    # --- Stage 5: Separate individual bones ---
    print("  Stage 5: Separating individual bones...")

    labeled, n_components = ndimage.label(bone_mask)
    props = regionprops(labeled, intensity_image=volume)

    bones = []
    small_count = 0

    for prop in props:
        vol_mm3 = prop.area * voxel_vol_mm3
        mean_hu = float(prop.intensity_mean)

        if vol_mm3 >= min_bone_volume_mm3:
            bones.append({
                'label': prop.label,
                'centroid': np.array(prop.centroid) * np.array(spacing),
                'mask': labeled == prop.label,
                'bbox': prop.bbox,
                'volume_mm3': vol_mm3,
                'mean_hu': mean_hu,
            })
        else:
            small_count += 1

    if small_count > 0:
        print(f"    Filtered {small_count} small components "
              f"(< {min_bone_volume_mm3} mm³)")

    # --- Stage 6: Tag-to-bone association ---
    for bone in bones:
        bone['tag_id'] = None
        bone['tag_dist'] = None

    if tag_components and bones:
        bone_centroids = np.array([b['centroid'] for b in bones])
        for tag in tag_components:
            dists = np.linalg.norm(bone_centroids - tag['centroid'], axis=1)
            nearest_idx = int(np.argmin(dists))
            nearest_dist = float(dists[nearest_idx])

            current = bones[nearest_idx].get('tag_dist')
            if current is None or nearest_dist < current:
                bones[nearest_idx]['tag_id'] = tag['label']
                bones[nearest_idx]['tag_dist'] = nearest_dist

    bones.sort(key=lambda b: b['volume_mm3'], reverse=True)

    print(f"\nFound {len(bones)} bones and {len(tag_components)} tags in scan")
    for i, bone in enumerate(bones):
        tag_str = f"tag {bone['tag_id']}" if bone['tag_id'] else "no tag"
        print(f"  Bone {i+1}: {bone['volume_mm3']:.1f} mm³, "
              f"mean HU {bone['mean_hu']:.0f}, {tag_str}")

    return {
        'volume': volume,
        'spacing': spacing,
        'bones': bones,
    }


def _compute_gradient_magnitude(volume, spacing):
    """Compute 3D gradient magnitude of the volume."""
    gz = ndimage.sobel(volume, axis=0) / spacing[0]
    gy = ndimage.sobel(volume, axis=1) / spacing[1]
    gx = ndimage.sobel(volume, axis=2) / spacing[2]
    return np.sqrt(gz**2 + gy**2 + gx**2)


def _morphological_reconstruct(seed, mask):
    """Binary morphological reconstruction: grow seed within mask.

    Equivalent to MATLAB's imreconstruct for binary images.
    Iteratively dilates the seed, intersecting with mask each step,
    until convergence.
    """
    result = seed & mask
    selem = ball(1)

    while True:
        expanded = ndimage.binary_dilation(result, structure=selem) & mask
        if np.array_equal(expanded, result):
            break
        result = expanded

    return result


def _refine_boundary(volume, bone_mask, grad, spacing, artifact_weight):
    """Refine bone boundary: protect interior, carve weak edges.

    Adapted from the scaphoid pipeline's multi-stage boundary refinement:
    - Deep interior (>1mm from surface) is always protected
    - Surface voxels removed only if BOTH low-HU AND low-gradient
    """
    mean_spacing = float(np.mean(spacing))

    # Distance from bone surface into interior
    dist_interior = ndimage.distance_transform_edt(bone_mask, sampling=spacing)
    deep_interior = dist_interior >= 1.0  # mm

    # Boundary band: 0-1mm from surface
    boundary_band = bone_mask & (~deep_interior)

    if not np.any(boundary_band):
        return bone_mask

    # Compute thresholds from boundary band statistics
    band_hu = volume[boundary_band]
    if len(band_hu) == 0:
        return bone_mask

    hu_carve_floor = max(100, float(np.percentile(band_hu, 15)))
    grad_carve_floor = float(np.percentile(grad[boundary_band], 25))

    # Remove boundary voxels that are BOTH low-HU AND low-gradient
    # AND near air (within 1 voxel of HU < -300)
    air_nearby = ndimage.binary_dilation(
        volume < -300, structure=ball(1))

    remove = (boundary_band &
              (volume < hu_carve_floor) &
              (grad < grad_carve_floor) &
              air_nearby)

    # Don't remove voxels near artifacts (could be metal-corrupted, not real air)
    if np.any(artifact_weight > 0.1):
        remove = remove & (artifact_weight < 0.1)

    if np.any(remove):
        bone_mask = bone_mask.copy()
        bone_mask[remove] = False

        # Restore deep interior (never touch it)
        bone_mask = bone_mask | deep_interior

        # Keep largest component
        bone_mask = _keep_largest_component(bone_mask)

        # Fill holes that may have been created
        bone_mask = ndimage.binary_fill_holes(bone_mask)

    return bone_mask


def _keep_largest_component(mask):
    """Keep only the largest connected component."""
    labeled, n = ndimage.label(mask)
    if n <= 1:
        return mask
    sizes = ndimage.sum(mask, labeled, range(1, n + 1))
    largest = int(np.argmax(sizes)) + 1
    return labeled == largest


def _remove_small_components(mask, min_voxels):
    """Remove connected components smaller than min_voxels."""
    labeled, n = ndimage.label(mask)
    if n == 0:
        return mask
    result = mask.copy()
    for i in range(1, n + 1):
        component = labeled == i
        if np.sum(component) < min_voxels:
            result[component] = False
    return result
