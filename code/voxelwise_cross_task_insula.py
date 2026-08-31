#!/usr/bin/env python3
"""Descriptive voxelwise TG--UG correspondence in anatomical left insula.

This analysis is separate from the existing FEAT workflows. It reads the final
subject-level inputs named in the N=132 L3 templates, never reruns FEAT, and
never extracts values from an L3 statistical map.

Stored contrast directions are retained without sign reversal:
  TG:  reciprocated > nonreciprocated
  UGR: positive fairness parametric modulation
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


# Loaded after argument parsing so that --help works even before dependencies
# have been installed on a new analysis host.
np: Any = None
nib: Any = None
resample_from_to: Any = None
odr: Any = None


def load_dependencies() -> None:
    global np, nib, resample_from_to, odr
    try:
        import numpy as numpy_module
        import nibabel as nibabel_module
        from nibabel.processing import resample_from_to as resample_function
        from scipy import odr as odr_module
        import matplotlib  # noqa: F401 - checked here, used later for plotting
    except ImportError as exc:  # pragma: no cover - depends on analysis host
        raise RuntimeError(
            "Missing Python dependency. Install numpy, nibabel, scipy, and "
            f"matplotlib. Original error: {exc}"
        ) from exc
    np = numpy_module
    nib = nibabel_module
    resample_from_to = resample_function
    odr = odr_module


TG_EXPECTED_COPE = 10
UG_EXPECTED_COPE = 11
GRID_ATOL = 1e-5


@dataclass(frozen=True)
class SubjectRoute:
    subject: str
    template_path: str
    cope: Path
    zstat: Path
    mask: Path

    @property
    def level(self) -> str:
        return "L1" if "/L1_" in self.template_path else "L2"

    @property
    def run(self) -> str:
        match = re.search(r"_run-(\d+)_", self.template_path)
        return match.group(1) if match else "both"


class RunningStats:
    """Vectorized Welford summaries without a subject-by-voxel output table."""

    def __init__(self, size: int) -> None:
        self.n = 0
        self.mean = np.zeros(size, dtype=np.float64)
        self.m2 = np.zeros(size, dtype=np.float64)

    def update(self, values: np.ndarray) -> None:
        self.n += 1
        delta = values - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (values - self.mean)

    def finish(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.n < 2:
            raise RuntimeError("At least two subjects are required.")
        sd = np.sqrt(np.maximum(self.m2 / (self.n - 1), 0))
        return self.mean, sd, sd / np.sqrt(self.n)


def parse_args() -> argparse.Namespace:
    script = Path(__file__).resolve()
    project_root = script.parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Create COPE and mean-subject-ZSTAT voxelwise TG--UG left-insula "
            "correspondence tables, ODR plots, and signed residual maps."
        )
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=project_root,
        help=f"rf1-betrayal root (default: {project_root})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Statistical output directory (default: PROJECT_ROOT/code/"
            "voxelwise_cross_task_insula_output)"
        ),
    )
    parser.add_argument(
        "--fsldir", type=Path, help="FSL installation root (default: $FSLDIR)"
    )
    parser.add_argument(
        "--atlas-file",
        type=Path,
        help="Optional explicit HarvardOxford-cort-maxprob-thr25 NIfTI",
    )
    parser.add_argument(
        "--atlas-xml",
        type=Path,
        help="Optional explicit Harvard-Oxford cortical atlas XML",
    )
    return parser.parse_args()


def read_subjects(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Subject list not found: {path}")
    subjects = [
        line.strip().removeprefix("sub-")
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]
    if not subjects:
        raise ValueError(f"Subject list is empty: {path}")
    if len(subjects) != len(set(subjects)):
        raise ValueError(f"Subject list contains duplicate IDs: {path}")
    return subjects


def parse_fsf_names(path: Path) -> dict[int, str]:
    pattern = re.compile(
        r'^set fmri\(conname_real\.(\d+)\) "([^"]+)"$', re.MULTILINE
    )
    return {
        int(number): name
        for number, name in pattern.findall(path.read_text(encoding="utf-8"))
    }


def parse_ev_titles(path: Path) -> dict[int, str]:
    pattern = re.compile(r'^set fmri\(evtitle(\d+)\) "([^"]+)"$', re.MULTILINE)
    return {
        int(number): name
        for number, name in pattern.findall(path.read_text(encoding="utf-8"))
    }


def parse_contrast_weights(path: Path, cope: int) -> dict[int, float]:
    pattern = re.compile(
        rf"^set fmri\(con_real{cope}\.(\d+)\)\s+([-+0-9.eE]+)$",
        re.MULTILINE,
    )
    return {
        int(ev): float(weight)
        for ev, weight in pattern.findall(path.read_text(encoding="utf-8"))
    }


def find_named_contrast(names: dict[int, str], exact_name: str, fsf: Path) -> int:
    matches = [number for number, name in names.items() if name == exact_name]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one contrast named {exact_name!r} in {fsf}; found {matches}."
        )
    return matches[0]


def validate_contrasts(project_root: Path) -> dict[str, Any]:
    tg_fsf = project_root / "templates" / "L1_task-trust_model-01_type-act.fsf"
    ug_fsf = project_root / "templates" / "L1_task-ugr_model-3_type-act.fsf"
    ug_generator = project_root / "code" / "a4_model-3.py"
    for path in (tg_fsf, ug_fsf, ug_generator):
        if not path.is_file():
            raise FileNotFoundError(f"Required contrast-definition file not found: {path}")

    tg_cope = find_named_contrast(parse_fsf_names(tg_fsf), "rec-def", tg_fsf)
    ug_cope = find_named_contrast(
        parse_fsf_names(ug_fsf), "offer (un)fairness (pmod)", ug_fsf
    )
    if tg_cope != TG_EXPECTED_COPE or ug_cope != UG_EXPECTED_COPE:
        raise ValueError(
            "Contrast numbering changed: expected TG cope 10 and UGR cope 11, "
            f"but found TG cope {tg_cope} and UGR cope {ug_cope}."
        )

    tg_titles = parse_ev_titles(tg_fsf)
    tg_weights = parse_contrast_weights(tg_fsf, tg_cope)
    defect_evs = sorted(ev for ev, title in tg_titles.items() if title.endswith("_def"))
    recip_evs = sorted(ev for ev, title in tg_titles.items() if title.endswith("_rec"))
    if defect_evs != [4, 6, 8] or recip_evs != [5, 7, 9]:
        raise ValueError(f"Unexpected TG outcome EV definitions in {tg_fsf}.")
    if not all(np.isclose(tg_weights.get(ev), -1.0) for ev in defect_evs) or not all(
        np.isclose(tg_weights.get(ev), 1.0) for ev in recip_evs
    ):
        raise ValueError(
            "TG cope 10 is no longer reciprocated > nonreciprocated; review the analysis."
        )

    ug_titles = parse_ev_titles(ug_fsf)
    ug_weights = parse_contrast_weights(ug_fsf, ug_cope)
    offer_pmods = [
        ev
        for ev, title in ug_titles.items()
        if title.endswith("_pmod") and not title.startswith("rt_")
    ]
    if offer_pmods != [2, 4, 6, 8] or not all(
        np.isclose(ug_weights.get(ev), 1.0) for ev in offer_pmods
    ):
        raise ValueError("UGR cope 11 no longer positively weights all offer pmods.")
    generator_text = ug_generator.read_text(encoding="utf-8")
    required_fragments = (
        'df["offer_amount"] = df["L_Option"] + df["R_Option"]',
        'group["demeaned_offer"] = group["offer_amount"] - group["offer_amount"].mean()',
        '"demeaned_offer"',
    )
    if not all(fragment in generator_text for fragment in required_fragments):
        raise ValueError(
            f"Could not verify positive demeaned-offer pmods in {ug_generator}."
        )

    return {
        "tg_cope": tg_cope,
        "tg_name": "rec-def",
        "tg_stored_direction": "reciprocated > nonreciprocated",
        "tg_applied_multiplier": 1.0,
        "tg_l1_fsf": str(tg_fsf),
        "ug_cope": ug_cope,
        "ug_name": "offer (un)fairness (pmod)",
        "ug_stored_direction": "fairness (positive demeaned offer amount)",
        "ug_applied_multiplier": 1.0,
        "ug_l1_fsf": str(ug_fsf),
        "ug_pmod_generator": str(ug_generator),
    }


def remap_template_path(template_path: str, project_root: Path) -> Path:
    normalized = template_path.replace("\\", "/")
    marker = "/derivatives/"
    if marker not in normalized:
        raise ValueError(f"L3 input is not under a derivatives directory: {template_path}")
    relative = normalized.split(marker, 1)[1]
    return project_root / "derivatives" / Path(relative)


def zstat_path_for(cope_path: Path) -> Path:
    match = re.fullmatch(r"cope(\d+)\.nii\.gz", cope_path.name)
    if not match:
        raise ValueError(f"Unexpected subject-level COPE filename: {cope_path}")
    return cope_path.with_name(f"zstat{match.group(1)}.nii.gz")


def parse_l3_routes(template: Path, cope: int, project_root: Path) -> list[SubjectRoute]:
    if not template.is_file():
        raise FileNotFoundError(f"L3 template not found: {template}")
    pattern = re.compile(r'^set feat_files\((\d+)\) "([^"]+)"$', re.MULTILINE)
    matches = pattern.findall(template.read_text(encoding="utf-8"))
    if not matches:
        raise ValueError(f"No feat_files entries found in {template}")

    routes: list[SubjectRoute] = []
    for expected_index, (raw_index, raw_path) in enumerate(matches, start=1):
        if int(raw_index) != expected_index:
            raise ValueError(
                f"Non-contiguous feat_files entries in {template} at {raw_index}."
            )
        routed_path = raw_path.replace("COPENUM", str(cope))
        subject_match = re.search(r"/sub-([^/]+)/", routed_path)
        if not subject_match:
            raise ValueError(f"Cannot parse subject ID from L3 input: {raw_path}")
        cope_path = remap_template_path(routed_path, project_root)
        routes.append(
            SubjectRoute(
                subject=subject_match.group(1),
                template_path=routed_path,
                cope=cope_path,
                zstat=zstat_path_for(cope_path),
                mask=cope_path.parent.parent / "mask.nii.gz",
            )
        )
    return routes


def assert_route_agreement(
    subjects: list[str],
    primary: list[SubjectRoute],
    alternate: list[SubjectRoute],
    task: str,
) -> None:
    if [route.subject for route in primary] != subjects:
        raise ValueError(
            f"{task} L3 template order does not exactly match sublist_n132.txt."
        )
    primary_paths = [(route.subject, route.template_path) for route in primary]
    alternate_paths = [(route.subject, route.template_path) for route in alternate]
    if primary_paths != alternate_paths:
        raise ValueError(f"{task} full and ones N=132 L3 templates use different inputs.")


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def require_complete_sample(
    subjects: list[str],
    tg_routes: list[SubjectRoute],
    ug_routes: list[SubjectRoute],
    output_dir: Path,
) -> list[tuple[SubjectRoute, SubjectRoute]]:
    complete: list[tuple[SubjectRoute, SubjectRoute]] = []
    excluded: list[dict[str, str]] = []
    for subject, tg, ug in zip(subjects, tg_routes, ug_routes, strict=True):
        required = {
            "TG COPE": tg.cope,
            "TG ZSTAT": tg.zstat,
            "TG FEAT mask": tg.mask,
            "UGR COPE": ug.cope,
            "UGR ZSTAT": ug.zstat,
            "UGR FEAT mask": ug.mask,
        }
        missing = [f"{label} missing: {path}" for label, path in required.items() if not path.is_file()]
        if missing:
            reason = "; ".join(missing)
            excluded.append({"subject": subject, "reason": reason})
            print(f"MISSING sub-{subject}: {reason}", file=sys.stderr)
        else:
            complete.append((tg, ug))

    write_tsv(output_dir / "excluded_subjects.tsv", excluded, ["subject", "reason"])
    included_rows = [
        {
            "subject": tg.subject,
            "tg_level": tg.level,
            "tg_run": tg.run,
            "tg_cope": str(tg.cope),
            "tg_zstat": str(tg.zstat),
            "tg_mask": str(tg.mask),
            "ug_level": ug.level,
            "ug_run": ug.run,
            "ug_cope": str(ug.cope),
            "ug_zstat": str(ug.zstat),
            "ug_mask": str(ug.mask),
        }
        for tg, ug in complete
    ]
    write_tsv(
        output_dir / "included_subjects.tsv",
        included_rows,
        [
            "subject",
            "tg_level",
            "tg_run",
            "tg_cope",
            "tg_zstat",
            "tg_mask",
            "ug_level",
            "ug_run",
            "ug_cope",
            "ug_zstat",
            "ug_mask",
        ],
    )
    if excluded:
        raise RuntimeError(
            f"{len(excluded)} of the required 132 subjects lack inputs. No analysis "
            f"was run and the sample was not reduced; see {output_dir / 'excluded_subjects.tsv'}"
        )
    if len(complete) != 132:
        raise RuntimeError(
            f"Internal sample validation failed: expected 132 complete subjects, found {len(complete)}."
        )
    return complete


def load_3d(path: Path) -> nib.spatialimages.SpatialImage:
    image = nib.load(str(path))
    if len(image.shape) != 3:
        raise ValueError(f"Expected a 3D image, found shape {image.shape}: {path}")
    return image


def same_grid(a: nib.spatialimages.SpatialImage, b: nib.spatialimages.SpatialImage) -> bool:
    return a.shape == b.shape and np.allclose(
        a.affine, b.affine, atol=GRID_ATOL, rtol=0
    )


def validate_subject_grids(
    included: list[tuple[SubjectRoute, SubjectRoute]],
) -> nib.spatialimages.SpatialImage:
    reference_path = included[0][0].cope
    reference = load_3d(reference_path)
    for tg, ug in included:
        for label, path in (
            ("TG COPE", tg.cope),
            ("TG ZSTAT", tg.zstat),
            ("TG mask", tg.mask),
            ("UGR COPE", ug.cope),
            ("UGR ZSTAT", ug.zstat),
            ("UGR mask", ug.mask),
        ):
            image = load_3d(path)
            if not same_grid(reference, image):
                raise ValueError(
                    f"Grid mismatch for sub-{tg.subject} {label}: {path}\n"
                    f"Reference {reference_path}: shape={reference.shape}, affine=\n"
                    f"{reference.affine}\nImage: shape={image.shape}, affine=\n"
                    f"{image.affine}\nSubject statistical images are never resampled."
                )
    return reference


def discover_xml(fsldir: Path, explicit: Path | None) -> Path:
    if explicit:
        if not explicit.is_file():
            raise FileNotFoundError(f"Atlas XML not found: {explicit}")
        return explicit.resolve()
    atlas_dir = fsldir / "data" / "atlases"
    preferred = atlas_dir / "HarvardOxford-Cortical.xml"
    if preferred.is_file():
        return preferred.resolve()
    candidates = sorted(atlas_dir.rglob("*HarvardOxford*Cortical*.xml"))
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Could not uniquely locate Harvard-Oxford cortical XML under {atlas_dir}; "
            f"found {[str(path) for path in candidates]}"
        )
    return candidates[0].resolve()


def atlas_candidates(fsldir: Path, explicit: Path | None) -> list[Path]:
    if explicit:
        if not explicit.is_file():
            raise FileNotFoundError(f"Atlas NIfTI not found: {explicit}")
        if "HarvardOxford-cort-maxprob-thr25-" not in explicit.name:
            raise ValueError(
                "--atlas-file must be a Harvard-Oxford cortical maxprob-thr25 atlas."
            )
        return [explicit.resolve()]
    atlas_dir = fsldir / "data" / "atlases"
    candidates = sorted(
        path.resolve()
        for path in atlas_dir.rglob("HarvardOxford-cort-maxprob-thr25-*.nii*")
        if path.is_file()
    )
    if not candidates:
        raise FileNotFoundError(
            f"No HarvardOxford-cort-maxprob-thr25 atlas found under {atlas_dir}"
        )
    return candidates


def choose_atlas(
    candidates: list[Path], reference: nib.spatialimages.SpatialImage
) -> tuple[Path, nib.spatialimages.SpatialImage, str]:
    loaded = [(path, load_3d(path)) for path in candidates]
    exact = [(path, image) for path, image in loaded if same_grid(image, reference)]
    if exact:
        path, image = sorted(exact, key=lambda item: str(item[0]))[0]
        return path, image, "exact statistical-grid match"
    path, image = min(
        loaded,
        key=lambda item: (
            abs(np.linalg.det(item[1].affine[:3, :3])),
            str(item[0]),
        ),
    )
    return (
        path,
        image,
        "no exact grid match; discovered atlas with smallest voxel volume",
    )


def insula_label_from_xml(xml_path: Path) -> tuple[int, int, str]:
    root = ET.parse(xml_path).getroot()
    matches: list[tuple[int, str]] = []
    for label in root.findall(".//label"):
        name = " ".join((label.text or "").split())
        if name.casefold() == "insular cortex":
            if "index" not in label.attrib:
                raise ValueError(f"Insular Cortex label lacks an index in {xml_path}")
            matches.append((int(label.attrib["index"]), name))
    if len(matches) != 1:
        raise ValueError(
            f"Expected one XML label named 'Insular Cortex' in {xml_path}; found {matches}"
        )
    xml_index, name = matches[0]
    # FSL XML indices are zero-based; value 0 in a max-probability atlas is background.
    return xml_index, xml_index + 1, name


def make_left_insula_mask(
    atlas: nib.spatialimages.SpatialImage,
    label_value: int,
    reference: nib.spatialimages.SpatialImage,
) -> tuple[np.ndarray, int, int, bool]:
    atlas_data = np.asanyarray(atlas.dataobj)
    assigned = np.isclose(atlas_data, label_value)
    indices = np.argwhere(assigned)
    left = np.zeros(atlas.shape, dtype=bool)
    if indices.size:
        xyz = nib.affines.apply_affine(atlas.affine, indices)
        left_indices = indices[xyz[:, 0] < 0]
        left[tuple(left_indices.T)] = True
    native_count = int(left.sum())
    if native_count == 0:
        raise ValueError(
            f"Insular Cortex atlas label {label_value} produced an empty left mask."
        )

    was_resampled = not same_grid(atlas, reference)
    if was_resampled:
        source = nib.Nifti1Image(left.astype(np.uint8), atlas.affine)
        resampled = resample_from_to(
            source,
            (reference.shape, reference.affine),
            order=0,
            mode="constant",
            cval=0,
        )
        left_reference = np.asanyarray(resampled.dataobj) > 0.5
    else:
        left_reference = left
    reference_count = int(left_reference.sum())
    if reference_count == 0:
        raise ValueError("Left-insula mask is empty on the statistical reference grid.")
    return left_reference, native_count, reference_count, was_resampled


def coverage_intersection(
    included: list[tuple[SubjectRoute, SubjectRoute]],
    reference_shape: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    tg_coverage = np.ones(reference_shape, dtype=bool)
    ug_coverage = np.ones(reference_shape, dtype=bool)
    for tg, ug in included:
        tg_mask = np.asanyarray(nib.load(str(tg.mask)).dataobj)
        ug_mask = np.asanyarray(nib.load(str(ug.mask)).dataobj)
        tg_coverage &= np.isfinite(tg_mask) & (tg_mask > 0)
        ug_coverage &= np.isfinite(ug_mask) & (ug_mask > 0)
    return tg_coverage, ug_coverage


def save_nifti(
    data: np.ndarray,
    path: Path,
    reference: nib.spatialimages.SpatialImage,
    dtype: np.dtype,
) -> None:
    header = reference.header.copy()
    header.set_data_dtype(dtype)
    image = nib.Nifti1Image(data.astype(dtype), reference.affine, header=header)
    image.set_qform(reference.affine, int(reference.header["qform_code"]))
    image.set_sform(reference.affine, int(reference.header["sform_code"]))
    nib.save(image, str(path))


def extract_summaries(
    included: list[tuple[SubjectRoute, SubjectRoute]], voxel_indices: np.ndarray
) -> dict[str, dict[str, np.ndarray]]:
    selection = tuple(voxel_indices.T)
    accumulators = {
        "cope_tg": RunningStats(len(voxel_indices)),
        "cope_ug": RunningStats(len(voxel_indices)),
        "zstat_tg": RunningStats(len(voxel_indices)),
        "zstat_ug": RunningStats(len(voxel_indices)),
    }
    for tg, ug in included:
        # Retain both stored contrast directions; no multiplication by -1.
        arrays = {
            "cope_tg": np.asanyarray(nib.load(str(tg.cope)).dataobj)[selection],
            "cope_ug": np.asanyarray(nib.load(str(ug.cope)).dataobj)[selection],
            "zstat_tg": np.asanyarray(nib.load(str(tg.zstat)).dataobj)[selection],
            "zstat_ug": np.asanyarray(nib.load(str(ug.zstat)).dataobj)[selection],
        }
        for name, values in arrays.items():
            values = np.asarray(values, dtype=np.float64)
            if not np.all(np.isfinite(values)):
                bad = int((~np.isfinite(values)).sum())
                raise ValueError(
                    f"sub-{tg.subject} has {bad} non-finite {name} values in the ROI."
                )
            accumulators[name].update(values)

    result: dict[str, dict[str, np.ndarray]] = {}
    for quantity in ("cope", "zstat"):
        tg_mean, tg_sd, tg_sem = accumulators[f"{quantity}_tg"].finish()
        ug_mean, ug_sd, ug_sem = accumulators[f"{quantity}_ug"].finish()
        result[quantity] = {
            "tg_mean": tg_mean,
            "tg_sd": tg_sd,
            "tg_sem": tg_sem,
            "ug_mean": ug_mean,
            "ug_sd": ug_sd,
            "ug_sem": ug_sem,
        }
    return result


def fit_odr(x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    if np.isclose(np.std(x), 0) or np.isclose(np.std(y), 0):
        raise ValueError("Cannot fit ODR because a voxelwise mean has zero variance.")
    initial_slope, initial_intercept = np.polyfit(x, y, 1)
    model = odr.Model(lambda beta, value: beta[0] + beta[1] * value)
    fit = odr.ODR(
        odr.RealData(x, y),
        model,
        beta0=[initial_intercept, initial_slope],
    ).run()
    intercept, slope = float(fit.beta[0]), float(fit.beta[1])
    residual = (y - (intercept + slope * x)) / np.sqrt(1.0 + slope**2)
    return {
        "pearson_r": float(np.corrcoef(x, y)[0, 1]),
        "odr_intercept": intercept,
        "odr_slope": slope,
        "odr_intercept_sd": float(fit.sd_beta[0]),
        "odr_slope_sd": float(fit.sd_beta[1]),
        "signed_orthogonal_residual": residual,
        "odr_stop_reason": list(fit.stopreason),
        "method": "unweighted orthogonal distance regression",
    }


def write_voxel_table(
    path: Path,
    quantity: str,
    indices: np.ndarray,
    coordinates: np.ndarray,
    summary: dict[str, np.ndarray],
    fit: dict[str, Any],
) -> None:
    prefix = "cope" if quantity == "cope" else "subject_zstat"
    fields = [
        "i",
        "j",
        "k",
        "mni_x_mm",
        "mni_y_mm",
        "mni_z_mm",
        f"tg_mean_{prefix}_recip_gt_nonrecip",
        f"tg_sd_{prefix}_recip_gt_nonrecip",
        f"tg_sem_{prefix}_recip_gt_nonrecip",
        f"ug_mean_{prefix}_fairness_pmod",
        f"ug_sd_{prefix}_fairness_pmod",
        f"ug_sem_{prefix}_fairness_pmod",
        "signed_odr_residual",
    ]
    residual = fit["signed_orthogonal_residual"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for position in range(len(indices)):
            writer.writerow(
                [
                    *[int(value) for value in indices[position]],
                    *[f"{value:.6f}" for value in coordinates[position]],
                    f"{summary['tg_mean'][position]:.10g}",
                    f"{summary['tg_sd'][position]:.10g}",
                    f"{summary['tg_sem'][position]:.10g}",
                    f"{summary['ug_mean'][position]:.10g}",
                    f"{summary['ug_sd'][position]:.10g}",
                    f"{summary['ug_sem'][position]:.10g}",
                    f"{residual[position]:.10g}",
                ]
            )


def make_plot(
    output_stem: Path,
    quantity: str,
    x: np.ndarray,
    y: np.ndarray,
    fit: dict[str, Any],
    n_subjects: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    quantity_label = (
        "mean COPE" if quantity == "cope" else "mean subject-level Z-statistic"
    )
    fig, ax = plt.subplots(figsize=(7.1, 6.0), constrained_layout=True)
    ax.scatter(
        x,
        y,
        s=18,
        alpha=0.58,
        color="#276FBF",
        edgecolors="white",
        linewidths=0.25,
        rasterized=True,
    )
    x_line = np.linspace(float(np.min(x)), float(np.max(x)), 250)
    ax.plot(
        x_line,
        fit["odr_intercept"] + fit["odr_slope"] * x_line,
        color="#B23A48",
        linewidth=2.2,
        label="ODR fit",
    )
    ax.set_xlabel(
        f"Trust Game: reciprocated > nonreciprocated\n({quantity_label})",
        fontsize=11,
    )
    ax.set_ylabel(
        f"Ultimatum Game: fairness parametric modulation\n({quantity_label})",
        fontsize=11,
    )
    ax.set_title(
        "Voxelwise cross-task correspondence in anatomical left insula",
        fontsize=12.5,
    )
    ax.text(
        0.03,
        0.97,
        (
            f"N = {n_subjects} subjects; V = {len(x)} voxels\n"
            f"Pearson r = {fit['pearson_r']:.3f}\n"
            f"ODR: y = {fit['odr_intercept']:.3g} + {fit['odr_slope']:.3g}x"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "alpha": 0.9,
            "edgecolor": "0.75",
        },
    )
    ax.text(
        0.5,
        -0.18,
        "Descriptive voxelwise correspondence; spatial dependence precludes treating voxels as independent.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.2,
        color="0.35",
    )
    ax.grid(True, color="0.9", linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    load_dependencies()
    project_root = args.project_root.expanduser().resolve()
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else project_root / "code" / "voxelwise_cross_task_insula_output"
    )
    masks_dir = project_root / "masks"
    output_dir.mkdir(parents=True, exist_ok=True)
    if not masks_dir.is_dir():
        raise FileNotFoundError(f"Repository masks directory not found: {masks_dir}")

    fsldir_value = args.fsldir or (
        Path(os.environ["FSLDIR"]) if os.environ.get("FSLDIR") else None
    )
    if fsldir_value is None:
        raise RuntimeError("FSLDIR is not set. Export FSLDIR or pass --fsldir.")
    fsldir = fsldir_value.expanduser().resolve()
    if not fsldir.is_dir():
        raise FileNotFoundError(f"FSL installation directory not found: {fsldir}")

    print("Validating contrast definitions and L3 subject-level routes...")
    contrast = validate_contrasts(project_root)
    subjects_path = project_root / "code" / "sublist_n132.txt"
    subjects = read_subjects(subjects_path)
    if len(subjects) != 132:
        raise ValueError(
            f"Expected sublist_n132.txt to contain 132 subjects; found {len(subjects)}."
        )

    templates = project_root / "templates"
    template_paths = {
        "tg_full": templates
        / "L3_task-trust_model-01_type-act_group-full_n132_flame1.fsf",
        "tg_ones": templates
        / "L3_task-trust_model-01_type-act_group-ones_n132_flame1.fsf",
        "ug_full": templates
        / "L3_task-ugr_model-3_type-act_group-full_n132_flame1.fsf",
        "ug_ones": templates
        / "L3_task-ugr_model-3_type-act_group-ones_n132_flame1.fsf",
    }
    tg_routes = parse_l3_routes(
        template_paths["tg_full"], contrast["tg_cope"], project_root
    )
    ug_routes = parse_l3_routes(
        template_paths["ug_full"], contrast["ug_cope"], project_root
    )
    tg_ones = parse_l3_routes(
        template_paths["tg_ones"], contrast["tg_cope"], project_root
    )
    ug_ones = parse_l3_routes(
        template_paths["ug_ones"], contrast["ug_cope"], project_root
    )
    assert_route_agreement(subjects, tg_routes, tg_ones, "TG")
    assert_route_agreement(subjects, ug_routes, ug_ones, "UGR")

    included = require_complete_sample(subjects, tg_routes, ug_routes, output_dir)
    print("Final subject N: 132 (all subjects in sublist_n132.txt are required)")
    reference = validate_subject_grids(included)
    reference_path = included[0][0].cope
    print(f"Statistical reference: {reference_path}")
    print(f"Reference shape: {reference.shape}")
    print("Reference affine:\n" + np.array2string(reference.affine, precision=6))

    xml_path = discover_xml(fsldir, args.atlas_xml)
    atlas_path, atlas, atlas_rule = choose_atlas(
        atlas_candidates(fsldir, args.atlas_file), reference
    )
    xml_index, label_value, label_name = insula_label_from_xml(xml_path)
    left_insula, native_voxels, reference_voxels, was_resampled = (
        make_left_insula_mask(atlas, label_value, reference)
    )
    print(f"Harvard-Oxford atlas: {atlas_path}")
    print(f"Atlas selection: {atlas_rule}")
    print(f"Atlas XML: {xml_path}")
    print(
        f"Atlas XML label: {label_name!r}, XML index={xml_index}, "
        f"NIfTI integer={label_value}"
    )
    print(f"Initial whole-left-insula voxels (atlas grid): {native_voxels}")
    if was_resampled:
        print(
            "Whole-left-insula voxels after nearest-neighbor mask resampling: "
            f"{reference_voxels}"
        )

    tg_coverage, ug_coverage = coverage_intersection(included, reference.shape)
    final_mask = left_insula & tg_coverage & ug_coverage
    final_voxels = int(final_mask.sum())
    if final_voxels == 0:
        raise RuntimeError(
            "The anatomical left-insula x TG coverage x UGR coverage mask is empty."
        )
    tg_covered = int((left_insula & tg_coverage).sum())
    ug_covered = int((left_insula & ug_coverage).sum())
    print(f"Left-insula voxels with all-subject TG coverage: {tg_covered}")
    print(f"Left-insula voxels with all-subject UGR coverage: {ug_covered}")
    print(f"Final dual-task coverage voxel count: {final_voxels}")

    anatomical_path = (
        masks_dir / "left_insula_maxprob-thr25_anatomical_refgrid.nii.gz"
    )
    final_mask_path = (
        masks_dir / "left_insula_maxprob-thr25_dualtask_coverage_mask.nii.gz"
    )
    save_nifti(left_insula, anatomical_path, reference, np.uint8)
    save_nifti(final_mask, final_mask_path, reference, np.uint8)

    indices = np.argwhere(final_mask)
    coordinates = nib.affines.apply_affine(reference.affine, indices)
    summaries = extract_summaries(included, indices)
    fit_results: dict[str, dict[str, Any]] = {}
    output_files: dict[str, str] = {
        "anatomical_mask_reference_grid": str(anatomical_path),
        "final_dualtask_coverage_mask": str(final_mask_path),
        "included_subjects": str(output_dir / "included_subjects.tsv"),
        "excluded_subjects": str(output_dir / "excluded_subjects.tsv"),
    }
    for quantity in ("cope", "zstat"):
        summary = summaries[quantity]
        fit = fit_odr(summary["tg_mean"], summary["ug_mean"])
        fit_results[quantity] = fit
        table_path = output_dir / f"voxelwise_cross_task_insula_{quantity}.tsv"
        plot_stem = output_dir / f"voxelwise_cross_task_insula_{quantity}_scatter"
        residual_path = (
            output_dir
            / f"voxelwise_cross_task_insula_{quantity}_signed_odr_residual.nii.gz"
        )
        write_voxel_table(table_path, quantity, indices, coordinates, summary, fit)
        make_plot(
            plot_stem,
            quantity,
            summary["tg_mean"],
            summary["ug_mean"],
            fit,
            len(included),
        )
        residual_map = np.zeros(reference.shape, dtype=np.float32)
        residual_map[tuple(indices.T)] = fit["signed_orthogonal_residual"]
        save_nifti(residual_map, residual_path, reference, np.float32)
        output_files.update(
            {
                f"{quantity}_voxel_table": str(table_path),
                f"{quantity}_scatter_pdf": str(plot_stem.with_suffix(".pdf")),
                f"{quantity}_scatter_png": str(plot_stem.with_suffix(".png")),
                f"{quantity}_signed_odr_residual_map": str(residual_path),
            }
        )
        label = "COPE" if quantity == "cope" else "mean subject-level ZSTAT"
        print(
            f"{label}: Pearson r={fit['pearson_r']:.6f}; "
            f"ODR slope={fit['odr_slope']:.6f}; "
            f"intercept={fit['odr_intercept']:.6f}"
        )

    correspondence = {
        quantity: {
            key: value
            for key, value in fit.items()
            if key != "signed_orthogonal_residual"
        }
        for quantity, fit in fit_results.items()
    }
    metadata_path = output_dir / "analysis_metadata.json"
    output_files["metadata"] = str(metadata_path)
    metadata = {
        "analysis": "descriptive voxelwise cross-task correspondence in anatomical left insula",
        "implementation": "Python",
        "package_versions": {
            "numpy": np.__version__,
            "nibabel": nib.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": __import__("matplotlib").__version__,
        },
        "project_root": str(project_root),
        "subject_list": str(subjects_path),
        "starting_n": len(subjects),
        "included_n": len(included),
        "excluded_n": 0,
        "contrast_validation": contrast,
        "l3_templates": {key: str(value) for key, value in template_paths.items()},
        "reference_image": str(reference_path),
        "reference_shape": list(reference.shape),
        "reference_affine": reference.affine.tolist(),
        "fsl_dir": str(fsldir),
        "harvard_oxford": {
            "atlas_file": str(atlas_path),
            "atlas_selection_rule": atlas_rule,
            "xml_file": str(xml_path),
            "xml_label": label_name,
            "xml_zero_based_index": xml_index,
            "nifti_integer_value": label_value,
            "left_hemisphere_rule": "voxel center MNI x < 0 mm",
            "initial_left_insula_voxels_atlas_grid": native_voxels,
            "left_insula_voxels_reference_grid": reference_voxels,
            "mask_resampled_nearest_neighbor": was_resampled,
        },
        "coverage": {
            "criterion": "intersection of every required subject's final TG and UGR FEAT masks",
            "left_insula_with_all_subject_tg_coverage": tg_covered,
            "left_insula_with_all_subject_ugr_coverage": ug_covered,
            "final_voxels": final_voxels,
        },
        "voxel_index_convention": (
            "i/j/k are zero-based NIfTI indices; MNI coordinates are millimetres"
        ),
        "correspondence": correspondence,
        "signed_residual_definition": (
            "(UG - (ODR intercept + ODR slope * TG)) / sqrt(1 + slope^2); "
            "positive means a stronger UG/fairness response than predicted from the "
            "TG reciprocity response, while negative means a weaker UG/fairness "
            "response than predicted. Residual maps are zero outside the saved final "
            "coverage mask."
        ),
        "inference_note": (
            "Pearson r and ODR are descriptive. No conventional voxelwise Pearson "
            "p-value is computed because neighboring voxels are spatially dependent."
        ),
        "outputs": output_files,
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Statistical outputs written to: {output_dir}")
    print(f"Masks written to: {masks_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
