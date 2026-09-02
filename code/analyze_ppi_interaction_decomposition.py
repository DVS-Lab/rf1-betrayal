#!/usr/bin/env python3
"""Decompose the UGR partner-by-betrayal aIns-vmPFC gPPI interaction.

This exploratory script uses contrasts already defined by the UGR gPPI model:

    COPE 13: nonsocial fairness-related PPI modulation
    COPE 14: social fairness-related PPI modulation
    COPE 12: social - nonsocial fairness-related PPI modulation

It does not create or estimate new first-level contrasts. By default, it
parses the exact subject order and covariates from the original Level 3 FSF,
calls FSL's ``fslmeants`` to extract ROI-average COPE values, fits a paired
mixed-effects model, and creates two exploratory figures.

Required Python packages: numpy, pandas, patsy, statsmodels, matplotlib.
Required for image extraction: FSL's fslmeants on PATH.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Any


DEFAULT_FSF = Path(
    "templates/L3_task-ugr_model-3_type-ppi_group-AIns_n132_flame1.fsf"
)
DEFAULT_OUT_DIR = Path("exploratory/reviewer1_comment2_ppi_decomposition")
COPE_NUMBERS = {
    "difference_cope": 12,
    "nonsocial_cope": 13,
    "social_cope": 14,
}
NUISANCE_COVARIATES = ["fd_mean", "tsnr", "age", "male", "female"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Decompose the existing social-minus-nonsocial fairness-related "
            "aIns-vmPFC PPI contrast into its social and nonsocial components."
        )
    )
    parser.add_argument(
        "--fsf",
        type=Path,
        default=DEFAULT_FSF,
        help=(
            "Level 3 FSF containing feat_files(), evtitle(), and evg() entries "
            f"(default: {DEFAULT_FSF})."
        ),
    )
    parser.add_argument(
        "--roi",
        type=Path,
        help=(
            "vmPFC ROI used by fslmeants. Required unless --effects-csv is "
            "provided. An ROI selected from the original interaction is "
            "descriptive; use an independent ROI for inferential tests."
        ),
    )
    parser.add_argument(
        "--effects-csv",
        type=Path,
        help=(
            "Optional pre-extracted CSV with sub_id, nonsocial_cope, "
            "social_cope, and optionally difference_cope. The FSF is still "
            "used for the original covariates."
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Exploratory output directory (default: {DEFAULT_OUT_DIR}).",
    )
    return parser.parse_args()


def load_dependencies() -> dict[str, Any]:
    missing: list[str] = []
    modules: dict[str, Any] = {}
    for name in ["numpy", "pandas", "patsy", "statsmodels.formula.api"]:
        try:
            modules[name] = __import__(name, fromlist=["*"])
        except ImportError:
            missing.append(name.split(".")[0])

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        modules["matplotlib.pyplot"] = plt
    except ImportError:
        missing.append("matplotlib")

    if missing:
        raise SystemExit(
            "Missing Python packages: "
            + ", ".join(sorted(set(missing)))
            + ". Install them in the Linux analysis environment before running."
        )
    return modules


def parse_group_fsf(path: Path, pd: Any) -> Any:
    """Return subject paths and the exact group covariates encoded in an FSF."""
    if not path.is_file():
        raise FileNotFoundError(f"Level 3 FSF not found: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    feat_re = re.compile(r'^set feat_files\(([0-9]+)\) "([^"]+)"$')
    title_re = re.compile(r'^set fmri\(evtitle([0-9]+)\) "([^"]+)"$')
    ev_re = re.compile(r"^set fmri\(evg([0-9]+)\.([0-9]+)\) ([^\s]+)$")

    feat_paths: dict[int, str] = {}
    titles: dict[int, str] = {}
    ev_values: dict[tuple[int, int], float] = {}

    for line in lines:
        if match := feat_re.match(line):
            feat_paths[int(match.group(1))] = match.group(2)
        elif match := title_re.match(line):
            titles[int(match.group(1))] = match.group(2)
        elif match := ev_re.match(line):
            ev_values[(int(match.group(1)), int(match.group(2)))] = float(
                match.group(3)
            )

    if not feat_paths:
        raise ValueError(f"No feat_files() entries found in {path}")
    if not titles or not ev_values:
        raise ValueError(f"No complete group design was found in {path}")

    rows: list[dict[str, Any]] = []
    for fsf_index in sorted(feat_paths):
        input_template = feat_paths[fsf_index]
        subject_match = re.search(r"[/\\]sub-([0-9]+)[/\\]", input_template)
        if not subject_match:
            raise ValueError(
                f"Could not recover a participant ID from: {input_template}"
            )
        row: dict[str, Any] = {
            "fsf_index": fsf_index,
            "sub_id": subject_match.group(1),
            "input_template": input_template,
        }
        for ev_index, title in sorted(titles.items()):
            key = (fsf_index, ev_index)
            if key not in ev_values:
                raise ValueError(f"Missing evg({fsf_index}.{ev_index}) in {path}")
            row[title] = ev_values[key]
        rows.append(row)

    result = pd.DataFrame(rows)
    if result["sub_id"].duplicated().any():
        duplicated = result.loc[result["sub_id"].duplicated(), "sub_id"].tolist()
        raise ValueError(f"Duplicate participant IDs in FSF: {duplicated}")
    return result


def cope_path(input_template: str, cope_number: int) -> Path:
    if "COPENUM" not in input_template:
        raise ValueError(
            f"Expected a COPENUM placeholder in feat_files() path: {input_template}"
        )
    return Path(input_template.replace("COPENUM", str(cope_number), 1))


def extract_roi_mean(image: Path, roi: Path) -> float:
    if not image.is_file():
        raise FileNotFoundError(f"COPE image not found: {image}")

    with tempfile.NamedTemporaryFile(suffix=".txt") as output:
        completed = subprocess.run(
            [
                "fslmeants",
                "-i",
                str(image),
                "-m",
                str(roi),
                "-o",
                output.name,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"fslmeants failed for {image}:\n{completed.stderr.strip()}"
            )
        values = Path(output.name).read_text(encoding="utf-8").split()

    if len(values) != 1:
        raise ValueError(f"Expected one ROI mean from {image}; received {values}")
    value = float(values[0])
    if not (-float("inf") < value < float("inf")):
        raise ValueError(f"Nonfinite ROI mean from {image}: {value}")
    return value


def extract_effects(group_design: Any, roi: Path, pd: Any) -> Any:
    if not roi.is_file():
        raise FileNotFoundError(f"ROI mask not found: {roi}")
    if shutil.which("fslmeants") is None:
        raise RuntimeError("fslmeants is not available on PATH")

    rows: list[dict[str, Any]] = []
    total = len(group_design)
    print(f"Extracting existing COPE 12/13/14 values for {total} participants...")
    for position, row in enumerate(group_design.itertuples(index=False), start=1):
        extracted: dict[str, Any] = {"sub_id": str(row.sub_id)}
        for column, number in COPE_NUMBERS.items():
            image = cope_path(row.input_template, number)
            extracted[column] = extract_roi_mean(image, roi)
        rows.append(extracted)
        if position % 10 == 0 or position == total:
            print(f"  {position}/{total}")
    return pd.DataFrame(rows)


def read_effects(path: Path, pd: Any, np: Any) -> Any:
    if not path.is_file():
        raise FileNotFoundError(f"Effects CSV not found: {path}")
    effects = pd.read_csv(path, dtype={"sub_id": str})
    required = {"sub_id", "nonsocial_cope", "social_cope"}
    missing = required.difference(effects.columns)
    if missing:
        raise ValueError(f"Effects CSV is missing columns: {sorted(missing)}")
    if effects["sub_id"].duplicated().any():
        duplicated = effects.loc[
            effects["sub_id"].duplicated(), "sub_id"
        ].tolist()
        raise ValueError(f"Effects CSV has duplicate participant IDs: {duplicated}")
    if "difference_cope" not in effects:
        effects["difference_cope"] = (
            effects["social_cope"] - effects["nonsocial_cope"]
        )
    numeric_columns = ["nonsocial_cope", "social_cope", "difference_cope"]
    effects[numeric_columns] = effects[numeric_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    if not np.isfinite(effects[numeric_columns].to_numpy()).all():
        raise ValueError("Effects CSV contains missing or nonnumeric COPE values")
    return effects[["sub_id", *numeric_columns]]


def fit_mixed_model(long_data: Any, smf: Any) -> Any:
    covariates = " + ".join(["trust_AIns", *NUISANCE_COVARIATES])
    formula = f"ppi_cope ~ C(partner) * ({covariates})"
    model = smf.mixedlm(
        formula,
        data=long_data,
        groups=long_data["sub_id"],
        re_formula="1",
    )

    attempts = [
        {"method": "lbfgs", "maxiter": 2000},
        {"method": "powell", "maxiter": 5000},
    ]
    last_error: Exception | None = None
    for attempt in attempts:
        try:
            result = model.fit(reml=False, disp=False, **attempt)
            if not result.converged:
                warnings.warn(
                    f"Mixed model did not converge with {attempt['method']}; "
                    "trying the next optimizer.",
                    stacklevel=2,
                )
                continue
            return result
        except Exception as error:  # statsmodels raises several optimizer errors
            last_error = error
    raise RuntimeError(f"Mixed model failed to converge: {last_error}")


def interaction_term_name(result: Any) -> str:
    candidates = [
        name
        for name in result.fe_params.index
        if "C(partner)[T.Social]" in name and "trust_AIns" in name
    ]
    if len(candidates) != 1:
        raise ValueError(
            "Could not identify one partner-by-trust_AIns coefficient; found "
            f"{candidates}"
        )
    return candidates[0]


def linear_combination(
    result: Any,
    terms: list[str],
    label: str,
    np: Any,
) -> dict[str, float | str]:
    names = list(result.fe_params.index)
    absent = set(terms).difference(names)
    if absent:
        raise ValueError(f"Terms absent from model: {sorted(absent)}")
    contrast = np.zeros(len(names), dtype=float)
    for term in terms:
        contrast[names.index(term)] = 1.0

    beta = result.fe_params.to_numpy()
    covariance = result.cov_params().loc[names, names].to_numpy()
    estimate = float(contrast @ beta)
    standard_error = float(np.sqrt(contrast @ covariance @ contrast))
    z_value = estimate / standard_error
    p_value = float(__import__("math").erfc(abs(z_value) / np.sqrt(2.0)))
    return {
        "partner": label,
        "estimate": estimate,
        "std_error": standard_error,
        "conf_low": estimate - 1.959963984540054 * standard_error,
        "conf_high": estimate + 1.959963984540054 * standard_error,
        "z_value": z_value,
        "p_value": p_value,
    }


def fixed_design(result: Any, data: Any, patsy: Any, np: Any) -> Any:
    design_info = result.model.data.design_info
    design = patsy.build_design_matrices([design_info], data)[0]
    design = np.asarray(design, dtype=float)
    if design.shape[1] != len(result.fe_params):
        raise ValueError(
            "Prediction design does not match the fitted fixed-effect coefficients"
        )
    return design


def fixed_prediction(result: Any, data: Any, patsy: Any, np: Any) -> tuple[Any, Any]:
    design = fixed_design(result, data, patsy, np)
    beta = result.fe_params.to_numpy()
    covariance = result.cov_params().loc[
        result.fe_params.index, result.fe_params.index
    ].to_numpy()
    estimate = design @ beta
    standard_error = np.sqrt(np.sum((design @ covariance) * design, axis=1))
    return estimate, standard_error


def model_coefficients(result: Any, pd: Any) -> Any:
    names = list(result.fe_params.index)
    standard_errors = pd.Series(result.bse_fe, index=names)
    z_values = pd.Series(result.tvalues, index=result.params.index).loc[names]
    p_values = pd.Series(result.pvalues, index=result.params.index).loc[names]
    return pd.DataFrame(
        {
            "term": names,
            "estimate": result.fe_params.loc[names].to_numpy(),
            "std_error": standard_errors.to_numpy(),
            "z_value": z_values.to_numpy(),
            "p_value": p_values.to_numpy(),
        }
    )


def make_plots(
    result: Any,
    long_data: Any,
    simple_slopes: Any,
    interaction_name: str,
    out_dir: Path,
    np: Any,
    pd: Any,
    patsy: Any,
    plt: Any,
) -> None:
    colors = {"Nonsocial": "#2C7FB8", "Social": "#D95F0E"}

    # Show model-adjusted lines at the zero point of the original centered
    # nuisance covariates. Partial-residualized points are placed on the same
    # scale so the displayed lines match the fitted interaction model.
    reference = long_data.copy()
    for covariate in NUISANCE_COVARIATES:
        reference[covariate] = 0.0
    observed_fit, _ = fixed_prediction(result, long_data, patsy, np)
    reference_fit, _ = fixed_prediction(result, reference, patsy, np)
    long_data = long_data.copy()
    long_data["adjusted_ppi_cope"] = (
        long_data["ppi_cope"].to_numpy() - observed_fit + reference_fit
    )

    x_values = np.linspace(
        long_data["trust_AIns"].min(), long_data["trust_AIns"].max(), 120
    )
    grid_parts = []
    for partner in ["Nonsocial", "Social"]:
        grid = pd.DataFrame(
            {
                "trust_AIns": x_values,
                "partner": partner,
                **{name: 0.0 for name in NUISANCE_COVARIATES},
            }
        )
        estimate, standard_error = fixed_prediction(result, grid, patsy, np)
        grid["estimate"] = estimate
        grid["conf_low"] = estimate - 1.959963984540054 * standard_error
        grid["conf_high"] = estimate + 1.959963984540054 * standard_error
        grid_parts.append(grid)
    prediction_grid = pd.concat(grid_parts, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7.5, 5.5), constrained_layout=True)
    ax.axhline(0, color="0.75", linewidth=0.8, zorder=0)
    for partner in ["Nonsocial", "Social"]:
        observed = long_data.loc[long_data["partner"] == partner]
        predicted = prediction_grid.loc[prediction_grid["partner"] == partner]
        color = colors[partner]
        ax.scatter(
            observed["trust_AIns"],
            observed["adjusted_ppi_cope"],
            s=25,
            alpha=0.48,
            color=color,
            edgecolor="none",
            label=partner,
        )
        ax.plot(
            predicted["trust_AIns"],
            predicted["estimate"],
            color=color,
            linewidth=2.2,
        )
        ax.fill_between(
            predicted["trust_AIns"],
            predicted["conf_low"],
            predicted["conf_high"],
            color=color,
            alpha=0.16,
            linewidth=0,
        )
    ax.set_xlabel(
        "TG social > nonsocial betrayal-related aIns response\n"
        "(covariate used in the group model)"
    )
    ax.set_ylabel(
        "Fairness-related aIns-vmPFC PPI COPE\n"
        "(adjusted to mean nuisance covariates)"
    )
    ax.set_title(
        "Decomposition of the partner x betrayal interaction",
        loc="left",
        fontweight="bold",
    )
    ax.legend(title="UG partner", frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(out_dir / "ppi_interaction_decomposition.png", dpi=300)
    plt.close(fig)

    interaction_estimate = float(result.fe_params[interaction_name])
    interaction_p = float(result.pvalues[interaction_name])
    x = np.arange(len(simple_slopes))
    y = simple_slopes["estimate"].to_numpy()
    error = np.vstack(
        [
            y - simple_slopes["conf_low"].to_numpy(),
            simple_slopes["conf_high"].to_numpy() - y,
        ]
    )
    fig, ax = plt.subplots(figsize=(6.2, 5.1), constrained_layout=True)
    ax.axhline(0, color="0.45", linewidth=0.9, zorder=0)
    ax.bar(
        x,
        y,
        color=[colors[label] for label in simple_slopes["partner"]],
        width=0.62,
        alpha=0.9,
    )
    ax.errorbar(x, y, yerr=error, fmt="none", ecolor="black", capsize=7, lw=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(simple_slopes["partner"])
    ax.set_ylabel("Association with fairness-related PPI modulation")
    p_text = "< .001" if interaction_p < 0.001 else f"= {interaction_p:.3f}".lstrip("0")
    ax.set_title(
        "Condition-specific betrayal slopes\n"
        f"Partner x betrayal: b = {interaction_estimate:.3f}, p {p_text}",
        loc="left",
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.savefig(out_dir / "ppi_simple_slopes_bar.png", dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    dependencies = load_dependencies()
    np = dependencies["numpy"]
    pd = dependencies["pandas"]
    patsy = dependencies["patsy"]
    smf = dependencies["statsmodels.formula.api"]
    plt = dependencies["matplotlib.pyplot"]

    group_design = parse_group_fsf(args.fsf, pd)
    expected_covariates = {"trust_AIns", *NUISANCE_COVARIATES}
    missing_covariates = expected_covariates.difference(group_design.columns)
    if missing_covariates:
        raise ValueError(
            f"FSF is missing expected covariates: {sorted(missing_covariates)}"
        )

    if args.effects_csv is None:
        if args.roi is None:
            raise ValueError("--roi is required unless --effects-csv is supplied")
        effects = extract_effects(group_design, args.roi, pd)
    else:
        effects = read_effects(args.effects_csv, pd, np)

    unknown_ids = sorted(set(effects["sub_id"]) - set(group_design["sub_id"]))
    if unknown_ids:
        raise ValueError(
            "Effects CSV contains participants absent from the FSF: "
            f"{unknown_ids}"
        )
    missing_ids = sorted(set(group_design["sub_id"]) - set(effects["sub_id"]))
    if missing_ids:
        raise ValueError(
            "No extracted effects were supplied for FSF participants: "
            f"{missing_ids}"
        )
    wide = group_design.merge(
        effects, on="sub_id", how="left", validate="one_to_one", sort=False
    )
    required = [
        "sub_id",
        "trust_AIns",
        *NUISANCE_COVARIATES,
        "nonsocial_cope",
        "social_cope",
        "difference_cope",
    ]
    wide = wide.dropna(subset=required).copy()
    if len(wide) < 10:
        raise ValueError("Fewer than 10 complete participants remain")

    wide["simple_difference"] = wide["social_cope"] - wide["nonsocial_cope"]
    wide["difference_discrepancy"] = (
        wide["difference_cope"] - wide["simple_difference"]
    )

    shared = [
        "sub_id",
        "trust_AIns",
        *NUISANCE_COVARIATES,
    ]
    nonsocial = wide[shared].copy()
    nonsocial["partner"] = "Nonsocial"
    nonsocial["ppi_cope"] = wide["nonsocial_cope"].to_numpy()
    social = wide[shared].copy()
    social["partner"] = "Social"
    social["ppi_cope"] = wide["social_cope"].to_numpy()
    long_data = pd.concat([nonsocial, social], ignore_index=True)
    long_data["partner"] = pd.Categorical(
        long_data["partner"], categories=["Nonsocial", "Social"], ordered=True
    )

    result = fit_mixed_model(long_data, smf)
    interaction_name = interaction_term_name(result)
    simple_slopes = pd.DataFrame(
        [
            linear_combination(result, ["trust_AIns"], "Nonsocial", np),
            linear_combination(
                result,
                ["trust_AIns", interaction_name],
                "Social",
                np,
            ),
        ]
    )

    nuisance = " + ".join(NUISANCE_COVARIATES)
    difference_model = smf.ols(
        f"difference_cope ~ trust_AIns + {nuisance}", data=wide
    ).fit()
    reconstructed_model = smf.ols(
        f"simple_difference ~ trust_AIns + {nuisance}", data=wide
    ).fit()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_columns = [
        "sub_id",
        "nonsocial_cope",
        "social_cope",
        "difference_cope",
        "simple_difference",
        "difference_discrepancy",
        "trust_AIns",
        *NUISANCE_COVARIATES,
    ]
    wide[output_columns].to_csv(
        args.out_dir / "participant_effects_and_covariates.csv", index=False
    )
    model_coefficients(result, pd).to_csv(
        args.out_dir / "model_coefficients.csv", index=False
    )
    simple_slopes.to_csv(args.out_dir / "simple_slopes.csv", index=False)

    correlation = wide[["difference_cope", "simple_difference"]].corr().iloc[0, 1]
    max_discrepancy = wide["difference_discrepancy"].abs().max()
    report = "\n".join(
        [
            "PPI interaction decomposition",
            "================================",
            f"N complete participants: {len(wide)}",
            f"FSF: {args.fsf}",
            f"ROI: {args.roi if args.roi else 'pre-extracted effects CSV'}",
            "",
            "Primary repeated-measures model:",
            str(result.model.formula),
            "",
            result.summary().as_text(),
            "",
            "Condition-specific trust_AIns slopes:",
            simple_slopes.to_string(index=False),
            "",
            "Original extracted COPE-12 difference model:",
            difference_model.summary().as_text(),
            "",
            "COPE-14 minus COPE-13 model:",
            reconstructed_model.summary().as_text(),
            "",
            "Consistency checks:",
            (
                "Correlation between extracted COPE 12 and COPE 14 - COPE 13: "
                f"{correlation:.6f}"
            ),
            f"Maximum absolute discrepancy: {max_discrepancy:.6g}",
            "",
            "Interpretation note:",
            (
                "COPE 13 and COPE 14 quantify condition-specific fairness-related "
                "modulation of aIns-vmPFC coupling. Their signs do not establish "
                "the sign of absolute coupling relative to fixation."
            ),
            (
                "If the ROI was selected from the original interaction map, the "
                "simple-effect tests are descriptive and non-independent. Use an "
                "anatomical or otherwise independent vmPFC ROI for inference."
            ),
        ]
    )
    (args.out_dir / "interaction_model.txt").write_text(report, encoding="utf-8")

    make_plots(
        result,
        long_data,
        simple_slopes,
        interaction_name,
        args.out_dir,
        np,
        pd,
        patsy,
        plt,
    )

    interaction_estimate = result.fe_params[interaction_name]
    interaction_p = result.pvalues[interaction_name]
    print(f"Analysis complete: {args.out_dir.resolve()}")
    print(
        "Partner x trust_AIns interaction: "
        f"b = {interaction_estimate:.6g}, p = {interaction_p:.6g}"
    )
    print(
        "Caution: inference from an ROI selected by this interaction is circular; "
        "treat it as descriptive or use an independent ROI."
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
