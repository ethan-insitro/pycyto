"""Cell Ranger-style alerts raised from the QC metrics."""

import statistics

# Alert thresholds. Tuples are (warn, error) where applicable.
THRESH = {
    "mapped_frac": (0.70, 0.50),  # below -> warn / error
    "frac_reads_in_cells": (0.70, 0.50),  # below -> warn / error
    "failed_umi_qual_of_total": 0.10,  # above -> warn
    "background_probe_read_frac": 0.05,  # mapped reads in probe barcodes w/o cells; above -> warn
    "median_umis_per_cell": 500,  # below -> warn (per probe barcode)
    "cells_cv": 0.5,  # CV of cells across probe barcodes with cells; above -> warn
}


def _examples(items: list[str], limit: int = 10) -> str:
    return ", ".join(items[:limit]) + (" …" if len(items) > limit else "")


def build_alerts(summary: dict, recs: list[dict]) -> list[dict]:
    alerts: list[dict] = []

    def add(level: str, title: str, detail: str) -> None:
        alerts.append({"level": level, "title": title, "detail": detail})

    mf = summary.get("mapped_reads_frac")
    if mf is not None:
        warn, err = THRESH["mapped_frac"]
        if mf < err:
            add(
                "error",
                "Low fraction of reads mapped",
                f"{mf:.1%} of reads mapped (expected ≥ {warn:.0%}). Check the probe set, "
                "whitelist, and geometry/preset.",
            )
        elif mf < warn:
            add(
                "warn",
                "Low fraction of reads mapped",
                f"{mf:.1%} of reads mapped (expected ≥ {warn:.0%}). The biggest unmapped "
                f"category is “{summary.get('top_unmapped_reason')}”.",
            )

    fu = summary.get("failed_umi_qual_of_total")
    if fu is not None and fu > THRESH["failed_umi_qual_of_total"]:
        add(
            "warn",
            "Many reads failed UMI quality",
            f"{fu:.1%} of all reads failed the UMI quality filter. This can point to low "
            "base quality in R1.",
        )

    fr = summary.get("frac_reads_in_cells")
    if fr is not None:
        warn, err = THRESH["frac_reads_in_cells"]
        if fr < err:
            add(
                "error",
                "Low fraction of reads in cells",
                f"{fr:.1%} of mapped reads are in cell barcodes (expected ≥ {warn:.0%}). "
                "Suggests high ambient RNA or failed cell calling.",
            )
        elif fr < warn:
            add(
                "warn",
                "Low fraction of reads in cells",
                f"{fr:.1%} of mapped reads are in cell barcodes (expected ≥ {warn:.0%}).",
            )

    bg = summary.get("background_probe_read_frac")
    if bg is not None and bg > THRESH["background_probe_read_frac"]:
        add(
            "warn",
            "Reads in probe barcodes without cells",
            f"{bg:.1%} of mapped reads went to {summary['n_probes_without_cells']} probe "
            "barcodes with no cells. Check for unexpected probe barcodes or barcode hopping.",
        )

    called = [r for r in recs if r["cells"] > 0]

    bad_bound = [
        r
        for r in called
        if r.get("cyto_retain_umis") is not None
        and r.get("cyto_reject_umis") is not None
        and r["cyto_retain_umis"] < r["cyto_reject_umis"]
    ]
    if bad_bound:
        lst = ", ".join(
            f"{r['probe']} ({r['cells']:,} cells, retain ≥{r['cyto_retain_umis']} UMIs)"
            for r in bad_bound
        )
        add(
            "warn",
            "cyto cell calling looks unreliable",
            f"For {len(bad_bound)} probe barcode(s), cyto auto-accepted barcodes below its "
            f"own rejection boundary, so almost every barcode was called a cell: {lst}. "
            "These are likely near-empty.",
        )

    low = [
        r
        for r in called
        if r["median_umis_per_cell"] is not None
        and r["median_umis_per_cell"] < THRESH["median_umis_per_cell"]
    ]
    if low:
        add(
            "warn",
            "Low median UMIs per cell",
            f"{len(low)} probe barcode(s) with cells have a median below "
            f"{THRESH['median_umis_per_cell']} UMIs per cell: "
            + _examples([f"{r['probe']} ({r['median_umis_per_cell']:.0f})" for r in low])
            + ".",
        )

    low_probes = {r["probe"] for r in low}
    err_frac = THRESH["frac_reads_in_cells"][1]
    lowfr = [
        r
        for r in called
        if r["probe"] not in low_probes
        and r["frac_reads_in_cells"] is not None
        and r["frac_reads_in_cells"] < err_frac
    ]
    if lowfr:
        add(
            "warn",
            "Probe barcodes with low reads in cells",
            f"{len(lowfr)} probe barcode(s) have <{err_frac:.0%} of reads in cells: "
            + _examples([f"{r['probe']} ({r['frac_reads_in_cells']:.0%})" for r in lowfr])
            + ".",
        )

    if len(called) >= 3:
        counts = [r["cells"] for r in called]
        cv = statistics.stdev(counts) / statistics.mean(counts)
        if cv > THRESH["cells_cv"]:
            add(
                "warn",
                "Uneven cell counts across probe barcodes",
                f"Coefficient of variation of cells per probe barcode is {cv:.2f}.",
            )

    if not alerts:
        add("ok", "No issues detected", "All checked metrics are within expected ranges.")
    return alerts
