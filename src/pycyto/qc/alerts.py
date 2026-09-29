"""Cell Ranger-style alerts for every workflow."""

import polars as pl

# Alert thresholds. Tuples are (warn, error).
THRESH = {
    # all workflows
    "mapped_frac": (0.70, 0.50),  # below -> warn / error
    "failed_umi_qual_of_total": 0.10,  # above -> warn
    # gex
    "frac_reads_in_cells": (0.70, 0.50),  # below -> warn / error
    "background_probe_read_frac": 0.05,  # reads in probe barcodes w/o cells; above -> warn
    "median_umis_per_cell": 500,  # below -> warn (per probe barcode)
    "cells_cv": 0.5,  # CV of cells across probe barcodes with cells; above -> warn
}

def _below(value: float | None, warn: float, error: float) -> str | None:
    if value is None:
        return None
    return "error" if value < error else "warn" if value < warn else None


def _above(value: float | None, warn: float) -> str | None:
    return "warn" if value is not None and value > warn else None


def _examples(rows: pl.DataFrame, fmt, limit: int = 10) -> str:
    items = [fmt(r) for r in rows.head(limit).iter_rows(named=True)]
    return ", ".join(items) + (" …" if rows.height > limit else "")


def build_alerts(workflow: str, summary: dict, probes: pl.DataFrame) -> list[dict]:
    alerts: list[dict] = []

    def add(level: str | None, title: str, detail: str) -> None:
        if level:
            alerts.append({"level": level, "title": title, "detail": detail})

    # --- all workflows --------------------------------------------------------------
    mf = summary["mapped_reads_frac"]
    add(
        _below(mf, *THRESH["mapped_frac"]),
        "Low fraction of reads mapped",
        f"{mf or 0:.1%} of reads mapped (expected ≥ {THRESH['mapped_frac'][0]:.0%}). "
        f"The biggest unmapped category is “{summary['top_unmapped_reason']}”.",
    )
    fu = summary["failed_umi_qual_of_total"]
    add(
        _above(fu, THRESH["failed_umi_qual_of_total"]),
        "Many reads failed UMI quality",
        f"{fu or 0:.1%} of all reads failed the UMI quality filter. "
        "This can point to low base quality in R1.",
    )
    if workflow == "gex":
        _gex_alerts(add, summary, probes)

    return alerts or [
        {"level": "ok", "title": "No issues detected", "detail": "All checked metrics are within expected ranges."}
    ]


def _gex_alerts(add, summary: dict, probes: pl.DataFrame) -> None:
    fr = summary["frac_reads_in_cells"]
    add(
        _below(fr, *THRESH["frac_reads_in_cells"]),
        "Low fraction of reads in cells",
        f"{fr or 0:.1%} of mapped reads are in cell barcodes "
        f"(expected ≥ {THRESH['frac_reads_in_cells'][0]:.0%}).",
    )
    bg = summary["background_probe_read_frac"]
    add(
        _above(bg, THRESH["background_probe_read_frac"]),
        "Reads in probe barcodes without cells",
        f"{bg or 0:.1%} of mapped reads went to {summary['n_probes_without_cells']} probe "
        "barcodes with no cells. Check for unexpected probe barcodes or barcode hopping.",
    )

    called = probes.filter(pl.col("cells") > 0)
    low = called.filter(pl.col("median_umis_per_cell") < THRESH["median_umis_per_cell"])
    add(
        "warn" if low.height else None,
        "Low median UMIs per cell",
        f"{low.height} probe barcode(s) with cells have a median below "
        f"{THRESH['median_umis_per_cell']} UMIs per cell: "
        + _examples(low, lambda r: f"{r['probe']} ({r['median_umis_per_cell']:.0f})")
        + ".",
    )
    min_frac = THRESH["frac_reads_in_cells"][1]
    low_frac = called.filter(
        (pl.col("frac_reads_in_cells") < min_frac) & ~pl.col("probe").is_in(low["probe"])
    )
    add(
        "warn" if low_frac.height else None,
        "Probe barcodes with low reads in cells",
        f"{low_frac.height} probe barcode(s) have <{min_frac:.0%} of reads in cells: "
        + _examples(low_frac, lambda r: f"{r['probe']} ({r['frac_reads_in_cells']:.0%})")
        + ".",
    )
    if called.height >= 3:
        cv = called["cells"].std() / called["cells"].mean()
        add(
            _above(cv, THRESH["cells_cv"]),
            "Uneven cell counts across probe barcodes",
            f"Coefficient of variation of cells per probe barcode is {cv:.2f}.",
        )
