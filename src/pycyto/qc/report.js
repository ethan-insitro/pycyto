"use strict";

// All numbers, flags and labels are computed in Python (pycyto.qc); this file only draws them.
// report.html calls main() once the payload and this script are loaded.
const D = JSON.parse(document.getElementById("data").textContent);
const S = D.summary;
const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => [...document.querySelectorAll(sel)];

// ============================================================================
// Formatting and HTML
// ============================================================================
const isNA = (v) => v === null || v === undefined || (typeof v === "number" && !Number.isFinite(v));
const fmt = {
  int: (v) => (isNA(v) ? "—" : Math.round(v).toLocaleString()),
  num: (v, digits = 1) => (isNA(v) ? "—" : v.toFixed(digits)),
  pct: (v, digits = 1) => (isNA(v) ? "—" : `${(v * 100).toFixed(digits)}%`),
  pct2: (v) => fmt.pct(v, 2),
  big: (v) => {
    if (isNA(v)) return "—";
    if (Math.abs(v) >= 1e9) return `${(v / 1e9).toFixed(2)}B`;
    if (Math.abs(v) >= 1e6) return `${(v / 1e6).toFixed(1)}M`;
    if (Math.abs(v) >= 1e4) return `${(v / 1e3).toFixed(1)}k`;
    return fmt.int(v);
  },
  of: (n, total) => (isNA(n) ? "—" : `${fmt.int(n)} / ${fmt.int(total)}`),
};

/** Tagged template that HTML-escapes every interpolated value unless wrapped in raw(). */
class Raw {
  constructor(s) { this.s = s; }
  toString() { return this.s; }
}
const raw = (s) => new Raw(s);
const escape = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
const render = (v) => (v instanceof Raw ? v.s : Array.isArray(v) ? v.map(render).join("") : v === null || v === undefined ? "" : escape(v));
const html = (strings, ...values) => raw(strings.reduce((out, s, i) => out + render(values[i - 1]) + s));
const setHTML = (target, content) => ((typeof target === "string" ? $(target) : target).innerHTML = render(content));

const kvTable = (sel, rows) => setHTML(sel, rows.map(([k, v, title]) => html`<tr title="${title ?? ""}"><td>${k}</td><td>${v}</td></tr>`));

/** Calls handler(match) when a click inside `container` lands on (a child of) `selector`. */
function onClick(container, selector, handler) {
  container.addEventListener("click", (e) => {
    const match = e.target.closest(selector);
    if (match && container.contains(match)) handler(match);
  });
}

// ============================================================================
// Header, tabs, summary tab
// ============================================================================
function showTab(name) {
  $$("nav button").forEach((b) => b.classList.toggle("active", b.dataset.tab === name));
  $$("section.tab").forEach((s) => s.classList.toggle("active", s.id === `tab-${name}`));
  history.replaceState(null, "", `#${name}`);
}

function drawHeader() {
  $("#title").textContent = D.title;
  $("#generated").textContent = `cyto workflow ${D.workflow} · generated ${D.generated}`;
  $("#footer").textContent = `pycyto ${D.version} · ${S.cyto_outdir}`;
  onClick($("nav"), "button", (b) => showTab(b.dataset.tab));
  if ($(`nav button[data-tab="${location.hash.slice(1)}"]`)) showTab(location.hash.slice(1));
}

function drawSummary() {
  setHTML("#hero", [
    [fmt.int(S.estimated_cells), "Estimated number of cells"],
    [fmt.int(S.mean_reads_per_cell), "Mean reads per cell"],
    [fmt.int(S.median_genes_per_cell), "Median genes per cell"],
    [fmt.int(S.median_umis_per_cell), "Median UMI counts per cell"],
  ].map(([v, label]) => html`<div class="card hero"><div class="value">${v}</div><div class="label">${label}</div></div>`));
  const inLibrary = S.probe_barcodes_in_library ? ` / ${fmt.int(S.probe_barcodes_in_library)} in library` : "";
  kvTable("#kv-sequencing", [
    ["Number of reads", fmt.int(S.total_reads)],
    ["Reads mapped to probes", fmt.pct(S.mapped_reads_frac), "Reads with a valid cell barcode, probe barcode, UMI and feature match"],
    ["Sequencing saturation", fmt.pct(S.seq_saturation), "1 − UMIs / mapped reads, over all barcodes"],
    ["UMIs corrected", fmt.pct2(S.umi_corrected_frac)],
    ["Probe barcodes with reads", `${fmt.int(S.probe_barcodes_with_reads)}${inLibrary}`],
  ]);
  kvTable("#kv-cells", [
    ["Probe barcodes with cells", fmt.of(S.probe_barcodes_with_cells, S.probe_barcodes_with_reads)],
    ["Median cells per probe barcode", fmt.int(S.cells_median_per_probe), "Among probe barcodes with cells"],
    ["Mean mapped reads per cell", fmt.int(S.mean_mapped_reads_per_cell)],
    ["Fraction reads in cells", fmt.pct(S.frac_reads_in_cells), "Mapped reads in cell barcodes / all mapped reads"],
    ["Mapped reads in probe barcodes without cells", fmt.pct(S.background_probe_read_frac)],
    ["Total genes detected", fmt.of(S.total_genes_detected, S.genes_in_reference)],
  ]);
}

// ============================================================================
// Entry point (called from report.html)
// ============================================================================
function main() {
  drawHeader();
  drawSummary();
}
