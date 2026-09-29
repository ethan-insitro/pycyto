"use strict";

// All numbers, flags and labels are computed in Python (pycyto.qc); this file only draws them.
// report.html calls main() once the payload and this script are loaded.
const D = JSON.parse(document.getElementById("data").textContent);
const S = D.summary;
const PROBES = D.probes;
const PLOTS = D.plots;
const BINS = D.log_bins; // log10 bin edges, width 0.05
const byProbe = Object.fromEntries(PROBES.map((p) => [p.probe, p]));
const cellsLabel = (p) => (p.cells ? ` · ${fmt.int(p.cells)} cells` : ""); // suffix in the probe barcode picker
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
// Charts (Observable Plot + d3, loaded from jsDelivr in report.html)
// ============================================================================
/** Resolved value of a CSS color variable (Plot needs concrete colors for scales). */
const color = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const BASE = { style: "background: transparent; color: var(--muted); font-size: 11px; overflow: visible" };
const message = (text) => Object.assign(document.createElement("p"), { className: "muted", textContent: text });

// Plot loads from a CDN; offline, charts show a note and everything else still works
const CHARTS_OK = typeof Plot !== "undefined" && typeof d3 !== "undefined";

/** Shows `plot()` in `host`; clicking a pointed-at datum calls onPick(datum). */
function show(host, plot, onPick = null) {
  if (!CHARTS_OK) return host.replaceChildren(message("Charts need an internet connection (d3 and Observable Plot load from cdn.jsdelivr.net)."));
  const figure = plot();
  if (onPick) figure.addEventListener("click", () => figure.value && onPick(figure.value));
  host.replaceChildren(figure);
}

const rankData = (probe) => PLOTS[probe].curve.map(([rank, umis, frac]) => ({ probe, rank, umis, frac }));
const RANK_AXES = {
  width: 560, height: 340, marginLeft: 54,
  x: { type: "log", label: "Barcodes (rank)", tickFormat: "~s", grid: true },
  y: { type: "log", label: "UMI counts", tickFormat: "~s", grid: true },
};

/** One probe barcode, shaded by the fraction of cells in each segment. */
function rankCurve(host, probe) {
  const data = rankData(probe);
  if (!data.length) return host.replaceChildren(message("No barcodes"));
  show(host, () => Plot.plot({
    ...BASE, ...RANK_AXES,
    color: { type: "linear", domain: [0, 1], range: [color("--bgbc"), color("--cell")] },
    marks: [
      Plot.line(data, { x: "rank", y: "umis", z: null, strokeWidth: 2.5, stroke: "frac" }),
      Plot.tip(data, Plot.pointerX({
        x: "rank", y: "umis",
        title: (d) => `Rank ${fmt.int(d.rank)}\n${fmt.int(d.umis)} UMIs\n${fmt.pct(d.frac, 0)} cells in segment`,
      })),
    ],
  }));
}

/** All probe barcodes (those with cells on top); click a curve to show that probe barcode on its own. */
function rankOverlay(host) {
  const [bg, active] = [PROBES.filter((p) => !p.cells), PROBES.filter((p) => p.cells > 0)]
    .map((ps) => ps.flatMap((p) => rankData(p.probe)));
  if (!bg.length && !active.length) return host.replaceChildren(message("No barcodes"));
  const line = (data, stroke, strokeOpacity) => Plot.line(data, { x: "rank", y: "umis", z: "probe", stroke, strokeOpacity });
  show(host, () => Plot.plot({
    ...BASE, ...RANK_AXES,
    marks: [
      line(bg, color("--bgbc"), 0.35),
      line(active, color("--cell"), 0.6),
      Plot.tip([...bg, ...active], Plot.pointer({
        x: "rank", y: "umis", title: (d) => `${d.probe}${cellsLabel(byProbe[d.probe])}\nclick to show it on its own`,
      })),
    ],
  }), (d) => {
    $("#rank-select").value = d.probe;
    drawSummaryPlots();
  });
}

/** Histogram of cells over the shared log10 bins. */
function histogram(host, counts, title) {
  const bins = (counts ?? []).flatMap((n, i) => (n ? [{ lo: 10 ** BINS[i], hi: 10 ** (BINS[i] + 0.05), n }] : []));
  if (!bins.length) return host.replaceChildren(message("No cells"));
  show(host, () => Plot.plot({
    ...BASE, width: 560, height: 230, marginLeft: 54,
    x: { type: "log", label: title, tickFormat: "~s" },
    y: { label: "Cells", tickFormat: "~s", grid: true },
    marks: [
      Plot.rectY(bins, {
        x1: "lo", x2: "hi", y: "n", fill: color("--cell"), insetLeft: 0.5, tip: true,
        title: (d) => `${fmt.int(d.lo)}–${fmt.int(d.hi)}\n${fmt.int(d.n)} cells`,
      }),
      Plot.ruleY([0]),
    ],
  }));
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
  setHTML("#rank-select", [
    html`<option value="">All probe barcodes (overlay)</option>`,
    ...PROBES.map((p) => html`<option value="${p.probe}">${p.probe}${cellsLabel(p)}</option>`),
  ]);
  $("#rank-select").onchange = drawSummaryPlots;
}

/** Rank plot and histograms for the probe barcode picked above the rank plot ("" = all). */
function drawSummaryPlots() {
  const probe = $("#rank-select").value;
  if (probe) rankCurve($("#rank-plot"), probe);
  else rankOverlay($("#rank-plot"));
  $$("#summary-panels .hint").forEach((node) => (node.textContent = probe || "all probe barcodes"));
  const hists = probe ? PLOTS[probe] : D.pooled;
  histogram($("#umi-hist"), hists.umi_hist, "UMIs per cell");
  histogram($("#gene-hist"), hists.gene_hist, "Genes per cell");
}

// ============================================================================
// Entry point (called from report.html)
// ============================================================================
function main() {
  drawHeader();
  drawSummary();
  const drawAll = () => drawSummaryPlots();
  drawAll();
  // chart colors are resolved from CSS variables, so redraw when light/dark mode flips
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
}
