"use strict";

// Shared report engine for every cyto workflow. All numbers, flags and labels are computed
// in Python (pycyto.qc); this file only draws them. What differs between workflows (headline
// metrics, table columns, extra panels) lives in report_<workflow>.js, which defines
// `WORKFLOW` (see the WorkflowConfig typedef below); report.html then calls main().
const D = JSON.parse(document.getElementById("data").textContent);
const S = D.summary;
const PROBES = D.probes;
const PLOTS = D.plots;
const BINS = D.log_bins; // log10 bin edges, width 0.05
const byProbe = Object.fromEntries(PROBES.map((p) => [p.probe, p]));
const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => [...document.querySelectorAll(sel)];

/**
 * @typedef {{title: string, draw: (host: Element, probe: string) => void}} Panel
 *   probe is "" when the summary shows all probe barcodes
 * @typedef {object} WorkflowConfig
 * @property {string} subtitle                  shown in the header
 * @property {() => Array} hero                 [value, label] headline numbers
 * @property {string} summaryTitle              title of the second summary table
 * @property {() => Array} summaryRows          [label, value, tooltip?] rows for it
 * @property {boolean} rankShaded               shade rank curves by fraction of cells
 * @property {object} rankLegend                html`` legend under rank plots
 * @property {(p: object) => boolean} isActive  probe barcodes drawn on top in the overlay
 * @property {(p: object) => string} probeLabel suffix in the probe barcode picker
 * @property {Panel[]} summaryPanels
 */

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

/** One titled card per panel inside `container`; returns the chart hosts. */
function panelHosts(container, panels) {
  setHTML(container, panels.map((p) => html`<div class="card">
    <h2>${p.title} <span class="hint"></span></h2><div class="panel"></div></div>`));
  return [...$(container).querySelectorAll(".panel")];
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

/** One probe barcode; if WORKFLOW.rankShaded, shaded by the fraction of cells in each segment. */
function rankCurve(host, probe) {
  const data = rankData(probe);
  if (!data.length) return host.replaceChildren(message("No barcodes"));
  const shaded = WORKFLOW.rankShaded;
  show(host, () => Plot.plot({
    ...BASE, ...RANK_AXES,
    color: { type: "linear", domain: [0, 1], range: [color("--bgbc"), color("--cell")] },
    marks: [
      Plot.line(data, { x: "rank", y: "umis", z: null, strokeWidth: 2.5, stroke: shaded ? "frac" : color("--cell") }),
      Plot.tip(data, Plot.pointerX({
        x: "rank", y: "umis",
        title: (d) => `Rank ${fmt.int(d.rank)}\n${fmt.int(d.umis)} UMIs` + (shaded ? `\n${fmt.pct(d.frac, 0)} cells in segment` : ""),
      })),
    ],
  }));
}

/** All probe barcodes; click a curve to show that probe barcode on its own. */
function rankOverlay(host) {
  const [bg, active] = [PROBES.filter((p) => !WORKFLOW.isActive(p)), PROBES.filter(WORKFLOW.isActive)]
    .map((ps) => ps.flatMap((p) => rankData(p.probe)));
  if (!bg.length && !active.length) return host.replaceChildren(message("No barcodes"));
  const line = (data, stroke, strokeOpacity) => Plot.line(data, { x: "rank", y: "umis", z: "probe", stroke, strokeOpacity });
  show(host, () => Plot.plot({
    ...BASE, ...RANK_AXES,
    marks: [
      line(bg, color("--bgbc"), 0.35),
      line(active, color("--cell"), 0.6),
      Plot.tip([...bg, ...active], Plot.pointer({
        x: "rank", y: "umis", title: (d) => `${d.probe}${WORKFLOW.probeLabel(byProbe[d.probe])}\nclick to show it on its own`,
      })),
    ],
  }), (d) => {
    $("#rank-select").value = d.probe;
    drawSummaryPlots();
  });
}

/** Histogram over the shared log10 bins; `unit` names what is counted (cells, guides). */
function histogram(host, counts, title, unit = "cells") {
  const bins = (counts ?? []).flatMap((n, i) => (n ? [{ lo: 10 ** BINS[i], hi: 10 ** (BINS[i] + 0.05), n }] : []));
  if (!bins.length) return host.replaceChildren(message(`No ${unit}`));
  show(host, () => Plot.plot({
    ...BASE, width: 560, height: 230, marginLeft: 54,
    x: { type: "log", label: title, tickFormat: "~s" },
    y: { label: unit[0].toUpperCase() + unit.slice(1), tickFormat: "~s", grid: true },
    marks: [
      Plot.rectY(bins, {
        x1: "lo", x2: "hi", y: "n", fill: color("--cell"), insetLeft: 0.5, tip: true,
        title: (d) => `${fmt.int(d.lo)}–${fmt.int(d.hi)}\n${fmt.int(d.n)} ${unit}`,
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
  $("#generated").textContent = `${WORKFLOW.subtitle} · generated ${D.generated}`;
  $("#footer").textContent = `pycyto ${D.version} · ${S.cyto_outdir}`;
  onClick($("nav"), "button", (b) => showTab(b.dataset.tab));
  if ($(`nav button[data-tab="${location.hash.slice(1)}"]`)) showTab(location.hash.slice(1));
  $$(".rank-legend").forEach((node) => setHTML(node, WORKFLOW.rankLegend));
}

function drawSummary() {
  setHTML("#alerts", D.alerts.map((a) => html`<div class="alert ${a.level}"><b>${a.title}</b>${a.detail}</div>`));
  setHTML("#hero", WORKFLOW.hero().map(([v, label]) => html`<div class="card hero"><div class="value">${v}</div><div class="label">${label}</div></div>`));
  const inLibrary = S.probe_barcodes_in_library ? ` / ${fmt.int(S.probe_barcodes_in_library)} in library` : "";
  kvTable("#kv-sequencing", [
    ["Number of reads", fmt.int(S.total_reads)],
    ["Reads mapped to probes", fmt.pct(S.mapped_reads_frac), "Reads with a valid cell barcode, probe barcode, UMI and feature match"],
    ["Sequencing saturation", fmt.pct(S.seq_saturation), "1 − UMIs / mapped reads, over all barcodes"],
    ["UMIs corrected", fmt.pct2(S.umi_corrected_frac)],
    ["Probe barcodes with reads", `${fmt.int(S.probe_barcodes_with_reads)}${inLibrary}`],
  ]);
  $("#kv-workflow-title").textContent = WORKFLOW.summaryTitle;
  kvTable("#kv-workflow", WORKFLOW.summaryRows());
  setHTML("#rank-select", [
    html`<option value="">All probe barcodes (overlay)</option>`,
    ...PROBES.map((p) => html`<option value="${p.probe}">${p.probe}${WORKFLOW.probeLabel(p)}</option>`),
  ]);
  $("#rank-select").onchange = drawSummaryPlots;
}

/** Rank plot and summary panels for the probe barcode picked above the rank plot ("" = all). */
function drawSummaryPlots() {
  const probe = $("#rank-select").value;
  if (probe) rankCurve($("#rank-plot"), probe);
  else rankOverlay($("#rank-plot"));
  panelHosts("#summary-panels", WORKFLOW.summaryPanels).forEach((host, i) => {
    host.parentElement.querySelector(".hint").textContent = probe || "all probe barcodes";
    WORKFLOW.summaryPanels[i].draw(host, probe);
  });
}

// ============================================================================
// Entry point (called from report.html once report_<workflow>.js has defined WORKFLOW)
// ============================================================================
function main() {
  drawHeader();
  drawSummary();
  const drawAll = () => drawSummaryPlots();
  drawAll();
  // chart colors are resolved from CSS variables, so redraw when light/dark mode flips
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
}
