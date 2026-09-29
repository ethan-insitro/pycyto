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
 * @typedef {[string, (p: object) => any]} DetailRow  label, value for one probe barcode
 * @typedef {[string, string, (v: any) => string]} Column  key, header, format
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
 * @property {string} sortKey                   probe barcode shown first (largest value)
 * @property {Panel[]} summaryPanels
 * @property {Panel[]} detailPanels
 * @property {object} plateMetrics              key -> [label, format, log scale when wide]
 * @property {DetailRow[]} detailRows
 * @property {Column[]} columns
 * @property {{label: string, test: (p: object) => boolean}} tableToggle
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
  text: (v) => (isNA(v) ? "—" : v),
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
const muted = (text) => html` <span class="muted">${text}</span>`;

/** Calls handler(match) when a click inside `container` lands on (a child of) `selector`. */
function onClick(container, selector, handler) {
  container.addEventListener("click", (e) => {
    const match = e.target.closest(selector);
    if (match && container.contains(match)) handler(match);
  });
}

/** One titled block per panel inside `container` (as cards if `card`); returns the chart hosts. */
function panelHosts(container, panels, card) {
  const heading = raw(card ? "h2" : "h3");
  setHTML(container, panels.map((p) => html`<div class="${card ? "card" : ""}">
    <${heading}>${p.title} <span class="hint"></span></${heading}><div class="panel"></div></div>`));
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

/** All probe barcodes; click a curve to open that probe barcode's details. */
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
        x: "rank", y: "umis", title: (d) => `${d.probe}${WORKFLOW.probeLabel(byProbe[d.probe])}\nclick for details`,
      })),
    ],
  }), (d) => selectProbe(d.probe, true));
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

/** Horizontal bars with a value label. rows: {label, value, color?, tip?, probe?}. */
function hbars(host, rows, { format, max, labelWidth = 170, onPick = null }) {
  const top = max ?? Math.max(...rows.map((r) => r.value), 1e-12);
  const at = { y: "label", insetTop: 5, insetBottom: 5 };
  show(host, () => Plot.plot({
    ...BASE, width: 560, height: rows.length * 26 + 10, marginTop: 5, marginBottom: 5, marginLeft: labelWidth, marginRight: 70,
    x: { domain: [0, top], axis: null },
    y: { domain: rows.map((r) => r.label), label: null, tickSize: 0 },
    color: { type: "identity" },
    marks: [
      Plot.barX(rows, { ...at, x: () => top, fill: color("--empty"), rx: 3 }),
      Plot.barX(rows, { ...at, x: "value", fill: (r) => r.color ?? color("--accent"), rx: 3, tip: true, title: (r) => r.tip ?? format(r.value) }),
      Plot.text(rows, { y: "label", x: () => top, text: (r) => format(r.value), textAnchor: "start", dx: 6, fill: "currentColor" }),
    ],
  }), onPick);
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
  panelHosts("#summary-panels", WORKFLOW.summaryPanels, true).forEach((host, i) => {
    host.parentElement.querySelector(".hint").textContent = probe || "all probe barcodes";
    WORKFLOW.summaryPanels[i].draw(host, probe);
  });
}

// ============================================================================
// Probe barcodes tab
// ============================================================================
let selected = null; // probe barcode shown in the detail panel

/** Moves the `.selected` highlight in the table without redrawing it. */
function highlightSelected() {
  $$("#probe-table .selected").forEach((n) => n.classList.remove("selected"));
  if (selected) $(`#probe-table [data-probe="${CSS.escape(selected)}"]`)?.classList.add("selected");
}

// Flex-V2 probe barcodes are <set>-<row><col>, e.g. A-A01 (underscore also accepted)
const WELL = /^([A-D])[-_]([A-H])(\d{2})$/;
const isPlate = PROBES.length > 0 && PROBES.every((p) => WELL.test(p.probe));

function drawPlates() {
  const host = $("#plates"), key = $("#plate-metric").value, [label, format, logCapable] = WORKFLOW.plateMetrics[key];
  if (!CHARTS_OK) return show(host, null); // the color scales below need d3
  const value = (p) => (isNA(p[key]) || p[key] <= 0 ? null : p[key]); // zero/missing: shown as empty
  const values = PROBES.map(value).filter((v) => v !== null);
  const type = logCapable && values.length && d3.max(values) / Math.max(d3.min(values), 1) > 50 ? "log" : "linear";
  const tip = (p) => `${p.probe}\n${label}: ${format(p[key])}\nMapped reads: ${fmt.big(p.mapped_reads)}`;
  const pick = (d) => d.probe && selectProbe(d.probe);

  if (!isPlate) { // e.g. Flex-V1 (BC001…): one bar per probe barcode
    const scale = (type === "log" ? d3.scaleSequentialLog : d3.scaleSequential)(d3.interpolateViridis).domain(d3.extent(values));
    hbars(host, PROBES.map((p) => ({ label: p.probe, value: value(p) ?? 0, color: value(p) === null ? null : scale(value(p)), tip: tip(p), probe: p.probe })),
      { format, labelWidth: 120, onPick: pick });
    return;
  }

  const wells = PROBES.map((p) => {
    const [, set, row, col] = p.probe.match(WELL);
    return { ...p, set, row, col: +col, value: value(p) };
  });
  const sets = [...new Set(wells.map((w) => w.set))].sort();
  const grid = sets.flatMap((set) => [..."ABCDEFGH"].flatMap((row) => d3.range(1, 13).map((col) => ({ set, row, col }))));
  const cell = { x: "col", y: "row", fx: "set", inset: 1, rx: 3 };
  const flagDot = (flag) => Plot.dot(wells.filter((w) => w.flag === flag), { x: "col", y: "row", fx: "set", r: 2.6, dx: 7, dy: -7, fill: color(`--${flag}`) });
  show(host, () => Plot.plot({
    ...BASE, width: 300 * sets.length, height: 250, marginLeft: 24,
    x: { domain: d3.range(1, 13), label: null, tickSize: 0 },
    y: { domain: [..."ABCDEFGH"], label: null, tickSize: 0 },
    fx: { label: null, tickFormat: (s) => `Set ${s}` },
    color: { type, scheme: "viridis", label, legend: values.length > 0, tickFormat: format },
    marks: [
      Plot.cell(grid, { ...cell, fill: color("--empty") }),
      Plot.cell(wells.filter((w) => w.value !== null), { ...cell, fill: "value" }),
      flagDot("warn"),
      Plot.cell(wells.filter((w) => w.probe === selected), { ...cell, fill: "none", stroke: color("--ink"), strokeWidth: 2 }),
      Plot.tip(wells, Plot.pointer({ x: "col", y: "row", fx: "set", title: tip })),
    ],
  }), pick);
  const flags = PROBES.some((p) => p.flag === "warn") ? html`<span><i class="flag warn"></i>warning</span>` : "";
  host.insertAdjacentHTML("beforeend", render(html`<div class="legend">${flags}<span><i class="swatch empty"></i>no reads / zero</span></div>`));
}

function selectProbe(probe, switchTab = false) {
  const p = byProbe[probe];
  if (!p) return;
  selected = probe;
  if (switchTab) showTab("probes");
  setHTML("#detail-title", html`Probe barcode <span class="mono">${probe}</span>`);
  kvTable("#kv-detail", WORKFLOW.detailRows.map(([label, value]) => [label, value(p)]));
  rankCurve($("#detail-rank"), probe);
  panelHosts("#detail-panels", WORKFLOW.detailPanels, false).forEach((host, i) => WORKFLOW.detailPanels[i].draw(host, probe));
  drawPlates();
  highlightSelected();
}

const sort = { key: "probe", asc: true };
function drawTable() {
  const query = $("#table-filter").value.toLowerCase(), only = $("#table-toggle").checked;
  const rows = PROBES
    .filter((p) => (!only || WORKFLOW.tableToggle.test(p)) && p.probe.toLowerCase().includes(query))
    .sort((a, b) => { // missing values last; probe names sort naturally for Flex-V1 and V2
      const [x, y] = [a[sort.key], b[sort.key]];
      if (isNA(x) || isNA(y)) return isNA(x) - isNA(y);
      return (x < y ? -1 : x > y ? 1 : 0) * (sort.asc ? 1 : -1);
    });
  const sortClass = (key) => (key === sort.key ? `sorted ${sort.asc ? "asc" : ""}` : "");
  setHTML("#probe-table", html`
    <thead><tr>${WORKFLOW.columns.map(([key, label]) => html`<th data-key="${key}" class="${sortClass(key)}">${label}</th>`)}</tr></thead>
    <tbody>${rows.map((p) => html`<tr data-probe="${p.probe}">${WORKFLOW.columns.map(([key, , format], i) =>
      html`<td>${i === 0 && p.flag ? html`<i class="flag ${p.flag}"></i>` : ""}${format(p[key])}</td>`)}</tr>`)}</tbody>`);
  highlightSelected();
}

function setUpProbesTab() {
  setHTML("#plate-metric", Object.entries(WORKFLOW.plateMetrics).map(([k, [label]]) => html`<option value="${k}">${label}</option>`));
  $("#plate-metric").value = Object.keys(WORKFLOW.plateMetrics)[0];
  $("#plate-metric").onchange = drawPlates;
  $("#table-toggle-label").textContent = WORKFLOW.tableToggle.label;
  $("#table-filter").oninput = drawTable;
  $("#table-toggle").onchange = drawTable;
  onClick($("#probe-table"), "th", (th) => {
    sort.asc = sort.key === th.dataset.key ? !sort.asc : th.dataset.key === "probe";
    sort.key = th.dataset.key;
    drawTable();
  });
  onClick($("#probe-table"), "tbody tr", (tr) => {
    selectProbe(tr.dataset.probe);
    $("#detail").scrollIntoView({ behavior: "smooth" });
  });
}

// ============================================================================
// Entry point (called from report.html once report_<workflow>.js has defined WORKFLOW)
// ============================================================================
function main() {
  drawHeader();
  drawSummary();
  setUpProbesTab();
  selected = [...PROBES].sort((a, b) => (b[WORKFLOW.sortKey] ?? 0) - (a[WORKFLOW.sortKey] ?? 0))[0]?.probe ?? null;
  const drawAll = () => {
    drawSummaryPlots();
    drawTable();
    if (selected) selectProbe(selected);
    else drawPlates();
  };
  drawAll();
  // chart colors are resolved from CSS variables, so redraw when light/dark mode flips
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", drawAll);
}
