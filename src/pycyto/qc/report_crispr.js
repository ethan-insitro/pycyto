"use strict";

// CRISPR report configuration (see WorkflowConfig in report.js). A CRISPR run has no cell
// calls, so the report covers guide capture per probe barcode and guide library coverage.
const WORKFLOW = {
  subtitle: "cyto workflow crispr · guide capture",

  hero: () => [
    [fmt.of(S.guides_detected, S.guides_in_library), "Guides detected"],
    [fmt.big(S.guide_umis), "Guide UMIs"],
    [fmt.int(S.median_umis_per_guide), "Median UMIs per guide"],
    [fmt.num(S.guide_skew_ratio), "Guide skew (90th / 10th percentile)"],
  ],

  summaryTitle: "Guide library",
  summaryRows: () => [
    ["Guides in library", fmt.int(S.guides_in_library)],
    ["Guides with ≥ 1 UMI", `${fmt.int(S.guides_detected)} (${fmt.pct(S.frac_guides_detected)})`],
    ["Median UMIs per guide", fmt.int(S.median_umis_per_guide), "Over every guide in the library, including guides with no UMIs"],
    ["Guide skew", fmt.num(S.guide_skew_ratio), "90th / 10th percentile of UMIs per guide; lower is more even (— if ≥ 10% of guides have no UMIs)"],
    ["Mean mapped reads per probe barcode", fmt.int(S.mean_reads_per_probe)],
  ],

  rankShaded: false,
  rankLegend: html`<span><i class="swatch cell"></i>Guide UMIs per barcode</span>`,
  isActive: () => true,
  probeLabel: (p) => ` · ${fmt.big(p.umis)} UMIs`,
  sortKey: "umis",

  summaryPanels: [
    {
      title: "UMIs per guide",
      draw: (host, probe) => histogram(host, (probe ? PLOTS[probe] : D.pooled).guide_hist, "UMIs per guide", "guides"),
    },
    {
      title: "Most abundant guides",
      draw: (host) => setHTML(host, html`<div class="tablewrap"><table class="data">
        <thead><tr><th>Guide</th><th>UMIs</th><th>% of guide UMIs</th></tr></thead>
        <tbody>${D.pooled.top_guides.map((g) => html`<tr><td class="mono">${g.guide}</td><td>${fmt.int(g.umis)}</td><td>${fmt.pct2(g.frac)}</td></tr>`)}</tbody>
        </table></div><p class="muted">Across all probe barcodes.</p>`),
    },
  ],
  detailPanels: [
    { title: "UMIs per guide", draw: (host, probe) => histogram(host, PLOTS[probe].guide_hist, "UMIs per guide", "guides") },
  ],

  plateMetrics: {
    umis: ["Guide UMIs", fmt.big, true],
    mapped_reads: ["Mapped reads", fmt.big, true],
    guides_detected: ["Guides detected", fmt.int, false],
    n_barcodes: ["Barcodes", fmt.int, true],
    seq_saturation: ["Sequencing saturation", fmt.pct, false],
  },

  detailRows: [
    ["Mapped reads", (p) => html`${fmt.int(p.mapped_reads)}${muted(`(${fmt.pct2(p.frac_of_mapped_reads)} of run)`)}`],
    ["Guide UMIs", (p) => fmt.int(p.umis)],
    ["Guides detected", (p) => fmt.int(p.guides_detected)],
    ["Barcodes observed", (p) => fmt.int(p.n_barcodes)],
    ["Sequencing saturation", (p) => fmt.pct(p.seq_saturation)],
    ["UMIs corrected", (p) => fmt.pct2(p.umi_corrected_frac)],
  ],

  columns: [
    ["probe", "Probe barcode", fmt.text],
    ["mapped_reads", "Mapped reads", fmt.int],
    ["frac_of_mapped_reads", "% of mapped", fmt.pct2],
    ["umis", "Guide UMIs", fmt.int],
    ["guides_detected", "Guides detected", fmt.int],
    ["n_barcodes", "Barcodes", fmt.int],
    ["seq_saturation", "Saturation", fmt.pct],
    ["umi_corrected_frac", "UMIs corrected", fmt.pct2],
  ],

  tableToggle: { label: "Only with guides detected", test: (p) => p.guides_detected > 0 },
};
