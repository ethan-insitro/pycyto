"use strict";

// GEX report configuration (see WorkflowConfig in report.js). Cells are cyto's filtered h5ad.
const WORKFLOW = {
  subtitle: "cyto workflow gex · cells from cyto filtered h5ad",

  hero: () => [
    [fmt.int(S.estimated_cells), "Estimated number of cells"],
    [fmt.int(S.mean_reads_per_cell), "Mean reads per cell"],
    [fmt.int(S.median_genes_per_cell), "Median genes per cell"],
    [fmt.int(S.median_umis_per_cell), "Median UMI counts per cell"],
  ],

  summaryTitle: "Cells",
  summaryRows: () => [
    ["Probe barcodes with cells", fmt.of(S.probe_barcodes_with_cells, S.probe_barcodes_with_reads)],
    ["Median cells per probe barcode", fmt.int(S.cells_median_per_probe), "Among probe barcodes with cells"],
    ["Mean mapped reads per cell", fmt.int(S.mean_mapped_reads_per_cell)],
    ["Fraction reads in cells", fmt.pct(S.frac_reads_in_cells), "Mapped reads in cell barcodes / all mapped reads"],
    ["Mapped reads in probe barcodes without cells", fmt.pct(S.background_probe_read_frac)],
    ["Total genes detected", fmt.of(S.total_genes_detected, S.genes_in_reference)],
  ],

  rankShaded: true,
  rankLegend: html`<span><i class="swatch cell"></i>Cells</span><span><i class="swatch bg"></i>Background</span>`,
  isActive: (p) => p.cells > 0,
  probeLabel: (p) => (p.cells ? ` · ${fmt.int(p.cells)} cells` : ""),

  summaryPanels: [
    { title: "UMIs per cell", draw: (host, probe) => histogram(host, (probe ? PLOTS[probe] : D.pooled).umi_hist, "UMIs per cell") },
    { title: "Genes per cell", draw: (host, probe) => histogram(host, (probe ? PLOTS[probe] : D.pooled).gene_hist, "Genes per cell") },
  ],

  plateMetrics: {
    cells: ["Cells", fmt.int, true],
    mapped_reads: ["Mapped reads", fmt.big, true],
    median_umis_per_cell: ["Median UMIs/cell", fmt.int, true],
    median_genes_per_cell: ["Median genes/cell", fmt.int, false],
    frac_reads_in_cells: ["Fraction reads in cells", fmt.pct, false],
    seq_saturation: ["Sequencing saturation", fmt.pct, false],
  },

  columns: [
    ["probe", "Probe barcode", fmt.text],
    ["cells", "Cells", fmt.int],
    ["mapped_reads", "Mapped reads", fmt.int],
    ["frac_of_mapped_reads", "% of mapped", fmt.pct2],
    ["mean_reads_per_cell", "Mean reads/cell", fmt.int],
    ["median_umis_per_cell", "Median UMIs/cell", fmt.int],
    ["median_genes_per_cell", "Median genes/cell", fmt.int],
    ["frac_reads_in_cells", "Reads in cells", fmt.pct],
    ["seq_saturation", "Saturation", fmt.pct],
    ["total_genes_detected", "Genes detected", fmt.int],
    ["n_barcodes", "Barcodes", fmt.int],
  ],

  tableToggle: { label: "Only with cells", test: (p) => p.cells > 0 },
};
