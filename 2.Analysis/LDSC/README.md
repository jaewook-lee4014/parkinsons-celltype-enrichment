# LDSC Analysis Directory

This directory contains the essential LDSC analysis files for the Parkinson's Disease GWAS enrichment study.

## Essential Files

### 📊 Final Analysis Results
- `final_analysis/` - Complete analysis results for all 8/8 datasets
  - `final_enrichment_results.csv` - Enrichment values and statistics
  - `final_enrichment_results.json` - JSON format results
  - `complete_8_datasets_report.md` - Comprehensive report

### 🎨 Visualization
- `visualize_8_datasets.py` - Main visualization script
- `visualizations/` - Generated plots
  - `enrichment_all_8_datasets.png` - Enrichment bar chart
  - `pvalue_all_8_datasets.png` - Statistical significance plot
  - `combined_all_8_datasets.png` - Combined analysis view
  - `ranking_all_8_datasets.png` - Cell type ranking plot

## Data Paths
All data files referenced in README.md at the repository root:
- GWAS data: `/cephfs/.../0.Data/GWAS/`
- Enhancer BED files: `/cephfs/.../0.Data/Enhancer/`
- LD scores: `/scratch/.../combined_ld_scores/`

## Note
Additional test files and experimental code have been moved to `backup_030825/`