# 🌟 Clean Repository Structure (2025-08-03)

## 📁 Essential Directory Structure

```
bomin/
├── 📊 0.Data/                     # All data files
│   ├── GWAS/                      # GWAS source data
│   ├── Enhancer/                  # Cell-type enhancer BED files
│   ├── Reference/                 # LDSC reference data
│   └── Results/                   # Processed data
│       ├── annotations/           # LDSC annotations (176 files)
│       └── combined_ld_scores/    # LD scores (176 files)
│
├── 💻 1.Scripts/                  # Essential scripts only
│   ├── LDSC/
│   │   ├── ldsc_analysis_system.py  # Main LDSC pipeline
│   │   └── ldsc/                    # Original LDSC software
│   ├── Utils/                     # Utility scripts
│   └── Visualization/             # Visualization scripts
│
├── 🔬 2.Analysis/                 # Analysis results
│   └── LDSC/
│       ├── final_analysis/        # Final results (8/8 datasets)
│       ├── visualizations/        # Generated plots
│       └── visualize_8_datasets.py # Main visualization script
│
├── 📈 4.figures/                  # All figures
│   ├── enrichment_all_8_datasets.png
│   ├── pvalue_all_8_datasets.png
│   ├── combined_all_8_datasets.png
│   ├── ranking_all_8_datasets.png
│   └── summary_statistics.json
│
├── 📚 README.md                   # Main documentation
└── 🗄️ backup_030825/              # All archived files

```

## 🎯 Key Files Locations

### Data Files
- **GWAS Data**: `0.Data/GWAS/GCST009325.h.tsv.gz`
- **Enhancer BED**: `0.Data/Enhancer/*.bed` (8 files)
- **LD Scores**: `0.Data/Results/combined_ld_scores/*.l2.ldscore.gz`

### Essential Scripts
- **LDSC Pipeline**: `1.Scripts/LDSC/ldsc_analysis_system.py`
- **Visualization**: `2.Analysis/LDSC/visualize_8_datasets.py`

### Results
- **Final Analysis**: `2.Analysis/LDSC/final_analysis/`
- **Visualizations**: `4.figures/`

## 🚀 Quick Start

1. **View Results**: Check `4.figures/` for all visualizations
2. **Read Analysis**: See `2.Analysis/LDSC/final_analysis/complete_8_datasets_report.md`
3. **Run Visualization**: 
   ```bash
   cd 2.Analysis/LDSC
   python visualize_8_datasets.py
   ```

## ⚠️ Important Notes

- All test files, experimental code, and intermediate results have been moved to `backup_030825/`
- The current structure contains only essential files for the LDSC analysis
- All data paths are documented in the main `README.md`

## 📊 Analysis Summary

- **8/8 datasets analyzed** (4 cell types × 2 processing methods)
- **All datasets significant** (Bonferroni corrected p < 0.00625)
- **Top enrichment**: Oligodendrocytes (Olig_cleaned: 4.041)
- **Mean enrichment**: 3.267 ± 0.611

---
*Repository cleaned and organized on 2025-08-03*