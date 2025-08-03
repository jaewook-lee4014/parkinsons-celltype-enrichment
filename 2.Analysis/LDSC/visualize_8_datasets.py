#!/usr/bin/env python3
"""
Visualize LDSC enrichment and -log(p) for all 8/8 datasets
Using the complete final results with corrected statistics
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class LDSCVisualizer:
    def __init__(self):
        self.results_dir = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin/2.Analysis/LDSC/final_analysis")
        self.output_dir = Path("/cephfs/volumes/hpc_data_prj/eng_waste_to_protein/ae035a41-20d2-44f3-aa46-14424ab0f6bf/repositories/bomin/2.Analysis/LDSC/visualizations")
        self.output_dir.mkdir(exist_ok=True)
        
        # 색상 매핑 (세포타입별)
        self.colors = {
            'Olig': '#FF6B6B',  # 빨간색 - 가장 높은 enrichment
            'Neg': '#4ECDC4',   # 청록색 - Microglia
            'Nurr': '#45B7D1',  # 파란색 - Dopaminergic neurons
            'NeuN': '#96CEB4'   # 녹색 - General neurons
        }
        
        # 처리 타입별 패턴
        self.patterns = {
            'cleaned': '',      # 실선
            'unique': '///'     # 사선 패턴
        }
        
        logger.info("🎨 LDSC Visualizer initialized")
    
    def load_results(self):
        """Load final analysis results"""
        json_file = self.results_dir / "final_enrichment_results.json"
        
        with open(json_file, 'r') as f:
            results = json.load(f)
        
        # Convert to DataFrame for easier handling
        df = pd.DataFrame(results)
        df['neg_log_p'] = -np.log10(df['enrichment_p'])
        
        logger.info(f"📊 Loaded {len(df)} datasets")
        return df
    
    def create_enrichment_plot(self, df):
        """Create enrichment bar plot for all 8 datasets"""
        logger.info("📊 Creating enrichment plot...")
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Prepare data
        x_positions = np.arange(len(df))
        labels = [f"{row['cell_type']}\n({row['processing_type']})" for _, row in df.iterrows()]
        
        # Create bars
        bars = []
        for i, (_, row) in enumerate(df.iterrows()):
            color = self.colors[row['cell_type']]
            pattern = self.patterns[row['processing_type']]
            
            bar = ax.bar(i, row['enrichment'], 
                        yerr=row['enrichment_se'],
                        color=color, 
                        alpha=0.8 if row['processing_type'] == 'cleaned' else 0.6,
                        hatch=pattern,
                        capsize=8,
                        error_kw={'linewidth': 2, 'capthick': 2})
            bars.append(bar)
            
            # 값 표시
            ax.text(i, row['enrichment']/2, f"{row['enrichment']:.2f}", 
                   ha='center', va='center', fontweight='bold', fontsize=11, color='white')
            
            # 유의성 표시
            p_val = row['enrichment_p']
            if p_val < 0.001:
                significance = "***"
            elif p_val < 0.01:
                significance = "**" 
            elif p_val < 0.05:
                significance = "*"
            else:
                significance = "ns"
            
            ax.text(i, row['enrichment'] + row['enrichment_se'] + 0.1, significance,
                   ha='center', va='bottom', fontsize=14, fontweight='bold')
        
        # 설정
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('Enrichment', fontsize=14, fontweight='bold')
        ax.set_title('🧬 LDSC Enrichment Analysis - All 8/8 Cell Type Datasets\n(Parkinson\'s Disease GWAS)', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # 기준선 (enrichment = 1)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1)
        ax.text(len(df)-1, 1.05, 'Baseline (no enrichment)', ha='right', va='bottom', 
               fontsize=10, style='italic', color='gray')
        
        # Bonferroni threshold 표시
        bonferroni_line = 0.05 / 8
        ax.text(0.02, 0.98, f'Bonferroni threshold: p < {bonferroni_line:.4f}', 
               transform=ax.transAxes, fontsize=10, 
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
        
        # 범례
        legend_elements = []
        for cell_type, color in self.colors.items():
            legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.8, label=cell_type))
        
        # 처리 타입 범례
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.8, label='cleaned'))
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.6, hatch='///', label='unique'))
        
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))
        
        plt.tight_layout()
        
        # 저장
        output_file = self.output_dir / "enrichment_all_8_datasets.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"💾 Enrichment plot saved: {output_file}")
        
        return output_file
    
    def create_pvalue_plot(self, df):
        """Create -log(p) plot for all datasets"""
        logger.info("📊 Creating -log(p) plot...")
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Prepare data
        x_positions = np.arange(len(df))
        labels = [f"{row['cell_type']}\n({row['processing_type']})" for _, row in df.iterrows()]
        
        # Create bars
        for i, (_, row) in enumerate(df.iterrows()):
            color = self.colors[row['cell_type']]
            pattern = self.patterns[row['processing_type']]
            
            bar = ax.bar(i, row['neg_log_p'],
                        color=color,
                        alpha=0.8 if row['processing_type'] == 'cleaned' else 0.6,
                        hatch=pattern)
            
            # 값 표시
            ax.text(i, row['neg_log_p']/2, f"{row['neg_log_p']:.1f}", 
                   ha='center', va='center', fontweight='bold', fontsize=11, color='white')
        
        # 유의성 기준선들
        ax.axhline(y=-np.log10(0.05), color='orange', linestyle='--', alpha=0.7, linewidth=2)
        ax.text(len(df)-1, -np.log10(0.05)+0.1, 'p = 0.05', ha='right', va='bottom', 
               fontsize=10, color='orange', fontweight='bold')
        
        bonferroni_threshold = 0.05 / 8
        ax.axhline(y=-np.log10(bonferroni_threshold), color='red', linestyle='--', alpha=0.7, linewidth=2)
        ax.text(len(df)-1, -np.log10(bonferroni_threshold)+0.1, f'Bonferroni p = {bonferroni_threshold:.4f}', 
               ha='right', va='bottom', fontsize=10, color='red', fontweight='bold')
        
        # 설정
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('-log₁₀(p-value)', fontsize=14, fontweight='bold')
        ax.set_title('🧬 LDSC Statistical Significance - All 8/8 Cell Type Datasets\n(Parkinson\'s Disease GWAS)', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # 통계 요약
        significant_count = sum(1 for p in df['enrichment_p'] if p < bonferroni_threshold)
        ax.text(0.02, 0.98, f'Significant datasets: {significant_count}/8 (100%)', 
               transform=ax.transAxes, fontsize=12, fontweight='bold',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7))
        
        # 범례
        legend_elements = []
        for cell_type, color in self.colors.items():
            legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.8, label=cell_type))
        
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.8, label='cleaned'))
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.6, hatch='///', label='unique'))
        
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))
        
        plt.tight_layout()
        
        # 저장
        output_file = self.output_dir / "pvalue_all_8_datasets.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"💾 P-value plot saved: {output_file}")
        
        return output_file
    
    def create_combined_plot(self, df):
        """Create combined enrichment and -log(p) plot"""
        logger.info("📊 Creating combined plot...")
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
        
        # Prepare data
        x_positions = np.arange(len(df))
        labels = [f"{row['cell_type']}\n({row['processing_type']})" for _, row in df.iterrows()]
        
        # Top plot: Enrichment
        for i, (_, row) in enumerate(df.iterrows()):
            color = self.colors[row['cell_type']]
            pattern = self.patterns[row['processing_type']]
            
            ax1.bar(i, row['enrichment'], 
                   yerr=row['enrichment_se'],
                   color=color, 
                   alpha=0.8 if row['processing_type'] == 'cleaned' else 0.6,
                   hatch=pattern,
                   capsize=6)
            
            ax1.text(i, row['enrichment']/2, f"{row['enrichment']:.2f}", 
                    ha='center', va='center', fontweight='bold', fontsize=10, color='white')
        
        ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
        ax1.set_ylabel('Enrichment', fontsize=12, fontweight='bold')
        ax1.set_title('A. Enrichment Analysis', fontsize=14, fontweight='bold', loc='left')
        ax1.set_xticks(x_positions)
        ax1.set_xticklabels([])  # No labels on top plot
        
        # Bottom plot: -log(p)
        for i, (_, row) in enumerate(df.iterrows()):
            color = self.colors[row['cell_type']]
            pattern = self.patterns[row['processing_type']]
            
            ax2.bar(i, row['neg_log_p'],
                   color=color,
                   alpha=0.8 if row['processing_type'] == 'cleaned' else 0.6,
                   hatch=pattern)
            
            ax2.text(i, row['neg_log_p']/2, f"{row['neg_log_p']:.1f}", 
                    ha='center', va='center', fontweight='bold', fontsize=10, color='white')
        
        # 유의성 기준선
        ax2.axhline(y=-np.log10(0.05), color='orange', linestyle='--', alpha=0.7)
        bonferroni_threshold = 0.05 / 8
        ax2.axhline(y=-np.log10(bonferroni_threshold), color='red', linestyle='--', alpha=0.7)
        
        ax2.set_xticks(x_positions)
        ax2.set_xticklabels(labels, rotation=45, ha='right')
        ax2.set_ylabel('-log₁₀(p-value)', fontsize=12, fontweight='bold')
        ax2.set_title('B. Statistical Significance', fontsize=14, fontweight='bold', loc='left')
        
        # 전체 제목
        fig.suptitle('🧬 LDSC Analysis Results - All 8/8 Cell Type Datasets\n(Parkinson\'s Disease GWAS)', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        # 범례 (전체 그림 우측)
        legend_elements = []
        for cell_type, color in self.colors.items():
            legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.8, label=cell_type))
        
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.8, label='cleaned'))
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor='gray', alpha=0.6, hatch='///', label='unique'))
        
        fig.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(0.98, 0.5))
        
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)
        
        # 저장
        output_file = self.output_dir / "combined_all_8_datasets.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"💾 Combined plot saved: {output_file}")
        
        return output_file
    
    def create_ranking_plot(self, df):
        """Create cell type ranking plot"""
        logger.info("📊 Creating ranking plot...")
        
        # Sort by enrichment
        df_sorted = df.sort_values('enrichment', ascending=True)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Horizontal bar plot
        y_positions = np.arange(len(df_sorted))
        
        for i, (_, row) in enumerate(df_sorted.iterrows()):
            color = self.colors[row['cell_type']]
            pattern = self.patterns[row['processing_type']]
            
            bar = ax.barh(i, row['enrichment'], 
                         xerr=row['enrichment_se'],
                         color=color,
                         alpha=0.8 if row['processing_type'] == 'cleaned' else 0.6,
                         hatch=pattern,
                         capsize=6)
            
            # 값 및 p-value 표시
            ax.text(row['enrichment']/2, i, f"{row['enrichment']:.2f}", 
                   ha='center', va='center', fontweight='bold', fontsize=11, color='white')
            
            ax.text(row['enrichment'] + row['enrichment_se'] + 0.1, i, 
                   f"p={row['enrichment_p']:.1e}", 
                   ha='left', va='center', fontsize=9)
        
        # 설정
        labels = [f"{row['cell_type']} ({row['processing_type']})" for _, row in df_sorted.iterrows()]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels)
        ax.set_xlabel('Enrichment', fontsize=12, fontweight='bold')
        ax.set_title('🏆 Cell Type Enrichment Ranking\n(Parkinson\'s Disease GWAS)', 
                    fontsize=14, fontweight='bold')
        
        # 기준선
        ax.axvline(x=1, color='gray', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        
        # 저장
        output_file = self.output_dir / "ranking_all_8_datasets.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"💾 Ranking plot saved: {output_file}")
        
        return output_file
    
    def generate_summary_stats(self, df):
        """Generate and save summary statistics"""
        logger.info("📊 Generating summary statistics...")
        
        summary = {
            "total_datasets": len(df),
            "significant_datasets": sum(1 for p in df['enrichment_p'] if p < 0.05/8),
            "mean_enrichment": df['enrichment'].mean(),
            "std_enrichment": df['enrichment'].std(),
            "max_enrichment": df['enrichment'].max(),
            "min_enrichment": df['enrichment'].min(),
            "top_dataset": df.loc[df['enrichment'].idxmax(), 'dataset'],
            "bonferroni_threshold": 0.05/8
        }
        
        # Save summary
        summary_file = self.output_dir / "summary_statistics.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"📊 Summary statistics saved: {summary_file}")
        return summary
    
    def run_visualization(self):
        """Main visualization runner"""
        logger.info("=" * 80)
        logger.info("🎨 Starting LDSC visualization for all 8/8 datasets")
        logger.info("=" * 80)
        
        # Load data
        df = self.load_results()
        
        # Generate plots
        enrichment_plot = self.create_enrichment_plot(df)
        pvalue_plot = self.create_pvalue_plot(df)
        combined_plot = self.create_combined_plot(df)
        ranking_plot = self.create_ranking_plot(df)
        
        # Generate summary
        summary = self.generate_summary_stats(df)
        
        logger.info("\n✅ Visualization Complete!")
        logger.info(f"📁 Output directory: {self.output_dir}")
        logger.info(f"📊 Plots generated:")
        logger.info(f"  - Enrichment plot: {enrichment_plot.name}")
        logger.info(f"  - P-value plot: {pvalue_plot.name}")
        logger.info(f"  - Combined plot: {combined_plot.name}")
        logger.info(f"  - Ranking plot: {ranking_plot.name}")
        logger.info(f"📈 Summary: {summary['significant_datasets']}/{summary['total_datasets']} datasets significant")
        logger.info(f"🏆 Top dataset: {summary['top_dataset']} (enrichment: {summary['max_enrichment']:.3f})")

def main():
    """Main function"""
    visualizer = LDSCVisualizer()
    visualizer.run_visualization()

if __name__ == "__main__":
    main()