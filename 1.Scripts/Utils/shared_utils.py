#!/usr/bin/env python3
"""
공통 유틸리티 클래스
=================
파킨슨병 GWAS 분석을 위한 공통 데이터 로딩 및 처리 유틸리티
"""

import pandas as pd
import numpy as np
from pathlib import Path
import pickle
from typing import Tuple, Dict, Any


class DataManager:
    """데이터 로딩 및 전처리 관리 클래스"""
    
    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        self.data_dir = self.base_dir / "0_data"
        self.cache_dir = self.base_dir / "shared_cache"
        self.cache_dir.mkdir(exist_ok=True)
        
        # 파일 경로
        self.gwas_file = self.data_dir / "raw" / "GCST009325.h.tsv.gz"
        self.enhancer_file = self.data_dir / "processed" / "Olig_cleaned_hg19_final_sorted.bed"
        
        # 캐시된 데이터
        self._gwas_df = None
        self._enhancers_df = None
        self._classified_gwas_df = None
        
    def load_gwas_data(self, force_reload: bool = False) -> pd.DataFrame:
        """GWAS 데이터 로딩 및 전처리"""
        cache_file = self.cache_dir / "gwas_processed.pkl"
        
        if not force_reload and cache_file.exists() and self._gwas_df is None:
            print("📁 캐시된 GWAS 데이터 로딩...")
            with open(cache_file, 'rb') as f:
                self._gwas_df = pickle.load(f)
            return self._gwas_df
        
        if self._gwas_df is not None and not force_reload:
            return self._gwas_df
            
        print("📊 GWAS 데이터 로딩 및 전처리...")
        
        # 원본 데이터 로딩
        gwas_df = pd.read_csv(self.gwas_file, sep='\t', compression='gzip')
        print(f"   원본 SNPs: {len(gwas_df):,}")
        
        # 데이터 정리
        gwas_df = gwas_df.dropna(subset=['p_value', 'chromosome', 'base_pair_location'])
        gwas_df = gwas_df[(gwas_df['p_value'] > 0) & (gwas_df['p_value'] <= 1)]
        gwas_df = gwas_df[gwas_df['chromosome'].isin(range(1, 23))]
        
        # 중복 제거 (동일 위치에서 가장 유의한 것만 유지)
        gwas_df = gwas_df.sort_values('p_value').drop_duplicates(
            subset=['chromosome', 'base_pair_location'], keep='first')
        
        # -log10(p) 계산
        gwas_df['neg_log10_p'] = -np.log10(gwas_df['p_value'] + 1e-300)
        
        print(f"   정리된 SNPs: {len(gwas_df):,}")
        
        # 캐시 저장
        with open(cache_file, 'wb') as f:
            pickle.dump(gwas_df, f)
        
        self._gwas_df = gwas_df
        return gwas_df
    
    def load_enhancer_data(self, force_reload: bool = False) -> pd.DataFrame:
        """Enhancer 데이터 로딩 및 전처리 (좌표계 변환 포함)"""
        cache_file = self.cache_dir / "enhancers_processed.pkl"
        
        if not force_reload and cache_file.exists() and self._enhancers_df is None:
            print("📁 캐시된 Enhancer 데이터 로딩...")
            with open(cache_file, 'rb') as f:
                self._enhancers_df = pickle.load(f)
            return self._enhancers_df
        
        if self._enhancers_df is not None and not force_reload:
            return self._enhancers_df
            
        print("🎯 Enhancer 데이터 로딩 및 전처리...")
        
        # 좌표 변환된 파일 우선 확인
        converted_file = self._find_converted_enhancer_file()
        
        if converted_file and converted_file.exists():
            print(f"   변환된 좌표 파일 사용: {converted_file}")
            enhancers_df = pd.read_csv(converted_file, sep='\t', header=None, 
                                      names=['CHR', 'START', 'END', 'NAME'])
        else:
            print("   ⚠️  변환된 좌표 파일이 없습니다. 원본 파일 사용 (부정확할 수 있음)")
            print("   정확한 분석을 위해 setup_liftover.py를 먼저 실행하세요.")
            enhancers_df = pd.read_csv(self.enhancer_file, sep='\t', header=None, 
                                      names=['CHR', 'START', 'END', 'NAME'])
        
        # 데이터 정리
        enhancers_df['CHR'] = enhancers_df['CHR'].str.replace('chr', '')
        numeric_mask = enhancers_df['CHR'].str.isnumeric()
        enhancers_df = enhancers_df[numeric_mask].copy()
        enhancers_df['CHR'] = enhancers_df['CHR'].astype(int)
        enhancers_df = enhancers_df[enhancers_df['CHR'].isin(range(1, 23))]
        
        print(f"   Enhancer 영역: {len(enhancers_df):,}")
        
        # 캐시 저장
        with open(cache_file, 'wb') as f:
            pickle.dump(enhancers_df, f)
        
        self._enhancers_df = enhancers_df
        return enhancers_df
    
    def _find_converted_enhancer_file(self) -> Path:
        """변환된 enhancer 파일 찾기"""
        # 변환된 파일 경로 생성
        base_dir = Path(".")
        hg19_dir = base_dir / "0_data" / "processed" / "hg19_coordinates"
        
        if not hg19_dir.exists():
            return None
        
        # 원본 파일 경로에서 변환된 파일명 추정
        # 예: 0_data/raw/cleaned_data/Olig_cleaned.bed → cleaned_data_Olig_cleaned_hg19.bed
        original_path = Path(self.enhancer_file)
        parent_name = original_path.parent.name  # cleaned_data 또는 unique_data
        file_stem = original_path.stem  # Olig_cleaned
        
        converted_filename = f"{parent_name}_{file_stem}_hg19.bed"
        converted_path = hg19_dir / converted_filename
        
        return converted_path if converted_path.exists() else None
    
    def classify_snps(self, force_reload: bool = False) -> pd.DataFrame:
        """SNP을 enhancer 내/외부로 분류"""
        cache_file = self.cache_dir / "classified_gwas.pkl"
        
        if not force_reload and cache_file.exists() and self._classified_gwas_df is None:
            print("📁 캐시된 분류 데이터 로딩...")
            with open(cache_file, 'rb') as f:
                self._classified_gwas_df = pickle.load(f)
            return self._classified_gwas_df
        
        if self._classified_gwas_df is not None and not force_reload:
            return self._classified_gwas_df
            
        print("🎯 SNP 분류 진행...")
        
        # 데이터 로딩
        gwas_df = self.load_gwas_data().copy()
        enhancers_df = self.load_enhancer_data()
        
        # 분류 초기화
        gwas_df['in_enhancer'] = False
        
        # Enhancer 영역과 교집합 찾기
        for _, enhancer in enhancers_df.iterrows():
            mask = ((gwas_df['chromosome'] == enhancer['CHR']) &
                   (gwas_df['base_pair_location'] >= enhancer['START']) &
                   (gwas_df['base_pair_location'] <= enhancer['END']))
            gwas_df.loc[mask, 'in_enhancer'] = True
        
        enhancer_count = gwas_df['in_enhancer'].sum()
        total_count = len(gwas_df)
        
        print(f"   Enhancer 내 SNPs: {enhancer_count:,} ({enhancer_count/total_count:.4f})")
        print(f"   Background SNPs: {total_count-enhancer_count:,}")
        
        # 캐시 저장
        with open(cache_file, 'wb') as f:
            pickle.dump(gwas_df, f)
        
        self._classified_gwas_df = gwas_df
        return gwas_df
    
    def get_all_data(self, force_reload: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """모든 데이터 반환 (GWAS, Enhancer, 분류된 GWAS)"""
        gwas_df = self.load_gwas_data(force_reload)
        enhancers_df = self.load_enhancer_data(force_reload)
        classified_gwas_df = self.classify_snps(force_reload)
        
        return gwas_df, enhancers_df, classified_gwas_df
    
    def clear_cache(self):
        """캐시 파일 삭제"""
        import shutil
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
            self.cache_dir.mkdir(exist_ok=True)
        print("🗑️ 캐시 초기화 완료")


class StatisticalAnalyzer:
    """통계 분석 유틸리티 클래스"""
    
    @staticmethod
    def calculate_enrichment_stats(classified_gwas_df: pd.DataFrame, 
                                 significance_threshold: float = 5e-8) -> Dict[str, Any]:
        """Enrichment 통계 계산 (LDSC 스타일 포함)"""
        from scipy.stats import mannwhitneyu, fisher_exact
        import numpy as np
        
        enhancer_snps = classified_gwas_df[classified_gwas_df['in_enhancer']]
        non_enhancer_snps = classified_gwas_df[~classified_gwas_df['in_enhancer']]
        
        # P-value 분포 비교 (Mann-Whitney U test)
        enhancer_neg_log_p = enhancer_snps['neg_log10_p']
        non_enhancer_neg_log_p = non_enhancer_snps['neg_log10_p']
        
        mw_stat, mw_pval = mannwhitneyu(enhancer_neg_log_p, non_enhancer_neg_log_p, 
                                       alternative='greater')
        
        # 유의한 SNP 비율 비교 (Fisher's exact test)
        enhancer_sig = (enhancer_snps['p_value'] < significance_threshold).sum()
        enhancer_total = len(enhancer_snps)
        non_enhancer_sig = (non_enhancer_snps['p_value'] < significance_threshold).sum()
        non_enhancer_total = len(non_enhancer_snps)
        
        contingency = [[enhancer_sig, enhancer_total - enhancer_sig],
                      [non_enhancer_sig, non_enhancer_total - non_enhancer_sig]]
        odds_ratio, fisher_pval = fisher_exact(contingency, alternative='greater')
        
        # Enrichment ratio (기존 방식)
        enhancer_sig_rate = enhancer_sig / enhancer_total if enhancer_total > 0 else 0
        non_enhancer_sig_rate = non_enhancer_sig / non_enhancer_total if non_enhancer_total > 0 else 0
        enrichment_ratio = enhancer_sig_rate / non_enhancer_sig_rate if non_enhancer_sig_rate > 0 else np.inf
        
        # LDSC 스타일 enrichment 계산
        ldsc_enrichment = StatisticalAnalyzer.calculate_ldsc_enrichment(classified_gwas_df)
        
        return {
            'total_snps': len(classified_gwas_df),
            'enhancer_snps': enhancer_total,
            'enhancer_sig_snps': enhancer_sig,
            'background_sig_snps': non_enhancer_sig,
            'enhancer_sig_rate': enhancer_sig_rate,
            'background_sig_rate': non_enhancer_sig_rate,
            'enrichment_ratio': enrichment_ratio,
            'ldsc_enrichment': ldsc_enrichment['enrichment'],
            'ldsc_enrichment_se': ldsc_enrichment['enrichment_se'],
            'ldsc_enrichment_p': ldsc_enrichment['enrichment_p'],
            'mann_whitney_pval': mw_pval,
            'fisher_pval': fisher_pval,
            'odds_ratio': odds_ratio,
            'mann_whitney_stat': mw_stat,
            'significance_threshold': significance_threshold
        }
    
    @staticmethod
    def calculate_ldsc_enrichment(classified_gwas_df: pd.DataFrame) -> Dict[str, float]:
        """LDSC 스타일 enrichment 계산"""
        
        enhancer_snps = classified_gwas_df[classified_gwas_df['in_enhancer']]
        all_snps = classified_gwas_df
        
        # 효과적인 annotation 수 계산
        M_annot = len(enhancer_snps)  # annotation에 속한 SNP 수
        M_total = len(all_snps)       # 전체 SNP 수
        
        # Chi-square 통계량 계산 (z-score의 제곱)
        enhancer_chisq = enhancer_snps['neg_log10_p'].apply(lambda x: (x / np.log10(np.e))**2 if x > 0 else 0)
        all_chisq = all_snps['neg_log10_p'].apply(lambda x: (x / np.log10(np.e))**2 if x > 0 else 0)
        
        # LDSC regression 근사 계산
        sum_chisq_annot = enhancer_chisq.sum()
        sum_chisq_total = all_chisq.sum()
        
        # Enrichment 계산 (annotation당 평균 chi-square / 전체 평균 chi-square)
        mean_chisq_annot = sum_chisq_annot / M_annot if M_annot > 0 else 0
        mean_chisq_total = sum_chisq_total / M_total if M_total > 0 else 0
        
        ldsc_enrichment = mean_chisq_annot / mean_chisq_total if mean_chisq_total > 0 else 0
        
        # Standard error 근사 계산
        var_chisq_annot = enhancer_chisq.var() if M_annot > 1 else 0
        var_chisq_total = all_chisq.var() if M_total > 1 else 0
        
        # SE 계산 (delta method 근사)
        se_numerator = np.sqrt(var_chisq_annot / M_annot) if M_annot > 0 else 0
        se_denominator = np.sqrt(var_chisq_total / M_total) if M_total > 0 else 0
        
        if mean_chisq_total > 0 and M_annot > 0:
            enrichment_se = ldsc_enrichment * np.sqrt(
                (se_numerator / mean_chisq_annot)**2 + (se_denominator / mean_chisq_total)**2
            )
        else:
            enrichment_se = 0
        
        # Z-score 및 p-value 계산
        if enrichment_se > 0:
            z_score = (ldsc_enrichment - 1) / enrichment_se
            from scipy.stats import norm
            enrichment_p = 2 * (1 - norm.cdf(abs(z_score)))  # two-tailed test
        else:
            enrichment_p = 1.0
        
        return {
            'enrichment': ldsc_enrichment,
            'enrichment_se': enrichment_se,
            'enrichment_p': enrichment_p,
            'mean_chisq_annot': mean_chisq_annot,
            'mean_chisq_total': mean_chisq_total,
            'n_snps_annot': M_annot,
            'n_snps_total': M_total
        }


class ManhattanPlotData:
    """Manhattan plot용 데이터 준비 클래스"""
    
    @staticmethod
    def prepare_plot_data(classified_gwas_df: pd.DataFrame, 
                         max_points: int = 200000) -> Dict[str, Any]:
        """Manhattan plot용 데이터 준비"""
        print("📊 Manhattan plot 데이터 준비...")
        
        # 염색체별로 누적 위치 계산
        plot_data = []
        cumulative_pos = 0
        chr_centers = {}
        chr_colors = {}
        
        # 염색체별 색상 (번갈아가며)
        colors = ['#4472C4', '#70AD47']  # 파란색, 녹색 번갈아
        
        for chrom in range(1, 23):
            chr_data = classified_gwas_df[classified_gwas_df['chromosome'] == chrom].copy()
            if len(chr_data) == 0:
                continue
            
            # 위치별 정렬
            chr_data = chr_data.sort_values('base_pair_location')
            
            # 플롯 위치 계산
            max_pos = chr_data['base_pair_location'].max()
            chr_data['plot_pos'] = chr_data['base_pair_location'] + cumulative_pos
            
            # 염색체 중심 위치 저장 (라벨용)
            chr_centers[chrom] = cumulative_pos + max_pos / 2
            chr_colors[chrom] = colors[(chrom - 1) % 2]
            
            cumulative_pos += max_pos + 5e6  # 5Mb 간격
            plot_data.append(chr_data)
        
        # 전체 데이터 결합
        plot_df = pd.concat(plot_data, ignore_index=True)
        
        # 너무 많은 점이면 샘플링 (시각화 성능을 위해)
        if len(plot_df) > max_points:
            # 유의한 SNP는 모두 유지, 나머지는 샘플링
            significant = plot_df[plot_df['neg_log10_p'] > -np.log10(1e-4)]
            non_significant = plot_df[plot_df['neg_log10_p'] <= -np.log10(1e-4)]
            
            remaining_points = max_points - len(significant)
            if len(non_significant) > remaining_points and remaining_points > 0:
                sampled_non_sig = non_significant.sample(remaining_points, random_state=42)
            else:
                sampled_non_sig = non_significant
            
            plot_df = pd.concat([significant, sampled_non_sig], ignore_index=True)
            print(f"   📊 시각화용 샘플링: {len(plot_df):,} SNPs")
        
        return {
            'plot_df': plot_df,
            'chr_centers': chr_centers,
            'chr_colors': chr_colors
        }


class ResultsManager:
    """결과 저장 및 로딩 관리 클래스"""
    
    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        self.results_dir = self.base_dir / "2_analysis" / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def save_enrichment_results(self, results: Dict[str, Any]) -> Path:
        """Enrichment 분석 결과 저장"""
        results_df = pd.DataFrame([results])
        output_file = self.results_dir / "enrichment_results.csv"
        results_df.to_csv(output_file, index=False)
        
        print(f"   💾 결과 저장: {output_file}")
        return output_file
    
    def load_enrichment_results(self) -> Dict[str, Any]:
        """Enrichment 분석 결과 로딩"""
        results_file = self.results_dir / "enrichment_results.csv"
        if not results_file.exists():
            raise FileNotFoundError("결과 파일이 없습니다. 먼저 분석을 실행하세요.")
        
        results_df = pd.read_csv(results_file)
        return results_df.iloc[0].to_dict()
    
    def save_detailed_snp_data(self, classified_gwas_df: pd.DataFrame):
        """상세한 SNP 데이터 저장 (시각화용)"""
        # Genome-wide significant SNPs 저장
        significant_snps = classified_gwas_df[classified_gwas_df['p_value'] < 5e-8]
        sig_file = self.results_dir / "genome_wide_significant_snps.csv"
        significant_snps.to_csv(sig_file, index=False)
        
        # Enhancer 내 유의한 SNPs 저장
        enhancer_sig_snps = significant_snps[significant_snps['in_enhancer']]
        enh_sig_file = self.results_dir / "enhancer_significant_snps.csv"
        enhancer_sig_snps.to_csv(enh_sig_file, index=False)
        
        print(f"   💾 상세 SNP 데이터 저장 완료")
        print(f"      - 유의한 SNPs: {len(significant_snps)}")
        print(f"      - Enhancer 내 유의한 SNPs: {len(enhancer_sig_snps)}")