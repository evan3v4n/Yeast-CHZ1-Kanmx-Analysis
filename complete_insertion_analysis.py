#!/usr/bin/env python3
"""
Complete analysis pipeline for CHZ1 kanamycin insertion project
"""

import os
import re
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# First, let's identify what we're looking for
EXPECTED_FEATURES = {
    'CHZ1_upstream': 'ATTACAACACTGCTTCACCAACATCTGCATC',  # Part of CHZ1 sequence
    'Kanamycin_junction': 'GTACGCTGCAGGTCGAC',  # Junction between CHZ1 and kan
    'Restriction_sites': {
        'BamHI': 'GGATCC',
        'SalI': 'GTCGAC',
        'PacI': 'TTAATTAA',
        'AscI': 'GGCGCGCC',
    },
    'Kan_promoter': 'GAGGCCCAGAATACCCTCCTTGAC',  # Kanamycin promoter region
}

def parse_reference_sequence():
    """Parse the reference sequence from the provided data"""
    # This is your construct sequence
    ref_seq = "GACATGGAGGCCCAGAATACCCTCCTTGACATTACAACACTGCTTCACCAACATCTGCATCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGAC"
    return ref_seq

def quick_fastq_analysis(filename, num_reads=1000):
    """Quick analysis of FASTQ file"""
    read_lengths = []
    quality_scores = []
    base_composition = Counter()
    
    with open(filename, 'r') as f:
        line_count = 0
        for line in f:
            if line_count % 4 == 1:  # Sequence line
                seq = line.strip()
                read_lengths.append(len(seq))
                base_composition.update(seq)
            elif line_count % 4 == 3:  # Quality line
                qual = line.strip()
                # Convert quality scores
                scores = [ord(c) - 33 for c in qual]  # Assuming Phred+33
                quality_scores.extend(scores)
            
            line_count += 1
            if line_count >= num_reads * 4:
                break
    
    return {
        'total_reads': line_count // 4,
        'mean_length': np.mean(read_lengths),
        'mean_quality': np.mean(quality_scores),
        'base_composition': dict(base_composition)
    }

def extract_junction_reads(filename, construct_seq, window=30):
    """Extract reads that might span the insertion junction"""
    junction_reads = []
    
    # Key sequences to look for
    chz1_end = "ATCTGCATC"  # End of CHZ1 sequence
    kan_start = "GTACGCTGCAGGTCGAC"  # Start of kanamycin cassette
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i in range(0, len(lines), 4):
        if i + 1 < len(lines):
            header = lines[i].strip()
            seq = lines[i+1].strip()
            
            # Check if read contains junction signatures
            if chz1_end in seq or kan_start in seq:
                if i + 3 < len(lines):
                    quality = lines[i+3].strip()
                    junction_reads.append({
                        'header': header,
                        'sequence': seq,
                        'quality': quality,
                        'has_chz1': chz1_end in seq,
                        'has_kan': kan_start in seq
                    })
    
    return junction_reads

def analyze_construct_coverage(filename, construct_seq, output_prefix="analysis"):
    """Comprehensive analysis of sequencing coverage on construct"""
    
    # Initialize coverage array
    coverage = np.zeros(len(construct_seq))
    mismatch_positions = defaultdict(Counter)
    
    reads_analyzed = 0
    reads_mapped = 0
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    print("Analyzing reads...")
    for i in range(0, min(len(lines), 10000), 4):  # Analyze first 2500 reads
        if i + 1 < len(lines):
            seq = lines[i+1].strip()
            
            # Try to find this read in the construct
            for start_pos in range(len(construct_seq) - len(seq) + 1):
                if construct_seq[start_pos:start_pos+20] == seq[:20]:  # Quick check
                    # Detailed alignment
                    matches = 0
                    for j in range(min(len(seq), len(construct_seq) - start_pos)):
                        if seq[j] == construct_seq[start_pos + j]:
                            matches += 1
                            coverage[start_pos + j] += 1
                        else:
                            mismatch_positions[start_pos + j][seq[j]] += 1
                    
                    if matches > len(seq) * 0.8:  # 80% match threshold
                        reads_mapped += 1
                        break
            
            # Also try reverse complement
            seq_rc = str(Seq(seq).reverse_complement())
            for start_pos in range(len(construct_seq) - len(seq_rc) + 1):
                if construct_seq[start_pos:start_pos+20] == seq_rc[:20]:
                    matches = 0
                    for j in range(min(len(seq_rc), len(construct_seq) - start_pos)):
                        if seq_rc[j] == construct_seq[start_pos + j]:
                            matches += 1
                            coverage[start_pos + j] += 1
                        else:
                            mismatch_positions[start_pos + j][seq_rc[j]] += 1
                    
                    if matches > len(seq_rc) * 0.8:
                        reads_mapped += 1
                        break
            
            reads_analyzed += 1
            if reads_analyzed % 100 == 0:
                print(f"  Analyzed {reads_analyzed} reads, mapped {reads_mapped}")
    
    # Create comprehensive plot
    fig, axes = plt.subplots(3, 1, figsize=(15, 12))
    
    # Plot 1: Coverage depth
    ax1 = axes[0]
    positions = np.arange(len(coverage))
    ax1.plot(positions, coverage, 'b-', linewidth=1)
    ax1.fill_between(positions, coverage, alpha=0.3)
    ax1.set_xlabel('Position in Construct (bp)')
    ax1.set_ylabel('Coverage Depth')
    ax1.set_title(f'Sequencing Coverage Across Kanamycin Insertion Construct\n({reads_mapped}/{reads_analyzed} reads mapped)')
    ax1.grid(True, alpha=0.3)
    
    # Mark important features
    feature_colors = ['red', 'green', 'orange', 'purple']
    feature_positions = []
    
    # Mark CHZ1-Kan junction
    junction_pos = construct_seq.find('GTACGCTGCAGGTCGAC')
    if junction_pos != -1:
        ax1.axvline(x=junction_pos, color='red', linestyle='--', linewidth=2, label='CHZ1-Kan Junction')
        feature_positions.append(('Junction', junction_pos))
    
    # Mark restriction sites
    for enzyme, site in EXPECTED_FEATURES['Restriction_sites'].items():
        pos = construct_seq.find(site)
        if pos != -1:
            ax1.axvline(x=pos, color='green', linestyle=':', alpha=0.7, label=f'{enzyme} site')
    
    ax1.legend(loc='upper right')
    
    # Plot 2: Coverage distribution
    ax2 = axes[1]
    coverage_nonzero = coverage[coverage > 0]
    if len(coverage_nonzero) > 0:
        ax2.hist(coverage_nonzero, bins=50, alpha=0.7, color='blue', edgecolor='black')
        ax2.axvline(x=np.mean(coverage_nonzero), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(coverage_nonzero):.1f}x')
        ax2.axvline(x=np.median(coverage_nonzero), color='green', linestyle='--', 
                   label=f'Median: {np.median(coverage_nonzero):.1f}x')
    ax2.set_xlabel('Coverage Depth')
    ax2.set_ylabel('Number of Positions')
    ax2.set_title('Coverage Depth Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Mismatch frequency
    ax3 = axes[2]
    mismatch_freq = []
    positions_with_mismatches = []
    
    for pos in range(len(construct_seq)):
        if pos in mismatch_positions and coverage[pos] > 0:
            total_mismatches = sum(mismatch_positions[pos].values())
            freq = total_mismatches / (coverage[pos] + total_mismatches)
            if freq > 0.05:  # Only show positions with >5% mismatch
                mismatch_freq.append(freq)
                positions_with_mismatches.append(pos)
    
    if positions_with_mismatches:
        ax3.bar(positions_with_mismatches, mismatch_freq, width=1, color='red', alpha=0.7)
        ax3.set_xlabel('Position in Construct (bp)')
        ax3.set_ylabel('Mismatch Frequency')
        ax3.set_title('Positions with High Mismatch Rates (>5%)')
    else:
        ax3.text(0.5, 0.5, 'No significant mismatches detected', 
                horizontalalignment='center', verticalalignment='center',
                transform=ax3.transAxes, fontsize=14)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_coverage_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Generate summary report
    report = f"""
Kanamycin Insertion Analysis Report
===================================

Construct Length: {len(construct_seq)} bp
Total Reads Analyzed: {reads_analyzed}
Reads Mapped to Construct: {reads_mapped} ({reads_mapped/reads_analyzed*100:.1f}%)

Coverage Statistics:
- Mean Coverage: {np.mean(coverage):.1f}x
- Positions with Coverage: {np.sum(coverage > 0)} ({np.sum(coverage > 0)/len(coverage)*100:.1f}%)
- Maximum Coverage: {np.max(coverage):.0f}x

Key Features Found:
"""
    
    for feature_name, feature_seq in [('CHZ1 region', 'ATCTGCATC'), 
                                      ('Kan junction', 'GTACGCTGCAGGTCGAC'),
                                      ('Kan promoter', 'GAGGCCCAGAATACCCTCCTTGAC')]:
        pos = construct_seq.find(feature_seq)
        if pos != -1:
            report += f"- {feature_name}: Position {pos}, Coverage {coverage[pos:pos+len(feature_seq)].mean():.1f}x\n"
    
    # Check for coverage gaps
    gaps = []
    in_gap = False
    gap_start = 0
    
    for i, cov in enumerate(coverage):
        if cov == 0 and not in_gap:
            in_gap = True
            gap_start = i
        elif cov > 0 and in_gap:
            in_gap = False
            if i - gap_start > 10:  # Only report gaps > 10bp
                gaps.append((gap_start, i))
    
    if gaps:
        report += f"\nCoverage Gaps (>10bp):\n"
        for start, end in gaps:
            report += f"- Position {start}-{end} ({end-start} bp)\n"
    
    with open(f'{output_prefix}_report.txt', 'w') as f:
        f.write(report)
    
    print(report)
    
    return coverage, mismatch_positions, report

def main():
    print("CHZ1 Kanamycin Insertion Analysis Pipeline")
    print("=" * 50)
    
    # Get construct sequence
    construct_seq = parse_reference_sequence()
    print(f"Construct sequence loaded: {len(construct_seq)} bp")
    
    # Check for expected features
    print("\nChecking for expected features in construct:")
    for feature, seq in [('CHZ1 end', 'ATCTGCATC'), 
                        ('Kan start', 'GTACGCTGCAGGTCGAC'),
                        ('BamHI', 'GGATCC'),
                        ('Kan promoter', 'GAGGCCCAGAATACCCTCCTTGAC')]:
        pos = construct_seq.find(seq)
        if pos != -1:
            print(f"  ✓ {feature} found at position {pos}")
        else:
            print(f"  ✗ {feature} not found")
    
    print("\n" + "="*50)
    print("To run the full analysis, execute:")
    print("python this_script.py matching1_2(3).fq")
    print("="*50)
    
    # If FASTQ file is provided as argument
    import sys
    if len(sys.argv) > 1:
        fastq_file = sys.argv[1]
        
        if os.path.exists(fastq_file):
            print(f"\nAnalyzing {fastq_file}...")
            
            # Quick stats
            print("\nGathering FASTQ statistics...")
            stats = quick_fastq_analysis(fastq_file)
            print(f"Total reads: {stats['total_reads']}")
            print(f"Mean read length: {stats['mean_length']:.1f} bp")
            print(f"Mean quality score: {stats['mean_quality']:.1f}")
            
            # Extract junction reads
            print("\nSearching for junction-spanning reads...")
            junction_reads = extract_junction_reads(fastq_file, construct_seq)
            print(f"Found {len(junction_reads)} potential junction reads")
            
            if junction_reads:
                with open('junction_reads.txt', 'w') as f:
                    for read in junction_reads[:10]:  # Save first 10
                        f.write(f"{read['header']}\n")
                        f.write(f"Sequence: {read['sequence']}\n")
                        f.write(f"Has CHZ1: {read['has_chz1']}, Has Kan: {read['has_kan']}\n\n")
                print("Junction reads saved to junction_reads.txt")
            
            # Full coverage analysis
            print("\nPerforming coverage analysis...")
            coverage, mismatches, report = analyze_construct_coverage(fastq_file, construct_seq)
            
            print("\nAnalysis complete! Generated files:")
            print("- analysis_coverage_analysis.png")
            print("- analysis_report.txt")
            print("- junction_reads.txt")
        else:
            print(f"Error: File {fastq_file} not found!")

if __name__ == "__main__":
    main()