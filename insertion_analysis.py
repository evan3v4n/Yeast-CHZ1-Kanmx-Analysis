#!/usr/bin/env python3
"""
Bioinformatics analysis for yeast CHZ1 gene knockout with kanamycin insertion
"""

import re
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np

def parse_fastq(filename):
    """Parse FASTQ file and return sequences"""
    sequences = []
    qualities = []
    headers = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for i in range(0, len(lines), 4):
        if i + 3 < len(lines):
            headers.append(lines[i].strip())
            sequences.append(lines[i+1].strip())
            qualities.append(lines[i+3].strip())
    
    return headers, sequences, qualities

def reverse_complement(seq):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def find_alignment_positions(reference, reads, min_match_length=20):
    """Find where reads align to reference sequence"""
    alignments = []
    
    for i, read in enumerate(reads):
        # Try forward orientation
        for j in range(len(reference) - min_match_length + 1):
            if reference[j:j+len(read)].startswith(read[:min_match_length]):
                # Check full alignment
                matches = sum(1 for k in range(min(len(read), len(reference)-j)) 
                             if read[k] == reference[j+k])
                if matches >= len(read) * 0.8:  # 80% similarity threshold
                    alignments.append({
                        'read_index': i,
                        'ref_start': j,
                        'ref_end': j + len(read),
                        'orientation': 'forward',
                        'matches': matches,
                        'length': len(read)
                    })
                    break
        
        # Try reverse complement
        rev_read = reverse_complement(read)
        for j in range(len(reference) - min_match_length + 1):
            if reference[j:j+len(rev_read)].startswith(rev_read[:min_match_length]):
                matches = sum(1 for k in range(min(len(rev_read), len(reference)-j)) 
                             if rev_read[k] == reference[j+k])
                if matches >= len(rev_read) * 0.8:
                    alignments.append({
                        'read_index': i,
                        'ref_start': j,
                        'ref_end': j + len(rev_read),
                        'orientation': 'reverse',
                        'matches': matches,
                        'length': len(rev_read)
                    })
                    break
    
    return alignments

def analyze_coverage(reference, alignments, read_length=150):
    """Analyze coverage across the reference sequence"""
    coverage = np.zeros(len(reference))
    
    for align in alignments:
        start = align['ref_start']
        end = min(align['ref_end'], len(reference))
        coverage[start:end] += 1
    
    return coverage

def identify_construct_features(reference):
    """Identify key features in the construct"""
    features = {}
    
    # Common restriction sites
    restriction_sites = {
        'BamHI': 'GGATCC',
        'SalI': 'GTCGAC',
        'AscI': 'GGCGCGCC',
        'PacI': 'TTAATTAA',
        'ApaI': 'GGGCCC'
    }
    
    for enzyme, site in restriction_sites.items():
        positions = [m.start() for m in re.finditer(site, reference)]
        if positions:
            features[enzyme] = positions
    
    # Look for kanamycin resistance gene signature
    # Common kanamycin promoter sequence
    kan_signatures = [
        'GAGGCCCAGAATACCCTCCTTGAC',  # Part of kanamycin promoter
        'CGTACGCTGCAGGTCGAC',         # Common in yeast knockout cassettes
    ]
    
    for i, sig in enumerate(kan_signatures):
        pos = reference.find(sig)
        if pos != -1:
            features[f'Kan_signature_{i+1}'] = [pos]
    
    return features

def plot_coverage(coverage, features, reference_length, output_file='coverage_plot.png'):
    """Plot coverage across the reference sequence"""
    plt.figure(figsize=(15, 8))
    
    # Plot coverage
    plt.subplot(2, 1, 1)
    positions = np.arange(len(coverage))
    plt.plot(positions, coverage, 'b-', linewidth=0.5)
    plt.fill_between(positions, coverage, alpha=0.3)
    plt.xlabel('Position in Reference')
    plt.ylabel('Coverage Depth')
    plt.title('Sequencing Coverage Across Construct')
    plt.grid(True, alpha=0.3)
    
    # Mark features
    ax = plt.gca()
    y_max = max(coverage) if max(coverage) > 0 else 1
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(features)))
    for i, (feature, positions) in enumerate(features.items()):
        for pos in positions:
            ax.axvline(x=pos, color=colors[i], linestyle='--', alpha=0.7, label=feature)
    
    # Remove duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    # Plot coverage histogram
    plt.subplot(2, 1, 2)
    plt.hist(coverage[coverage > 0], bins=50, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel('Coverage Depth')
    plt.ylabel('Number of Positions')
    plt.title('Coverage Depth Distribution')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def analyze_insertion_junction(reference, alignments, headers, sequences):
    """Analyze reads spanning the insertion junction"""
    # The junction would be where yeast genomic DNA meets the kanamycin cassette
    # Look for reads that don't fully align - these might span junctions
    
    partial_alignments = []
    for align in alignments:
        if align['matches'] < align['length'] * 0.95:  # Less than 95% match
            partial_alignments.append(align)
    
    return partial_alignments

def main():
    # Load reference sequence from the provided data
    reference_lines = [
        "GACATGGAGGCCCAGAATACCC",
        "ATTACAACACTGCTTCACCAACATCTGCATCGTACGCTGCAGGTCGACGGATCCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCCCAGAATACCCTCCTTGAC",
        # Add more lines as needed...
    ]
    
    # Find the longest continuous sequence
    reference = max(reference_lines, key=len)
    print(f"Reference sequence length: {len(reference)} bp")
    
    # Analyze construct features
    print("\nAnalyzing construct features...")
    features = identify_construct_features(reference)
    for feature, positions in features.items():
        print(f"  {feature}: {positions}")
    
    # Instructions for running with actual FASTQ file
    print("\n" + "="*60)
    print("INSTRUCTIONS FOR RUNNING WITH YOUR DATA:")
    print("="*60)
    print("\n1. Save this script as 'analyze_kan_insertion.py'")
    print("\n2. Run the following command:")
    print("   python analyze_kan_insertion.py matching1_2(3).fq")
    print("\n3. The script will generate:")
    print("   - Coverage plot showing sequencing depth")
    print("   - Analysis of insertion junctions")
    print("   - Summary statistics")
    
    # If running with actual file
    import sys
    if len(sys.argv) > 1:
        fastq_file = sys.argv[1]
        print(f"\nAnalyzing {fastq_file}...")
        
        headers, sequences, qualities = parse_fastq(fastq_file)
        print(f"Loaded {len(sequences)} reads")
        
        # Find alignments
        print("Finding alignments (this may take a while)...")
        alignments = find_alignment_positions(reference, sequences[:1000])  # First 1000 reads for speed
        print(f"Found {len(alignments)} alignments")
        
        # Analyze coverage
        coverage = analyze_coverage(reference, alignments)
        mean_coverage = np.mean(coverage[coverage > 0]) if any(coverage > 0) else 0
        print(f"Mean coverage: {mean_coverage:.1f}x")
        
        # Plot results
        plot_coverage(coverage, features, len(reference))
        print("Coverage plot saved as 'coverage_plot.png'")
        
        # Analyze junctions
        junction_reads = analyze_insertion_junction(reference, alignments, headers, sequences)
        print(f"\nFound {len(junction_reads)} potential junction-spanning reads")

if __name__ == "__main__":
    main()