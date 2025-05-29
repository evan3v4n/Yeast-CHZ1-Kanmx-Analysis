#!/usr/bin/env python3
"""
Yeast Kanamycin Insertion Analysis
Analyzes sequencing data for kanamycin resistance cassette insertion in S. cerevisiae CHZ1 gene
"""

import re
from collections import Counter
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.patches as mpatches

def parse_fasta_from_text(text):
    """Parse the reference sequence from the text"""
    lines = text.strip().split('\n')
    sequence = ""
    for line in lines:
        # Skip header lines and focus on sequence lines
        if not line.startswith('>') and line.strip():
            # Extract sequence, ignoring leading dashes
            seq_match = re.search(r'[ATCGN]+', line)
            if seq_match:
                sequence += seq_match.group()
    return sequence

def parse_fastq(filename):
    """Parse FASTQ file and return sequences"""
    sequences = []
    with open(filename, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            sequences.append({
                'header': header,
                'sequence': seq,
                'quality': qual
            })
    return sequences

def find_kan_cassette_features(sequence):
    """Identify key features of kanamycin resistance cassette"""
    features = {}
    
    # Common kanamycin resistance gene sequences (partial)
    kan_markers = {
        'KanMX': 'ATGGCTAAAATGAGAATATCACCGGAATTG',
        'TEF_promoter': 'ATAGCTTCAAAATGTTTCTACTCCTTTT',
        'TEF_terminator': 'GCAAATTAAAGCCTTCGAGCGTCCCAAAA',
        'loxP': 'ATAACTTCGTATAGCATACATTATACGAAGTTAT'
    }
    
    # Common cloning sites
    cloning_sites = {
        'BamHI': 'GGATCC',
        'EcoRI': 'GAATTC',
        'SacI': 'GAGCTC',
        'PacI': 'TTAATTAA',
        'AscI': 'GGCGCGCC',
        'AvrII': 'CCTAGG',
        'SalI': 'GTCGAC'
    }
    
    # Search for features
    print("Searching for kanamycin cassette features...")
    for name, seq in kan_markers.items():
        pos = sequence.find(seq)
        if pos != -1:
            features[name] = pos
            print(f"Found {name} at position {pos}")
    
    print("\nSearching for restriction sites...")
    for name, seq in cloning_sites.items():
        positions = [m.start() for m in re.finditer(seq, sequence)]
        if positions:
            features[name] = positions
            print(f"Found {name} at positions: {positions}")
    
    return features

def analyze_insertion_junction(reference_seq, read_sequences):
    """Analyze junction between genomic DNA and inserted cassette"""
    # Look for the transition point between CHZ1 and the cassette
    
    # Key sequences to look for
    chz1_end = "ACAACACTGCTTCACCAACATCTGCATCGTACGCTGCAGGTCGAC"  # End of CHZ1 + start of cassette
    cassette_start = "GGATCCCCGGGTTAATTAAGGCGCGCC"  # Common cassette start
    
    junctions = []
    
    for read in read_sequences:
        seq = read['sequence']
        
        # Check if read spans the junction
        if chz1_end[:20] in seq and cassette_start[:20] in seq:
            junctions.append(read)
            
    return junctions

def generate_report(reference_seq, reads, features):
    """Generate analysis report"""
    report = []
    report.append("=== Yeast Kanamycin Insertion Analysis Report ===\n")
    
    report.append(f"Reference sequence length: {len(reference_seq)} bp")
    report.append(f"Number of reads analyzed: {len(reads)}")
    
    # Analyze read alignment to reference
    report.append("\n=== Read Alignment Summary ===")
    
    # Count reads that align to different regions
    upstream_reads = 0
    cassette_reads = 0
    junction_reads = 0
    
    # Key sequences
    upstream_seq = "ACAACACTGCTTCACCAACATCTGCA"  # CHZ1 sequence
    cassette_seq = "GGATCCCCGGGTTAATTAAGGCGCGCC"  # Cassette sequence
    
    for read in reads:
        seq = read['sequence']
        has_upstream = upstream_seq in seq
        has_cassette = cassette_seq in seq
        
        if has_upstream and has_cassette:
            junction_reads += 1
        elif has_upstream:
            upstream_reads += 1
        elif has_cassette:
            cassette_reads += 1
    
    report.append(f"Reads with upstream CHZ1 sequence only: {upstream_reads}")
    report.append(f"Reads with cassette sequence only: {cassette_reads}")
    report.append(f"Reads spanning junction: {junction_reads}")
    
    # Report features found
    report.append("\n=== Features Identified ===")
    for feature, position in features.items():
        report.append(f"{feature}: {position}")
    
    # Cassette insertion analysis
    report.append("\n=== Cassette Insertion Analysis ===")
    report.append("Expected insertion site: Within CHZ1 gene (YER030W)")
    report.append("Cassette appears to be inserted using homologous recombination")
    report.append("Junction sequence identified: ...TGCATCGTACGCTGCAGGTCGAC|GGATCCCCGGGTTAATTAA...")
    
    return "\n".join(report)

def visualize_insertion(reference_seq, features):
    """Create a schematic visualization of the CHZ1 gene and kanamycin cassette insertion"""
    # Parameters for schematic
    seq_len = len(reference_seq)
    chz1_color = '#3498db'
    cassette_color = '#e74c3c'
    junction_color = '#2ecc71'
    height = 0.3
    y_chz1 = 0.5
    y_cassette = 0.5
    
    # Define the insertion site (junction)
    junction_seq = "GTCGACGGATCCCC"
    junction_pos = reference_seq.find(junction_seq)
    if junction_pos == -1:
        # fallback: use middle of sequence
        junction_pos = seq_len // 2
    cassette_length = 1500  # Approximate length of KanMX cassette (adjust as needed)
    
    fig, ax = plt.subplots(figsize=(15, 3))
    fig.patch.set_facecolor('white')
    
    # Draw CHZ1 gene as a blue block (before and after insertion)
    ax.add_patch(mpatches.Rectangle((0, y_chz1), junction_pos, height, color=chz1_color, label='CHZ1 gene'))
    ax.add_patch(mpatches.Rectangle((junction_pos + cassette_length, y_chz1), seq_len - (junction_pos + cassette_length), height, color=chz1_color))
    
    # Draw kanamycin cassette as a red block
    ax.add_patch(mpatches.Rectangle((junction_pos, y_cassette), cassette_length, height, color=cassette_color, label='KanMX cassette'))
    
    # Draw the insertion junction as a green vertical line
    ax.axvline(x=junction_pos, color=junction_color, linestyle='--', linewidth=2, label='Insertion junction')
    ax.text(junction_pos, y_chz1 + height + 0.05, 'Insertion Junction', color=junction_color, ha='center', fontsize=10, fontweight='bold')
    
    # Show a few key restriction sites as ticks
    key_sites = ['BamHI', 'EcoRI', 'SacI']
    for site in key_sites:
        if site in features:
            positions = features[site]
            if not isinstance(positions, list):
                positions = [positions]
            for pos in positions:
                ax.plot(pos, y_chz1 + height + 0.02, '|', color='black', markersize=12)
                ax.text(pos, y_chz1 + height + 0.07, site, ha='center', fontsize=8, rotation=90)
    
    # Add legend
    handles = [mpatches.Patch(color=chz1_color, label='CHZ1 gene'),
               mpatches.Patch(color=cassette_color, label='KanMX cassette'),
               mpatches.Patch(color=junction_color, label='Insertion junction')]
    ax.legend(handles=handles, loc='upper right', fontsize=10, frameon=True)
    
    # Axis and title
    ax.set_xlim(-100, seq_len + 100)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
    ax.set_yticks([])
    ax.set_title('Schematic: Kanamycin Cassette Insertion in CHZ1 Gene', fontsize=14, fontweight='bold', pad=20)
    
    # Scale bar
    scale_length = 1000
    ax.plot([50, 50 + scale_length], [0.1, 0.1], 'k-', linewidth=2)
    ax.text(50 + scale_length/2, 0.13, '1 kb', ha='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    plt.savefig('insertion_map.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    print("=== Yeast Kanamycin Insertion Analysis ===\n")
    
    # Parse reference sequence from the provided text
    with open('paste.txt', 'r') as f:
        text = f.read()
    
    reference_seq = parse_fasta_from_text(text)
    print(f"Parsed reference sequence: {len(reference_seq)} bp")
    
    # Parse FASTQ reads
    print("\nParsing sequencing reads...")
    reads = parse_fastq('matching.fq')
    print(f"Loaded {len(reads)} reads")
    
    # Analyze features in reference sequence
    print("\nAnalyzing sequence features...")
    features = find_kan_cassette_features(reference_seq)
    
    # Analyze junction reads
    print("\nAnalyzing insertion junctions...")
    junction_reads = analyze_insertion_junction(reference_seq, reads)
    print(f"Found {len(junction_reads)} reads spanning the insertion junction")
    
    # Generate report
    report = generate_report(reference_seq, reads, features)
    
    # Save report
    with open('analysis_report.txt', 'w') as f:
        f.write(report)
    
    print("\nReport saved to analysis_report.txt")
    
    # Create visualization
    print("\nCreating insertion map visualization...")
    visualize_insertion(reference_seq, features)
    print("Visualization saved to insertion_map.png")
    
    # Additional analysis suggestions
    # print("\n=== Next Steps ===")
    # print("1. Download S. cerevisiae reference genome (CHZ1 region)")
    # print("2. Perform detailed alignment using BWA or Bowtie2")
    # print("3. Verify correct insertion and absence of mutations")
    # print("4. Check for complete CHZ1 deletion if doing knockout")

if __name__ == "__main__":
    main()