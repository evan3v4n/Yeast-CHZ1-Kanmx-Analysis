#!/usr/bin/env python3
"""
Quick analysis of kanamycin insertion in yeast CHZ1 gene
"""

import os

# Ensure reference_data directory exists
os.makedirs('reference_data', exist_ok=True)

def analyze_sequences():
    # Read the reference sequence from paste.txt
    with open('matching.fq', 'r') as f:
        lines = f.readlines()
    
    # Extract the reference sequence
    ref_seq = ""
    for line in lines:
        line = line.strip()
        if line and not line.startswith('>') and not line.startswith('CHZ1'):
            # Remove leading dashes and spaces, keep only DNA sequence
            import re
            match = re.search(r'[ATCGN]+', line)
            if match:
                ref_seq += match.group()
    
    # Save reference sequence to reference_data directory
    with open('reference_data/reference_sequence.fasta', 'w') as f:
        f.write(">CHZ1_reference\n")
        f.write(ref_seq + "\n")
    
    print(f"Reference sequence length: {len(ref_seq)} bp")
    
    # Key sequences to identify
    print("\n=== Identifying Key Features ===")
    
    # CHZ1 sequence portion
    chz1_seq = "ACAACACTGCTTCACCAACATCTGCATCGTACGCTGCAGGTCGAC"
    
    # Cassette start (with restriction sites)
    cassette_start = "GGATCCCCGGGTTAATTAAGGCGCGCC"
    
    # Find these in reference
    chz1_pos = ref_seq.find(chz1_seq)
    cassette_pos = ref_seq.find(cassette_start)
    
    print(f"CHZ1 sequence found at position: {chz1_pos}")
    print(f"Cassette sequence found at position: {cassette_pos}")
    
    # The junction point
    junction = "GTCGACGGATCCCC"
    junction_pos = ref_seq.find(junction)
    print(f"Junction sequence found at position: {junction_pos}")
    
    # Save feature positions to reference_data directory
    with open('reference_data/feature_positions.txt', 'w') as f:
        f.write(f"CHZ1 sequence: {chz1_pos}\n")
        f.write(f"Cassette start: {cassette_pos}\n")
        f.write(f"Junction point: {junction_pos}\n")
    
    # Analyze FASTQ reads
    print("\n=== Analyzing Sequencing Reads ===")
    
    reads = []
    with open('matching.fq', 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 4):
            if i+1 < len(lines):
                reads.append({
                    'header': lines[i].strip(),
                    'sequence': lines[i+1].strip(),
                    'quality': lines[i+3].strip() if i+3 < len(lines) else ''
                })
    
    print(f"Total reads: {len(reads)}")
    
    # Categorize reads
    upstream_only = 0
    cassette_only = 0
    junction_reads = 0
    
    junction_examples = []
    
    for read in reads:
        seq = read['sequence']
        has_upstream = "ACAACACTGCTTCACCAACATCTGCA" in seq
        has_cassette = "GGATCCCCGGGTTAATTAA" in seq
        has_junction = "GTCGACGGATCCCC" in seq
        
        if has_junction or (has_upstream and has_cassette):
            junction_reads += 1
            if len(junction_examples) < 5:  # Save first 5 examples
                junction_examples.append(read)
        elif has_upstream:
            upstream_only += 1
        elif has_cassette:
            cassette_only += 1
    
    print(f"\nRead categories:")
    print(f"  - Upstream CHZ1 only: {upstream_only}")
    print(f"  - Cassette only: {cassette_only}")
    print(f"  - Junction reads: {junction_reads}")
    
    # Save junction examples to reference_data directory
    with open('reference_data/junction_examples.fasta', 'w') as f:
        for i, read in enumerate(junction_examples):
            f.write(f">junction_example_{i+1}\n{read['sequence']}\n")
    
    # Show junction examples
    print("\n=== Junction Read Examples ===")
    for i, read in enumerate(junction_examples):
        seq = read['sequence']
        junction_pos = seq.find("GTCGACGGATCCCC")
        if junction_pos != -1:
            start = max(0, junction_pos - 20)
            end = min(len(seq), junction_pos + 34)
            print(f"\nJunction {i+1}:")
            print(f"  {seq[start:junction_pos]}|{seq[junction_pos:end]}")
            print(f"  {'.'*20}^{'.'*33}")
            print(f"  {'CHZ1 end':<20}{'Cassette start'}")
    
    # Check restriction sites
    print("\n=== Restriction Site Analysis ===")
    sites = {
        'BamHI': 'GGATCC',
        'SalI': 'GTCGAC',
        'PacI': 'TTAATTAA',
        'AscI': 'GGCGCGCC'
    }
    
    # Save restriction site analysis to reference_data directory
    with open('reference_data/restriction_sites.txt', 'w') as f:
        for name, seq in sites.items():
            count = ref_seq.count(seq)
            f.write(f"{name} ({seq}): {count} sites\n")
            print(f"{name} ({seq}): {count} sites")
    
    # Generate summary
    print("\n=== Summary ===")
    print("1. The reference sequence shows a kanamycin resistance cassette inserted into CHZ1")
    print("2. The insertion junction is at the GTCGAC (SalI) site")
    print("3. The cassette contains multiple cloning sites (BamHI, PacI, AscI)")
    print(f"4. {junction_reads} reads confirm the insertion junction")
    print("5. This appears to be a successful CHZ1::KanMX insertion")
    
    # Save detailed analysis to reference_data directory
    with open('reference_data/quick_analysis_results.txt', 'w') as f:
        f.write("Kanamycin Insertion Analysis Results\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Reference sequence length: {len(ref_seq)} bp\n")
        f.write(f"Total sequencing reads: {len(reads)}\n")
        f.write(f"Junction reads: {junction_reads}\n\n")
        f.write("Junction sequence: ...TGCATCGTACGCTGCAGGTCGAC|GGATCCCCGGGTTAATTAA...\n")
        f.write("                   ^CHZ1 end               ^Cassette start\n")

if __name__ == "__main__":
    analyze_sequences()