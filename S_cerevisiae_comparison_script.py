#!/usr/bin/env python3
"""
S. cerevisiae Genome Comparison Script
Compares sequencing data against S. cerevisiae reference genome
"""

import subprocess
import os
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import matplotlib.pyplot as plt
import numpy as np

# Set email for NCBI queries
Entrez.email = "your_email@example.com"

def download_chz1_reference():
    """Download CHZ1 gene sequence from SGD"""
    print("Downloading CHZ1 reference sequence...")
    
    # SGD REST API endpoint for CHZ1
    url = "https://www.yeastgenome.org/webservice/locus/YER030W/sequence"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            # Extract genomic sequence
            genomic_seq = data.get('genomic_sequence', '')
            
            # Save to file
            with open('chz1_reference.fasta', 'w') as f:
                f.write(">YER030W CHZ1 S. cerevisiae\n")
                f.write(genomic_seq + "\n")
            
            print(f"Downloaded CHZ1 sequence: {len(genomic_seq)} bp")
            return genomic_seq
        else:
            print("Failed to download from SGD, using NCBI...")
            return download_from_ncbi()
    except Exception as e:
        print(f"Error downloading from SGD: {e}")
        return download_from_ncbi()

def download_from_ncbi():
    """Alternative: Download from NCBI"""
    handle = Entrez.efetch(db="nucleotide", id="NC_001137.3", 
                          rettype="fasta", retmode="text",
                          seq_start=98734, seq_stop=99477)  # CHZ1 coordinates
    record = SeqIO.read(handle, "fasta")
    handle.close()
    
    with open('chz1_reference.fasta', 'w') as f:
        SeqIO.write(record, f, "fasta")
    
    return str(record.seq)

def prepare_blast_database():
    """Create BLAST database from reference sequence"""
    print("\nCreating BLAST database...")
    
    # Create a combined reference with CHZ1 and known cassette sequences
    with open('combined_reference.fasta', 'w') as f:
        # Add CHZ1 reference
        f.write(">CHZ1_genomic\n")
        if os.path.exists('chz1_reference.fasta'):
            with open('chz1_reference.fasta', 'r') as ref:
                lines = ref.readlines()
                f.write(''.join(lines[1:]))
        
        # Add common kanamycin cassette sequences
        f.write("\n>KanMX4_cassette\n")
        f.write("ATGGCTAAAATGAGAATATCACCGGAATTGAAGGCGATATAGTCCTGACGCAGATCGCCG\n")
        
        f.write("\n>TEF_promoter\n")
        f.write("ATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCA\n")
    
    # Create BLAST database
    cmd = ["makeblastdb", "-in", "combined_reference.fasta", 
           "-dbtype", "nucl", "-out", "reference_db"]
    
    try:
        subprocess.run(cmd, check=True)
        print("BLAST database created successfully")
    except subprocess.CalledProcessError:
        print("Note: makeblastdb not found. Install BLAST+ for local alignment")
        print("Continuing with sequence comparison...")

def compare_sequences(query_file, reference_seq):
    """Compare sequences using various methods"""
    print("\n=== Sequence Comparison Analysis ===")
    
    # Extract unique sequences from FASTQ
    unique_seqs = set()
    with open(query_file, 'r') as f:
        lines = f.readlines()
        for i in range(1, len(lines), 4):  # Every 4th line starting from 1
            unique_seqs.add(lines[i].strip())
    
    print(f"Found {len(unique_seqs)} unique sequences")
    
    # Analyze each unique sequence
    results = {
        'perfect_matches': 0,
        'partial_matches': 0,
        'junction_sequences': 0,
        'cassette_only': 0,
        'genomic_only': 0
    }
    
    # Key sequences to look for
    cassette_markers = [
        "GGATCCCCGGGTTAATTAA",  # Cassette start
        "GGCGCGCCAGATCTGTTTA",  # Within cassette
        "GTCGACGGATCCCCGGGTT"   # Junction region
    ]
    
    genomic_markers = [
        "ACAACACTGCTTCACCAACATCTGCA",  # CHZ1 sequence
        "AAATTGGCGTATTATTACAACACTGC"   # CHZ1 upstream
    ]
    
    junction_sequences = []
    
    for seq in unique_seqs:
        has_cassette = any(marker in seq for marker in cassette_markers)
        has_genomic = any(marker in seq for marker in genomic_markers)
        
        if has_cassette and has_genomic:
            results['junction_sequences'] += 1
            junction_sequences.append(seq)
        elif has_cassette:
            results['cassette_only'] += 1
        elif has_genomic:
            results['genomic_only'] += 1
        
        # Check for perfect match to reference
        if seq in reference_seq:
            results['perfect_matches'] += 1
        elif any(seq[i:i+50] in reference_seq for i in range(0, len(seq)-50, 10)):
            results['partial_matches'] += 1
    
    # Save junction sequences for further analysis
    with open('junction_sequences.fasta', 'w') as f:
        for i, seq in enumerate(junction_sequences[:10]):  # Save first 10
            f.write(f">junction_{i+1}\n{seq}\n")
    
    return results

def identify_insertion_site(reference_seq, junction_sequences):
    """Identify exact insertion site"""
    print("\n=== Insertion Site Analysis ===")
    
    # Expected junction based on homologous recombination
    expected_upstream = "ATCTGCATCGTACGCTGCAGGTCGAC"
    expected_downstream = "GGATCCCCGGGTTAATTAA"
    
    insertion_sites = []
    
    for seq in junction_sequences:
        # Find the junction point
        if expected_upstream in seq and expected_downstream in seq:
            pos1 = seq.find(expected_upstream)
            pos2 = seq.find(expected_downstream)
            
            if pos2 > pos1:
                junction_point = pos1 + len(expected_upstream)
                insertion_sites.append({
                    'sequence': seq,
                    'junction_position': junction_point,
                    'upstream_context': seq[max(0, pos1-20):pos1],
                    'downstream_context': seq[pos2+len(expected_downstream):pos2+len(expected_downstream)+20]
                })
    
    return insertion_sites

def generate_visualizations(results, insertion_sites):
    """Generate visualizations for the analysis results"""
    # Create a figure with multiple subplots
    fig = plt.figure(figsize=(15, 10))
    
    # 1. Pie chart for sequence classification
    plt.subplot(2, 2, 1)
    labels = ['Genomic Only', 'Cassette Only', 'Junction Sequences']
    sizes = [results['genomic_only'], results['cassette_only'], results['junction_sequences']]
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.title('Sequence Classification Distribution')
    
    # 2. Bar chart for match types
    plt.subplot(2, 2, 2)
    match_types = ['Perfect Matches', 'Partial Matches']
    match_counts = [results['perfect_matches'], results['partial_matches']]
    plt.bar(match_types, match_counts)
    plt.title('Sequence Match Types')
    plt.ylabel('Count')
    
    # 3. Junction sequence length distribution
    if insertion_sites:
        plt.subplot(2, 2, 3)
        lengths = [len(site['sequence']) for site in insertion_sites]
        plt.hist(lengths, bins=20)
        plt.title('Junction Sequence Length Distribution')
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Count')
    
    # 4. Junction position distribution
    if insertion_sites:
        plt.subplot(2, 2, 4)
        positions = [site['junction_position'] for site in insertion_sites]
        plt.hist(positions, bins=20)
        plt.title('Junction Position Distribution')
        plt.xlabel('Position in Read')
        plt.ylabel('Count')
    
    plt.tight_layout()
    plt.savefig('analysis_visualizations.png')
    plt.close()

def generate_comparison_report(results, insertion_sites):
    """Generate detailed comparison report"""
    report = []
    report.append("=== S. cerevisiae Genome Comparison Report ===\n")
    
    # Generate visualizations
    generate_visualizations(results, insertion_sites)
    
    report.append("Sequence Classification:")
    report.append(f"  - Genomic only: {results['genomic_only']}")
    report.append(f"  - Cassette only: {results['cassette_only']}")
    report.append(f"  - Junction sequences: {results['junction_sequences']}")
    report.append(f"  - Perfect matches to reference: {results['perfect_matches']}")
    report.append(f"  - Partial matches: {results['partial_matches']}")
    
    report.append("\nVisualizations have been saved to 'analysis_visualizations.png'")
    
    report.append("\n=== Insertion Site Details ===")
    if insertion_sites:
        site = insertion_sites[0]  # Use first identified site
        report.append(f"Junction position in read: {site['junction_position']}")
        report.append(f"Upstream context: ...{site['upstream_context']}")
        report.append(f"Downstream context: {site['downstream_context']}...")
        report.append("\nInsertion appears to be at the expected site in CHZ1")
    else:
        report.append("No clear junction sequences identified")
    
    report.append("\n=== Recommendations ===")
    report.append("1. Perform PCR verification with primers flanking the insertion")
    report.append("2. Sequence the PCR product to confirm correct insertion")
    report.append("3. Check for any unwanted mutations in flanking regions")
    report.append("4. Verify loss of CHZ1 function if doing knockout")
    
    return "\n".join(report)

def main():
    print("=== S. cerevisiae Genome Comparison Tool ===\n")
    
    # Download reference sequence
    # reference_seq = download_chz1_reference()
    
    # For now, use the provided reference sequence
    with open('paste.txt', 'r') as f:
        text = f.read()
    
    # Extract reference sequence
    import re
    lines = text.strip().split('\n')
    reference_seq = ""
    for line in lines:
        seq_match = re.search(r'[ATCGN]+', line)
        if seq_match and len(seq_match.group()) > 50:
            reference_seq += seq_match.group()
    
    print(f"Reference sequence length: {len(reference_seq)} bp")
    
    # Prepare BLAST database
    prepare_blast_database()
    
    # Compare sequences
    results = compare_sequences('matching.fq', reference_seq)
    
    # Identify insertion sites
    with open('matching.fq', 'r') as f:
        lines = f.readlines()
        junction_seqs = []
        for i in range(1, len(lines), 4):
            seq = lines[i].strip()
            if "GTCGACGGATCCCC" in seq:  # Junction marker
                junction_seqs.append(seq)
    
    insertion_sites = identify_insertion_site(reference_seq, junction_seqs)
    
    # Generate report
    report = generate_comparison_report(results, insertion_sites)
    
    # Save report
    with open('genome_comparison_report.txt', 'w') as f:
        f.write(report)
    
    print("\nReport saved to genome_comparison_report.txt")
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()