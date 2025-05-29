# Yeast CHZ1 Kanamycin Insertion Analysis

This project provides visualization and analysis tools for studying the insertion of a Kanamycin resistance cassette (KanMX) into the CHZ1 gene in Saccharomyces cerevisiae.

## Overview

The CHZ1 gene (YER030W) encodes a histone chaperone for Htz1p/H2A-H2B dimer in yeast. This project focuses on analyzing and visualizing the targeted disruption of this gene using a Kanamycin resistance cassette insertion.

## Features

- Generation of publication-quality visualizations:
  - Schematic diagrams of the CHZ1 locus before and after Kanamycin insertion
  - Detailed junction analysis with restriction sites
  - Sequencing read coverage visualization
  - Restriction map of the integrated cassette
- Comprehensive data analysis and reporting
- Generation of paper-ready figures and summaries

## Requirements

- Python 3.x
- Required Python packages:
  - matplotlib
  - numpy
  - BioPython

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/Yeast-CHZ1-Kanmx-Analysis.git
cd Yeast-CHZ1-Kanmx-Analysis
```

2. Install required dependencies:
```bash
pip install matplotlib numpy biopython
```

## Usage

Run the main script to generate all visualizations and data summaries:

```bash
python insert_visualization.py
```

This will generate the following output files:
- `chz1_insertion_visualization.png/pdf`: Main visualization of the insertion
- `restriction_map.png/pdf`: Restriction map of the integrated cassette
- `paper_data_summary.txt`: Comprehensive data summary for publication

## Output Files

### Visualizations
- **chz1_insertion_visualization.png/pdf**
  - Panel A: CHZ1 locus before and after Kanamycin insertion
  - Panel B: Detailed junction view with restriction sites
  - Panel C: Sequencing read coverage analysis

- **restriction_map.png/pdf**
  - Detailed restriction map showing positions of:
    - SalI sites
    - BamHI sites
    - PacI sites
    - AscI sites

### Data Summary
- **paper_data_summary.txt**
  - Title suggestion
  - Key results
  - Methods summary
  - Figure legends
  - Recommended next experiments

## Analysis Results

The analysis shows:
- Successful integration confirmed by 32 junction-spanning reads
- Clean junction without indels
- Total of 265 sequencing reads analyzed
- 178 upstream reads
- 32 junction-spanning reads
- 13 cassette-only reads

## Next Steps

Recommended follow-up experiments:
1. PCR verification with external primers
2. Southern blot to confirm single integration
3. Phenotypic analysis (MMS and benomyl sensitivity)
4. Complementation test with wild-type CHZ1

## License

[Add your license information here]

## Contact

[Add your contact information here] 