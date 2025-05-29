#!/usr/bin/env python3
"""
Create a flowchart of the bioinformatics pipeline for CHZ1::KanMX analysis
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import matplotlib.lines as mlines

def create_pipeline_flowchart():
    """Create a detailed flowchart of the bioinformatics analysis pipeline"""
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 14))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 20)
    ax.axis('off')
    
    # Color scheme
    colors = {
        'input': '#E3F2FD',      # Light blue
        'process': '#FFF3E0',    # Light orange
        'analysis': '#E8F5E9',   # Light green
        'output': '#FCE4EC',     # Light pink
        'decision': '#FFFDE7'    # Light yellow
    }
    
    # Title
    ax.text(5, 19.5, 'Bioinformatics Pipeline for CHZ1::KanMX Verification', 
            ha='center', va='center', fontsize=16, fontweight='bold')
    
    # Helper function to create boxes
    def create_box(x, y, width, height, text, color, style='round'):
        if style == 'round':
            box = FancyBboxPatch((x-width/2, y-height/2), width, height,
                                boxstyle="round,pad=0.1", 
                                facecolor=color, edgecolor='black', linewidth=1.5)
        else:  # Diamond for decision
            box = FancyBboxPatch((x-width/2, y-height/2), width, height,
                                boxstyle="round,pad=0.1", 
                                facecolor=color, edgecolor='black', linewidth=1.5,
                                transform=ax.transData)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=10, 
                fontweight='bold', wrap=True)
        return box
    
    # Helper function to create arrows
    def create_arrow(x1, y1, x2, y2, label=''):
        arrow = FancyArrowPatch((x1, y1), (x2, y2),
                               connectionstyle="arc3,rad=0", 
                               arrowstyle='->', lw=1.5, color='black')
        ax.add_patch(arrow)
        if label:
            mid_x, mid_y = (x1+x2)/2, (y1+y2)/2
            ax.text(mid_x+0.3, mid_y, label, fontsize=8, style='italic')
    
    # Input data
    y_pos = 18
    create_box(2.5, y_pos, 4, 1, 'Raw FASTQ File\n(265 reads)', colors['input'])
    create_box(7.5, y_pos, 4, 1, 'Reference Sequence\n(paste.txt)', colors['input'])
    
    # Parse and QC
    y_pos = 16
    create_box(5, y_pos, 5, 1, 'Data Parsing & Quality Control', colors['process'])
    create_arrow(2.5, 17.5, 5, 16.5)
    create_arrow(7.5, 17.5, 5, 16.5)
    
    # Sequence extraction
    y_pos = 14.5
    create_box(2.5, y_pos, 3.5, 0.8, 'Extract Read\nSequences', colors['process'])
    create_box(7.5, y_pos, 3.5, 0.8, 'Parse Reference\nSequence', colors['process'])
    create_arrow(5, 15.5, 2.5, 15)
    create_arrow(5, 15.5, 7.5, 15)
    
    # Feature identification
    y_pos = 13
    create_box(5, y_pos, 5, 1, 'Feature Identification', colors['analysis'])
    create_arrow(2.5, 14.1, 5, 13.5)
    create_arrow(7.5, 14.1, 5, 13.5)
    
    # Parallel analyses
    y_pos = 11.5
    create_box(2, y_pos, 3, 0.8, 'Restriction Site\nMapping', colors['analysis'])
    create_box(5, y_pos, 3, 0.8, 'Read\nClassification', colors['analysis'])
    create_box(8, y_pos, 3, 0.8, 'Junction\nDetection', colors['analysis'])
    
    create_arrow(5, 12.5, 2, 11.9)
    create_arrow(5, 12.5, 5, 11.9)
    create_arrow(5, 12.5, 8, 11.9)
    
    # Sub-processes
    y_pos = 10
    ax.text(2, y_pos, '• BamHI: 16 sites\n• SalI: 14 sites\n• PacI: 21 sites\n• AscI: 28 sites',
            ha='center', va='top', fontsize=8)
    ax.text(5, y_pos, '• Upstream: 178\n• Cassette: 13\n• Junction: 32\n• Other: 42',
            ha='center', va='top', fontsize=8)
    ax.text(8, y_pos, '• GTCGAC|GGATCC\n• 32 confirmed\n  junctions',
            ha='center', va='top', fontsize=8)
    
    # Integration point
    y_pos = 8
    create_box(5, y_pos, 5, 1, 'Data Integration & Validation', colors['analysis'])
    create_arrow(2, 10, 5, 8.5)
    create_arrow(5, 10, 5, 8.5)
    create_arrow(8, 10, 5, 8.5)
    
    # Decision point
    y_pos = 6.5
    create_box(5, y_pos, 4, 1, 'Junction Reads\nPresent?', colors['decision'])
    create_arrow(5, 7.5, 5, 7)
    
    # Outcomes
    y_pos = 5
    create_box(3, y_pos, 3, 0.8, 'Failed\nIntegration', colors['output'])
    create_box(7, y_pos, 3, 0.8, 'Successful\nIntegration', colors['output'])
    
    create_arrow(5, 6, 3, 5.4, 'No')
    create_arrow(5, 6, 7, 5.4, 'Yes (32 reads)')
    
    # Coverage analysis
    y_pos = 3.5
    create_box(7, y_pos, 3.5, 0.8, 'Coverage\nAnalysis', colors['analysis'])
    create_arrow(7, 4.6, 7, 3.9)
    
    # Final outputs
    y_pos = 2
    create_box(2, y_pos, 3, 0.8, 'Troubleshooting\nRecommendations', colors['output'])
    create_box(5, y_pos, 3, 0.8, 'Visualization\nGeneration', colors['output'])
    create_box(8, y_pos, 3, 0.8, 'Statistical\nReport', colors['output'])
    
    create_arrow(3, 4.6, 2, 2.4)
    create_arrow(7, 3.1, 5, 2.4)
    create_arrow(7, 3.1, 8, 2.4)
    
    # Final integration
    y_pos = 0.5
    create_box(5, y_pos, 6, 0.8, 'Comprehensive Analysis Report', colors['output'])
    create_arrow(2, 1.6, 5, 0.9)
    create_arrow(5, 1.6, 5, 0.9)
    create_arrow(8, 1.6, 5, 0.9)
    
    # Add legend
    legend_x = 0.2
    legend_y = 4
    ax.text(legend_x, legend_y+0.5, 'Legend:', fontsize=10, fontweight='bold')
    
    for i, (label, color) in enumerate(colors.items()):
        rect = Rectangle((legend_x, legend_y-i*0.4-0.3), 0.3, 0.3, 
                        facecolor=color, edgecolor='black')
        ax.add_patch(rect)
        ax.text(legend_x+0.4, legend_y-i*0.4-0.15, label.capitalize(), 
                fontsize=8, va='center')
    
    # Add annotations
    ax.text(9.5, 12, 'Key Finding:\n32 junction reads\nconfirm integration', 
            fontsize=9, ha='center', 
            bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.3))
    
    ax.text(9.5, 9, 'Total:\n79 restriction\nsites identified', 
            fontsize=9, ha='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.3))

    plt.tight_layout()
    plt.savefig('bioinformatics_pipeline_flowchart.png', dpi=300, bbox_inches='tight')
    plt.savefig('bioinformatics_pipeline_flowchart.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    
    print("Flowchart saved as 'bioinformatics_pipeline_flowchart.png' and '.pdf'")

if __name__ == "__main__":
    create_pipeline_flowchart()