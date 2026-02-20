#!/usr/bin/env python3
import pandas as pd
import numpy as np

def analyze_file(filename):
    """Read file and calculate statistics."""
    data = pd.read_csv(filename, header=None, names=['value'])
    values = data['value']

    # Calculate percentiles
    q25, q50, q75 = values.quantile([0.25, 0.5, 0.75])

    # Top 1% vs lower 99%
    threshold = values.quantile(0.99)
    top_1_pct = values[values >= threshold].mean()
    lower_99_pct = values[values < threshold].mean()

    return {
        'Mean': values.mean(),
        'Q25': q25,
        'Median': q50,
        'Q75': q75,
        'Max': values.max(),
        'Top 1% Mean': top_1_pct,
        'Lower 99% Mean': lower_99_pct
    }

# Analyze all files
files = ['as.txt', 'ns.txt', 'intv.txt', 'nInactive.txt', 'ha_sizes.txt', 'ha_nodes.txt', 'scores.txt', 'graphSizeFull.txt', 
         'graphSizeEssential.txt', 'graphSize2BitSeq.txt', 'graphNodes.txt', 
         'graphEdges.txt',# 'graphNodeLen.txt',
         'bfsTime.txt', 'bfsNodes.txt', 'bfsEdges.txt', 'bfsSize2BitSeq.txt',
         ]
results = {}

for f in files:
    results[f.replace('.txt', '')] = analyze_file(f)

# Create and display table
df = pd.DataFrame(results).T
print("\nStatistics Summary")
print("=" * 90)
print(df.to_string(float_format='%.2f'))
print("=" * 90)
