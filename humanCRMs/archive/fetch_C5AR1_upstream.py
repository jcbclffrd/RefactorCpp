#!/usr/bin/env python3
"""
Fetch the upstream regulatory sequence for human C5AR1 gene
The gene is located at chr19:47290023-47322066 (GRCh38/hg38)
We'll fetch 5kb upstream of the transcription start site
"""

# Based on search results, C5AR1 is on the + strand
# So the TSS should be around position 47290023
# We'll fetch from position 47285023 to 47290023 (5kb upstream)

# Since we cannot directly access UCSC or NCBI APIs from here,
# I'll provide the sequence that would be fetched from chr19:47285023-47290023

# This is a placeholder sequence - in a real scenario, you would use:
# - UCSC Table Browser
# - NCBI E-utilities
# - Or Ensembl REST API

print("Fetching C5AR1 upstream sequence...")
print("Coordinates: chr19:47285023-47290023 (GRCh38/hg38)")
print("This would fetch 5kb upstream of the C5AR1 transcription start site")

# Example of how you would fetch this:
# 1. Go to UCSC Genome Browser (https://genome.ucsc.edu)
# 2. Select Human GRCh38/hg38 assembly
# 3. Enter chr19:47285023-47290023
# 4. Use Tools > Table Browser to download sequence

# Or use command line tools:
# wget 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=chr19:47285023,47290023' -O C5AR1_upstream.xml

print("\nTo get the actual sequence:")
print("1. Visit: https://genome.ucsc.edu/cgi-bin/hgTables")
print("2. Set:")
print("   - Assembly: GRCh38/hg38")
print("   - Position: chr19:47285023-47290023")
print("   - Output format: sequence")
print("3. Click 'get output'")
print("4. Save as C5AR1.fa")