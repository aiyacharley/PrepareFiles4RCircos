# PrepareFiles4RCircos
When we want to plot a circos figures to represent several genes or other labels's count and their relation, I write a python script to process data to prepare files for RCircos.

# Input files
1. Germline gene fasta. *human_gl_v*, *human_gl_d*,*human_gl_j*

2. Expression data

3. Recombination data

# Output files
1. human.VDJ.CytoBand.txt

2. human.VDJ.Label.txt

3. human.VDJ.Heatmap.txt

4. human.VDJ.Link.txt

# Result figure
RCircosDemo.pdf

# Scripts
**1. CreateRcircosFiles.py**

**2. RCircos_run.r**

# Run
```
python CreateRcircosFiles.py
Rscript RCircos_run.r
```
