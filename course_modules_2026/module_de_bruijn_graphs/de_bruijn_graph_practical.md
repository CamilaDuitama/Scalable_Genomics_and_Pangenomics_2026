## Prerequisites

Download the data from the de bruijn graph data folder. 

Software used this session:
 - python

## Test out some graphs! 

de bruijn graphs are an intuitive graph structure. This module is designed to get everyone comfortable with how graphs work 
and the relationships between graphs and 

### building the first de bruijn graph

copy and paste the first set of reads into the "reads" box on this website: https://debruijn.mortality.watch/

Ask yourself, how does the k length change the graph? 
How do the contigs change with different K? 
Why do you think that is? 

### building the more de bruijn graphs 

Use the rest of the given sets to explore how read length impacts the graph structure. 
Again, for each read length, change the k sliding bar to see how k impacts the graph and the contig number. 
For larger reads, how does the graph change? 

## Make your own reads 

Use the python script to simulate reads with different parameters. 

'''
python simulate_reads.py \
    --length 100 \
    --coverage 2 \
    --read-length 5 \
    --error-rate 0.01 \
    --heterozygosity 0.001 \
    --seed 42
'''

Copy and paste these reads into the de bruijn graph interface. 


