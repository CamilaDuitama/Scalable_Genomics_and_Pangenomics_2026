# Practical: Exploring de Bruijn Graphs

## Learning objectives

By the end of this practical, you should be able to:

* Explain how a de Bruijn graph is constructed from sequencing reads.
* Describe how the choice of *k*-mer length affects graph structure.
* Explain how read length, sequencing coverage, sequencing errors, and heterozygosity influence graph complexity.
* Explore the relationship between graph structure, contigs, and genome assembly.
* Use Python to simulate sequencing reads and investigate how changing simulation parameters affects assembly.

## Prerequisites

### Data

Download the example datasets from the **de Bruijn graph data** folder.

### Software and resources

* Python
* The supplied `simulate_reads.py` script
* The [de Bruijn graph visualisation tool](https://debruijn.mortality.watch/)

No previous experience with graph theory is required.

---

## Part 1: Exploring de Bruijn graphs

A de Bruijn graph is a graph representation of sequence data. In genome assembly, sequencing reads are broken into overlapping *k*-mers, which are then used to construct a graph.

The choice of *k*, the amount of sequencing data, and the characteristics of the genome can all influence the graph. These changes can affect how easily the original genome can be reconstructed.

In this practical, you will explore these relationships using a graph visualisation tool.

### 1.1 Build your first de Bruijn graph

1. Open the [de Bruijn graph visualisation tool](https://debruijn.mortality.watch/).
2. Copy the first set of reads from the supplied dataset.
3. Paste the reads into the `reads` box.
4. Generate the graph using the default settings.
5. Examine the graph before changing any parameters.

**Questions**

* What do the nodes and edges represent in this graph?
* How are overlapping sequence fragments represented?
* Does the graph contain any branches? If so, what might they represent?
* How might you use this graph to reconstruct the original sequence?

Write down your initial observations before moving on.

### 1.2 Investigate the effect of *k*-mer length

The parameter *k* determines the length of the sequence fragments used to construct the graph.

Using the same set of reads, change the *k*-mer length several times. Keep the reads unchanged so that you can isolate the effect of *k*.

For each value of *k*:

1. Predict how the graph will change.
2. Generate the graph.
3. Record the number of nodes, branches, and contigs, where these are reported by the tool.
4. Compare the result with your prediction.

**Questions**

* What happens to the graph as *k* increases?
* What happens as *k* decreases?
* How does the number of contigs change?
* Why might a small *k* make it difficult to distinguish between repeated regions of a genome?
* Why might a large *k* make it more difficult to assemble a genome when sequencing coverage is low?

**Think about it:** A smaller *k* can connect more sequence fragments, but may introduce ambiguity. A larger *k* can provide more sequence context, but requires sufficient read length and coverage to observe the necessary *k*-mers.

### 1.3 Investigate the effect of read length

Now repeat the investigation using the other supplied read datasets, which contain different read lengths.

For each dataset:

1. Load the reads into the visualisation tool.
2. Test several values of *k* that are shorter than the read length.
3. Record the graph structure and number of contigs.
4. Compare the results across datasets.

Keep *k* constant when comparing different read lengths, then vary *k* to investigate the interaction between read length and *k*-mer length.

**Questions**

* How does read length affect the graph when *k* is held constant?
* Which datasets produce the simplest or most complex graphs?
* Do longer reads always produce fewer contigs? Explain your observations.
* Why must *k* be shorter than the read length for a read to contain a *k*-mer?
* How might longer sequencing reads help resolve repeated regions?

### Observation table

Use this table to organise your results.

| Dataset   | Read length | *k* | Graph observations | Number of contigs |
| --------- | ----------: | --: | ------------------ | ----------------: |
| Dataset 1 |             |     |                    |                   |
| Dataset 1 |             |     |                    |                   |
| Dataset 2 |             |     |                    |                   |
| Dataset 3 |             |     |                    |                   |


---

## Part 2: Simulating your own sequencing reads

The supplied Python script allows you to simulate sequencing reads while changing the properties of the genome and the sequencing process.

This gives you more control over the data than the example datasets.

### 2.1 Generate a baseline dataset

Run the following command in your terminal from the directory containing `simulate_reads.py`:

```bash
python simulate_reads.py \
    --length 100 \
    --coverage 2 \
    --read-length 5 \
    --error-rate 0.01 \
    --heterozygosity 0.001 \
    --seed 42
```

The parameters control different aspects of the simulation:

| Parameter          | Meaning                                                                         |
| ------------------ | ------------------------------------------------------------------------------- |
| `--length`         | Length of the simulated genome                                                  |
| `--coverage`       | Sequencing coverage                                                             |
| `--read-length`    | Length of each simulated read                                                   |
| `--error-rate`     | Rate of sequencing errors                                                       |
| `--heterozygosity` | Rate of differences between the simulated haplotypes, if modelled by the script |
| `--seed`           | Random seed for reproducibility                                                 |

Consult the script's help output or source code if you want to confirm exactly how each parameter is implemented.

**Important:** These are VERY small example parameters for exploration, not realistic values for a genome assembly.

Once the reads have been generated:

1. Inspect the output from the script.
2. Copy the reads from the terminal into the visualisation tool.
3. Generate a graph.
4. Explore several *k*-mer lengths.
5. Record your observations.

This is the baseline dataset, now change one parameter at a time in the following experiments so that you can identify its effect.

### 2.2 Experiment A: Change the read length

Generate datasets using shorter and longer reads while keeping the other parameters fixed.

For example, compare read lengths of 5, 10, and 20 bases, provided these values are supported by the script and are appropriate for your genome length.

**Questions**

* How does read length affect the number of distinct *k*-mers?
* How does it affect graph connectivity and the number of contigs?
* Which read lengths work best for the *k*-mer values you tested?
* Can longer reads provide more information about the order of sequence fragments?

### 2.3 Experiment B: Change the genome length

Increase and decrease `--length` while keeping the other parameters fixed.

**Questions**

* How does genome length affect the number of reads and the size of the graph?
* Does a longer genome necessarily produce more contigs?
* How would you distinguish an effect caused by genome length from one caused by insufficient sequencing coverage?

### 2.4 Experiment C: Change sequencing coverage

Coverage describes how much sequencing data is available relative to the genome length. Higher coverage means that sequence positions are sampled more frequently on average.

Generate datasets with low, intermediate, and high coverage. For example, try coverage values of 1, 2, 5, and 10 if supported by the script.

**Questions**

* How does increasing coverage change the graph?
* What happens to regions that are poorly represented in the reads?
* Does higher coverage always improve assembly?
* Why might some *k*-mers be absent from a low-coverage dataset?
* How might coverage interact with *k*-mer length?

### 2.5 Experiment D: Change the sequencing error rate

Compare datasets generated with no sequencing errors and with increasing error rates.

For example, try error rates of 0, 0.01, and 0.05.

**Questions**

* What happens to the graph when the error rate is zero?
* How do sequencing errors affect the number of distinct *k*-mers?
* Do errors create additional branches or disconnected parts of the graph?
* Why can sequencing errors make genome assembly more difficult?
* Would increasing coverage help compensate for sequencing errors? Why or why not?

### 2.6 Experiment E: Change heterozygosity

Heterozygosity describes genetic differences between homologous chromosomes or haplotypes in a diploid genome.

Generate datasets with different heterozygosity values, keeping the other parameters fixed.

For example, compare low and high heterozygosity using values supported by the simulation script.

**Questions**

* How do sequence differences between haplotypes affect the graph?
* Can heterozygous sites introduce branches into the graph?
* How might you distinguish a branch caused by heterozygosity from one caused by sequencing errors?
* What might happen when heterozygosity is high but sequencing coverage is low?

**Further investigation:** Combine high heterozygosity with low coverage. Compare the resulting graph with one generated using low heterozygosity and high coverage.

---

## Part 3: Putting it all together

Choose two parameters to investigate together. For example:

* Read length and *k*-mer length
* Coverage and *k*-mer length
* Sequencing error rate and coverage
* Heterozygosity and coverage
* Genome length and coverage

Before running the simulation, write down a prediction. Then generate the data, build the graph, and record the result.

Use the following table to document your experiment.

| Experiment | Parameters changed | Prediction | Observation | Interpretation |
| ---------- | ------------------ | ---------- | ----------- | -------------- |
| 1          |                    |            |             |                |
| 2          |                    |            |             |                |
| 3          |                    |            |             |                |

### Final discussion

Use your observations to answer the following questions:

1. Which parameter had the largest effect on graph structure in your experiments?
2. Which conditions produced the most fragmented assembly?
3. What is the trade-off between small and large *k*-mer lengths?
4. How do read length and coverage influence your choice of *k*?
5. How can sequencing errors and heterozygosity produce similar-looking graph structures?
6. What additional information would you need to decide whether a branch represents a repeat, a sequencing error, or a genuine genetic difference?
7. If you were planning a real sequencing experiment, what factors would you consider when choosing read length and sequencing depth?

---

## Summary

De Bruijn graphs provide a useful way to understand how sequencing data are assembled into longer sequences. Their structure depends on both the genome and the sequencing process.

The key lesson is that **there is no single best *k*-mer length or sequencing configuration for every dataset**. Assembly depends on the interaction between read length, coverage, genome complexity, sequencing errors, and genetic variation.

Your task is to use the simulations to develop an evidence-based explanation of these relationships.




