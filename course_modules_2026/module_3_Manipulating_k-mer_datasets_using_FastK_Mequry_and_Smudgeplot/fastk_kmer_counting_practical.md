# Practical: K-mer Counting and Operations Using FastK

## Prerequisites

This practical introduces k-mer counting with **FastK**, inspecting the resulting k-mer database, viewing k-mer histograms and read profiles, and combining k-mer tables with **Logex**.

You should have access to:
- FastK and its companion tools (`Histex`, `Tabex`, `Profex`, and `Logex`)
- The course dataset, including the main trio dataset if supplied
- A terminal in the environment where the FastK tools are installed

## Download the input data

Before counting k-mers, make a directory for this excercise, download the input files and check that they contain sequence data.

```bash
mkdir -p course_data_2026/trio && cd course_data_2026/trio
wget https://kj11.cog.sanger.ac.uk/kmer_course_data/kmer_course_data/{Child,Father,Mother}.Reads.fasta
wget https://kj11.cog.sanger.ac.uk/kmer_course_data/kmer_course_data/Hap{1,2}.fasta
```

these are uncompressed fasta sequences, you can check what is in those files:

```bash
head Child.Reads.fasta
```

## Build k-mer tables with FastK

FastK counts k-mers, generates a sorted table of k-mers & their counts, histogram and for each read/contig/sequence, a k-mer profile if `-p` argument is used. The most essential parameters are `-k` that sets the k-mer length, and `-t` that sets the minimum count to retain in the table. The `-N` option specifies the output k-mer database. Resources are specified by `-T` (the number of threads) and `-M` (for the memory limit). For full list of parameters run `FastK` without any arguments.

Run FastK on the child reads using 21-mers, retaining k-mers with count at least 1 and requesting profiles:

```bash
FastK -v -k21 -t1 -p -T4 Child.Reads.fasta -NChild
```

Repeat for the mother and father reads:

```bash
FastK -v -k21 -t1 -p -T4 Mother.Reads.fasta -NMother
FastK -v -k21 -t1 -p -T4 Father.Reads.fasta -NFather
```

FastK accepts FASTA, FASTQ, SAM, BAM, and CRAM inputs; input files in a single run must be of the same type.

Inspect the output files:

```bash
ls -lh Child*
```

Fair warning, there are also some hidden files in here, you can try `ls -a` to see them

### Looking at the k-mer databases

For that, you use `Histex` too. It can give you number of unique k-mers (default) or k-mer instances (when `-k` is specified).

<details> 
  <summary>Who has the most unique k-mers? Mother, father or Child? </summary>
    ```
    Histex Child | head
    Histex Father | head
    Histex Mother | head
    ```

    Father and Mother have about the same, but father wins with 5,964,873. Child has... a lot fewer.
</details>

To figure out how comes child has so much fewer unique k-mers, let's wait to see the k-mer histograms.

<details> 
  <summary>Who has the most k-mer instnces? Mother, father or Child? </summary>
    ```
    Histex -k Child | head
    Histex -k Father | head
    Histex -k Mother | head
    ```
    Perhaps 
</details>


<details> 
  <summary>Question? </summary>
    Answer
</details>


## Inspect k-mer histograms with Histex

Use `Histex` to print the histogram. In any case, the histogram includes at least two columns: k-mer coverage/count and the number of distinct k-mers observed at that count (frequency). The `-A` option requests tab-delimited ASCII output, while `-k` reports k-mer instances rather than the number of unique k-mers. You can let the output be print, piped or redirected to a file. `-h` enables printing coverages in a certain range. Examples:

```bash

Histex -h1:100 Child | head # yeah, that's the default
Histex -A Child > Child_ascii.hist
Histex Child | tail # unique k-mers
Histex -k Child | tail # k-mer instances; look at the difference in coverages 2, 3, or 10, the pattern should be clear 
```

## k-mer histograms for GenomeScope

The `-G` option produces a GenomeScope-ready histogram, meaning that the histogram provides numbers of unique k-mers for all coverages but the last one, where the agregate of all repretitive k-mers is stored. 

<details> 
  <summary>Create yourself a directory for genomescope estiamtes, get histograms using `Histex` and fit a genomescope model for each of the datasets (Child, Mother, Father). </summary>

    ```bash
    mkdir -p genomescope
    Histex -G Child > genomescope/child.hist.txt
    genomescope2 -i genomescope/child.hist.txt -o gs_child -n child

    Histex -G Father > genomescope/Father.hist.txt
    genomescope2 -i genomescope/Father.hist.txt -o genomescope -n Father

    Histex -G Mother > genomescope/Mother.hist.txt
    genomescope2 -i genomescope/Mother.hist.txt -o genomescope -n Mother
    ```
</details>

**BONUS**

Now, looking at all these GenomeScope profiles, you see that their have genome sizes. 

<details> 
  <summary>How comes Child has so many fewer unique k-mers? Any ideas how to look at it with the tools you were shown already? </summary>

    They have similar genome sizes, so assuming similar repetitiveness they are expected to have a very similar number of unique k-mers in the genome. Meaning, the difference must come from k-mers that are NOT faithfuly representing the genomic template i.e. sequencing errors! The coverage of all datasets is >20x, therefore we can take a shortcut and say, that all k-mers with coverages 1-5x are just error k-mers. Let's compare how many they are in the three datasets. `Histex -h5:6` will do the trick

    ```bash
    Histex -h5:6 Child
    Histex -h5:6 Father 
    Histex -h5:6 Mother 
    ```
    In Child, 38.3% of unique k-mers are error k-mers, in Mother and Father, it's close to 90%! But of course, here if we want to estimate error rate, we should be counting k-mer instances, not unique k-mers.

    ```bash
    Histex -k -h5:6 Child
    Histex -k -h5:6 Father 
    Histex -k -h5:6 Mother 
    ```
    Only 2% of all k-mers are errors in the Child and ~14% in parents! Probability of a k-mer being an genomic k-mer is (1 - e)^k, which is then corresponding to error rate of 0.1% for child (evey 1000th base is an error) and 0.7% for the parents.

</details>

## Inspect k-mer tables with Tabex

`Tabex` can list entries in a k-mer table, check that the table is ordered correctly, or look up a particular k-mer.

```bash
Tabex -C Child
Tabex Child 0-20
```

To suppress entries whose counts are below a threshold:

```bash
Tabex -t10 Child 0-20
```

To look up a particular k-mer, provide a sequence of the correct k-mer length for the table:

```bash
Tabex Child aaaaaaaaaaaaatgacattt
```

The exact result depends on the dataset. A k-mer that is absent from the table is reported as not found. FastK tables contain canonical k-mers, try to search for `aaatgtcattttttttttttt`. 

<details> 
  <summary>Compare which k-mer index is reported in the case of the previous example and this one? Is it the same? </summary>

Ya, they are the reverse complements, meaning they are considered the same sequence if we disregard strand information, which is the default behaviour of FastK.

</details>

## Inspect read profiles with Profex

When FastK is run with `-p`, it creates a profile for each sequence. The profile gives the count of each successive k-mer along that sequence. `Profex` can display selected reads; read numbering starts at 1. The `-z` option compresses consecutive runs of equal profile values, and `-A` requests tab-delimited output.

```bash
Profex Child 1 | head # hard to look at
Profex -z Child 1 | head  # easier to see what's going on
```

<details> 
  <summary>What's going on 21 bases 8-28? What about 38 bases 101-138?</summary>

That will be likely a heterozygous SNP on that location and this read contains only one of the two variants, meaning only 1/2 of coverage of homozygous sites. The 38 bases of "haploid" coverage likely indicates two SNPs 17 bases from each other.

</details>

<details> 
  <summary>Look at the 21 positions 17-37 on the 7th read of the Mother read set. How can be k-mer coverage so low, any ideas?</summary>

```bash
Profex -z Mother 7 | head -20
```

That will be a sequencing error!
</details>


## Perform k-mer table operations with Logex

`Logex` combines k-mer tables using logical expressions. In the examples below, `A` refers to the first input table, `B` to the second, and `C` to the third.

| Operator | Meaning |
|---|---|
| `A \| B` | Union: k-mers present in A or B |
| `A & B` | Intersection: k-mers present in both A and B |
| `A - B` | Difference: k-mers present in A but not B |
| `A ^ B` | Exclusive OR: k-mers present in one table but not both |

`|` and `&` must be followed by `+`, `<`, `>`, `*` or `.` to say which count to keep (e.g. `A&<B`); plain `A&B` gives a "Modeless operator" error. `-` and `^` can be used as they are. For union (`|`) and intersection (`&`), count operators specify how counts from k-mers present in both tables are combined: `|+` (sum), `|<` (minimum), `|>` (maximum), `|*` (average), `|.` (use the left count when present, otherwise the right), and `|-` (left count minus right count, with a minimum of zero). The same count modifiers can be used with intersection.

A postfix range such as `[5-10]` retains k-mers whose counts are in that range. `[-10,20-]` retains counts from 1–10 or 20 and above. The prefix `#` changes all k-mer counts in a table to 1.

Logex takes `-h` arugment to make both table `.ktab` and histogram `.hist` and `-H` to compute only the histogram. Example command structure:

```bash
Logex -H 'CM=(A|+B)' Child Mother
Histex CM
```

This writes a table named `CM` using the sum-count union of the `Child` and `Mother` tables. Use quoted expressions so the shell passes the expression to Logex as one argument. If the output table already exists, Logex asks `Output table ./out.ktab already exists, continue?` and waits for an answer, so the terminal looks frozen. Type `y` to overwrite, or use a new output name for each task.

```bash
FastK -v -k21 -t1 -p -T4 Hap1.fasta -NHap1
FastK -v -k21 -t1 -p -T4 Hap2.fasta -NHap2

Logex -H1 -T8 'N1=A&.C[5-]' 'N2=B&.C[5-]' 'NU=(A|+B)&.C[5-]' 'D1=A' 'D2=B' 'DU=A|>B' Hap1 Hap2 Child

Histex -A N1
Histex -A D1

# to divide the two numbers in R (this is the way how to do R calculation as oneliners)
R --slave -e "$(Histex -A N1 | cut -f 2) / $(Histex -A D1 | cut -f 2)"
```

hap1 is 99.98% k-mer complete assembly (counting only k-mers >5x).

## Relative profiles

FastK can profile input sequences using k-mer counts from a different table. This is useful, for example, to profile child reads against a table of k-mers specific to the mother or father. The `-k` value must match the k-mer size of the table given to `-p:` (here 21); FastK's default is 40, but we used 21.

The `-p:<table>` option specifies the table used for profiling. For example:

```bash
FastK -k21 -p:Mother Child.Reads.fasta -NChild_Mother
FastK -k21 -p:Father Child.Reads.fasta -NChild_Father
```

If it's still Monday and you got it here already, STOP. Keep the rest for tomorrow.

## Merqury.FK

### using trio datas

Start by making haplotype files using `HAPmaker`

```bash
HAPmaker -v -T4 Mother.ktab Father.ktab Child.ktab # generated Mother.hap.ktab and Father.hap.ktab
mkdir -p merqury
MerquryFK Child.ktab Father.ktab Mother.ktab Hap1.fasta Hap2.fasta merqury/merqury
```

Now, all these calculations and plotted histograms are just series of `Logex` expressions passed into simple plotting scripts. Try to go through all the plots (explained in the lecture), and perhaps try to rebuild some of the `Logex` expressions (for instance for k-mer completeness) or even you could try to recreate some of the figures in R. 

## Smudgeplot

