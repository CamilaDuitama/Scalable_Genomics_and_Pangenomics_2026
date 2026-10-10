# Practical 2 - What do we do with a pangenome?

In this practical session, we will try a few visualizations and queries that let us explore our pangenome collection. Most of these steps are follow ups to the previous section, where we built various types of indexes and representations of a "pangenome" for small datasets of genome assemblies. The next question is: what can we do with these pangenomes?

**Docker Image with all the tools pre-installed:** [npmalfoy/scalable:2026](https://hub.docker.com/r/npmalfoy/scalable) (`linux/amd64`)

**Course data:** [https://kj11.cog.sanger.ac.uk/#kmer_course_data/kmer_course_data/](https://kj11.cog.sanger.ac.uk/#kmer_course_data/kmer_course_data/)

> Expected runtimes and memory below are rough guidelines from a dry-run on an Apple Silicon M5 Pro MacBook Pro (Docker `linux/amd64`). Your machine may differ.

---

## Learning objectives

In this practical you will:

1. Explore how to visualize the outputs for pangenome exploration
2. Learn what types of queries each type of pangenome index is suitable for

---


## 0) Pull the docker container and datasets

See the course informatics guide for full setup instructions. This practical assumes you already downloaded and unpacked the course archives in Practical 1 (`datasets_agc.tar.gz`, `course_data.tar.gz`, `alignments_athaliana.tar.gz`). As a reminder, your working directory should look like:

```text
datasets/          # AGC assemblies (from datasets_agc.tar.gz)
data/              # pre-built indexes, reads, MHC graph (from course_data.tar.gz)
  athaliana/
    alignments/    # pairwise PAFs (from alignments_athaliana.tar.gz)
    tool_outputs/  # mumemto / ropebwt3 / syng / …
  human/           # MHC graph, Giraffe indexes, HiFi reads
  work/            # writable outputs from Practical 1 (and this practical)
```

Start an interactive session with the working directory mounted at `/course`:

**Docker:**

```bash
docker pull npmalfoy/scalable:2026
docker run --rm -it --platform linux/amd64 -v "$PWD:/course" -w /course npmalfoy/scalable:2026
```

**Singularity / Apptainer:**

```bash
apptainer pull scalable.sif docker://npmalfoy/scalable:2026

apptainer shell \
  --bind "$PWD:/course" \
  --pwd /course \
  scalable.sif
```

> Use `singularity` in place of `apptainer` if that is what your cluster provides. On Apple Silicon (or other arm64 hosts), keep `--platform linux/amd64` for Docker.

Inside the container, AGCs are under `/course/datasets`, pre-built outputs under `/course/data`, and writable outputs under `/course/data/work`.

To test that the tools are installed and available:

```bash
which agc ropebwt3 syng impg BandageNG panacus vg mumemto shredtools minimap2
```


---


## 1) How to use each pangenome representation

Following the previous practical, we will use the shipped *A. thaliana* indexes under `data/athaliana/` for the steps below (Practical 1’s yeast builds are a different scale). Writable outputs go under `$OUT`. Set shared paths once:

```bash
DATASETS=/course/datasets
DATA=/course/data
OUT=/course/data/work
ATHALIANA_OUT=$DATA/athaliana/tool_outputs
ATHALIANA_AGC=$DATASETS/athaliana_all.agc
ATHALIANA_FASTA_DIR=$OUT/athaliana_fastas
ANALYSIS_DIR=$OUT/analysis
mkdir -p "$OUT" "$ATHALIANA_FASTA_DIR" "$ANALYSIS_DIR"
```

First, decompress the 68 *A. thaliana* assemblies (assemblies from Lian et al., *Nat Genet* 2024) from the AGC archive into a local FASTA directory:

**Expect:** ~30 s; ~0.5–1 GB RAM.

<details>
<summary>Show AGC decompress command</summary>

```bash
agc info "$ATHALIANA_AGC"
agc listset "$ATHALIANA_AGC" | head
agc getcol -o "$ATHALIANA_FASTA_DIR" "$ATHALIANA_AGC"
ls "$ATHALIANA_FASTA_DIR" | head
```

</details>


---


### Mumemto / Shredtools

We've computed the set of multi-MUMs across the pangenome. In general for useful downstream analysis, we want the multi-MUM coverage to be >50% across the collection. Lower coverage likely indicates more dissimilar assemblies or high repeat content.

Quick check — compute MUM coverage on the shipped *A. thaliana* collection (`-i` takes the mumemto **prefix**, not the `.bumbl` path; `-s` selects which assembly index to report):

<details>
<summary>Show mumemto coverage command</summary>

```bash
mumemto coverage -i "$ATHALIANA_OUT/mumemto/mumemto" -s 0
```

</details>

<details>
<summary>Example output (coverage)</summary>

```text
seq0: 44.786%
```

Coverage around ~45% here is below the rough >50% guideline above — useful to keep in mind when interpreting synteny plots for this collection.

</details>

First, we will visualize the multi-MUM synteny of *A. thaliana* genomes. Since each assembly has the same number of contigs, we can split the visualization by chromosome using `mumemto viz --mode gapped`.

**Expect:** `mumemto viz` ~15–30 s and ~4 GB RAM; `shredtools extract` a few seconds and <0.5 GB.

<details>
<summary>Show mumemto command</summary>

```bash
mumemto viz -o "$ANALYSIS_DIR/mumemto.pdf" -i "$ATHALIANA_OUT/mumemto/mumemto" --mode gapped
```

</details>

<details>
<summary>Example output (multi-MUM gapped synteny)</summary>

![Multi-MUM gapped synteny plot across *A. thaliana* assemblies (course dry-run)](figures/p2_mumemto_gapped.png)

</details>


Next, we can query a region of interest and extract syntenic regions across the pangenome using `shredtools extract`. For this exercise, we will extract the FLC gene involved in flowering using the following region: `chr5:3,173,000–3,179,000`. Precomputed outputs live under `$ATHALIANA_OUT/mumemto/` (`mumemto.bumbl`, `mumemto.bumbl.bi`, `mumemto.lengths`).

<details>
<summary>Show shredtools command</summary>

```bash
shredtools extract -o "$ANALYSIS_DIR/flc" -s 0 -r CP138175.1:3173000-3179000 --plot "$ATHALIANA_OUT/mumemto/mumemto.bumbl"
```

</details>

<details>
<summary>Example output (FLC locus synteny extract)</summary>

![Synteny plot for the FLC region (chr5:3,173,000–3,179,000) across the pangenome (course dry-run)](figures/p2_flc_extract_synteny.png)

</details>

<details>
<summary>Example output (flc.bed)</summary>

```bash
head -n 10 "$ANALYSIS_DIR/flc.bed"
```

```text
CP138175.1	3173000	3178999	/course/data/work/athaliana_fastas/GCA_036942435.1_ASM3694243v1_genomic.fa
CM072661.1	3200189	3206197	/course/data/work/athaliana_fastas/GCA_036926925.1_ASM3692692v1_genomic.fa
CM072641.1	3205318	3211334	/course/data/work/athaliana_fastas/GCA_036926965.1_ASM3692696v1_genomic.fa
CM072656.1	3206991	3213009	/course/data/work/athaliana_fastas/GCA_036926975.1_ASM3692697v1_genomic.fa
CM072646.1	3184353	3190390	/course/data/work/athaliana_fastas/GCA_036927025.1_ASM3692702v1_genomic.fa
CM072651.1	3187017	3193043	/course/data/work/athaliana_fastas/GCA_036927045.1_ASM3692704v1_genomic.fa
CM072631.1	3170844	3176858	/course/data/work/athaliana_fastas/GCA_036927085.1_ASM3692708v1_genomic.fa
CM072626.1	3205639	3211652	/course/data/work/athaliana_fastas/GCA_036927245.1_ASM3692724v1_genomic.fa
CM072636.1	3169581	3175586	/course/data/work/athaliana_fastas/GCA_036927255.1_ASM3692725v1_genomic.fa
CM072621.1	3194921	3200928	/course/data/work/athaliana_fastas/GCA_036927265.1_ASM3692726v1_genomic.fa
```

</details>

We will use the MHC locus later in this practical for the Giraffe mapping comparison. For now, you can explore synteny of that locus interactively in the [shredtools HPRC browser](https://vikshiv.github.io/shredtools/hprc/) with precomputed MUMs from HPRC2 — try querying CHM13 `chr6:28381448-33301940` to see the homologous interval on I002C (`chr6:28477481-33546197`) and across other haplotypes.

---


### ropebwt3

We previously generated an FM-index using ropebwt3. There are a few things we can do, all centered around finding exact matches between a query and the index.

The first command is `mem`. We will query a set of reads against the index and find all the MEMs (maximal exact match) that appear between a read and the index.

**Expect:** `mem` ~5–15 s and ~1 GB RAM; hapdiv prep + `hapdiv` on the short FLC extract a few seconds and ~1 GB.

<details>
<summary>Show ropebwt3 mem command</summary>

```bash
ropebwt3 mem "$ATHALIANA_OUT/ropebwt3/rb3.fmd" "$DATA/athaliana/reads/tanz1_1k.fq" > "$ANALYSIS_DIR/read_mems.txt"
```

</details>

<details>
<summary>Example output (read_mems.txt)</summary>

Column legend ([ropebwt3 MEM output](https://github.com/lh3/ropebwt3#mem)):

| Col | Meaning |
| --- | ------- |
| 1 | query sequence name |
| 2 | match start on the query (0-based) |
| 3 | match end on the query |
| 4 | number of hits in the index |

```bash
head -n 10 "$ANALYSIS_DIR/read_mems.txt"
```

```text
m64079_211006_084742/18/ccs	0	452	12
m64079_211006_084742/18/ccs	439	458	2
m64079_211006_084742/18/ccs	448	470	47
m64079_211006_084742/18/ccs	450	471	40
m64079_211006_084742/18/ccs	453	724	1
m64079_211006_084742/18/ccs	467	983	2
m64079_211006_084742/18/ccs	979	1564	2
m64079_211006_084742/18/ccs	1565	2235	1
m64079_211006_084742/18/ccs	1603	2301	6
m64079_211006_084742/18/ccs	2300	3218	5
```

</details>


Next, we can query an assembly against the index and compute kmer diversity with respect to the pangenome. This is helpful to identify regions that are highly similar or dissimilar in a query assembly with respect to the population. For this, we can use the `ropebwt3 hapdiv` command.

We've held out Tanz-1 from the ropebwt3 index (`$DATA/athaliana/holdout/Tanz-1.fa`), but it is included in the mumemto collection so we can recover its FLC locus from the BED produced above. Pull that sequence with `shredtools fasta`, then run `hapdiv` on the short FLC FASTA (querying the full Tanz-1 genome is slow on a laptop).

<details>
<summary>Show ropebwt3 hapdiv command</summary>

```bash
grep Tanz-1 "$ANALYSIS_DIR/flc.bed" > "$ANALYSIS_DIR/flc_tanz.bed"
shredtools fasta -o "$ANALYSIS_DIR/flc_fa" "$ANALYSIS_DIR/flc_tanz.bed"
ropebwt3 hapdiv "$ATHALIANA_OUT/ropebwt3/rb3.fmd" "$ANALYSIS_DIR/flc_fa/"*extract*.fa > "$ANALYSIS_DIR/hapdiv.txt"
```

</details>

<details>
<summary>Example output (hapdiv.txt)</summary>

Column legend ([ropebwt3 hapdiv](https://github.com/lh3/ropebwt3#e2e)): each row is a 101-mer window on the query.

| Col | Meaning |
| --- | ------- |
| 1 | query name |
| 2 | query start |
| 3 | query end |
| 4 | number of distinct alleles the 101-mer matches |
| 5 | maximum edit distance observed |
| 6 | haplotypes with a perfect match |
| 7–11 | haplotypes at edit distance 1–5 |
| 12 | haplotypes at edit distance ≥6 |

```bash
head -n 10 "$ANALYSIS_DIR/hapdiv.txt"
```

```text
CP138120.1:3204073-3210103	0	101	2	2	0	67	1	0	0	0	0
CP138120.1:3204073-3210103	50	151	2	3	0	66	0	2	0	0	0
CP138120.1:3204073-3210103	100	201	6	2	54	13	1	0	0	0	0
CP138120.1:3204073-3210103	150	251	8	2	48	18	2	0	0	0	0
CP138120.1:3204073-3210103	200	301	4	2	59	8	1	0	0	0	0
CP138120.1:3204073-3210103	250	351	2	1	65	3	0	0	0	0	0
CP138120.1:3204073-3210103	300	451	3	1	64	4	0	0	0	0	0
CP138120.1:3204073-3210103	400	501	3	3	63	1	0	4	0	0	0
CP138120.1:3204073-3210103	450	551	5	5	0	18	44	2	0	4	0
CP138120.1:3204073-3210103	500	601	6	4	0	19	24	20	5	0	0
```

</details>

---


### syng

We previously built a syncmer dictionary (`.1khash`) and GBWT (`.1gbwt`) with syng. Analogous to ropebwt3's MEMs over bases, `syngmap` finds MEMs over syncmers between a query read set and the pangenome GBWT.

**Expect:** ~5–15 s; ~4 GB RAM.

<details>
<summary>Show syngmap command</summary>

```bash
syngmap -o "$ANALYSIS_DIR/syngmap" -outputIds "$ATHALIANA_OUT/syng/syng.1khash" "$ATHALIANA_OUT/syng/syng.1gbwt" "$DATA/athaliana/reads/tanz1_1k.fq"
```

</details>

<details>
<summary>Example output — header schema (syngmap.1map)</summary>

```bash
ONEview "$ANALYSIS_DIR/syngmap.1map" | head -n 25
```

```text
1 3 map 2 1
! 4 7 syngmap 3 2.0 217 /opt/tools/bin/syngmap -o /course/data/work/analysis/syngmap -outputIds /course/data/athaliana/tool_outputs/syng/syng.1khash /course/data/athaliana/tool_outputs/syng/syng.1gbwt /course/data/athaliana/reads/tanz1_1k.fq 19 2026-10-06_18:55:56
! 4 7 ONEview 3 0.0 83 /opt/conda/envs/scalable_course/bin/ONEview /course/data/work/analysis/syngmap.1map 19 2026-10-10_23:03:49
.
< 52 /course/data/athaliana/tool_outputs/syng/syng.1khash 1
< 51 /course/data/athaliana/tool_outputs/syng/syng.1gbwt 2
< 40 /course/data/athaliana/reads/tanz1_1k.fq 3
.
~ O S 2 3 INT 3 INT         query sequence: index in source file (1-based) length
~ D I 1 6 STRING            identifier from source file (if requested)
~ D F 1 4 CHAR              filter: Z zero-length, Q quality, G poly-G (Illumina bad read)
~ D M 3 3 INT 3 INT 3 INT   mem: start, end (0-based), count
~ D X 2 3 INT 3 DNA         missing syncmer not found in graph: start coordinate, sequence
~ D U 3 3 INT 3 INT 3 INT   unique mapping: file, path, offset (negative offset = reverse)
```

</details>

<details>
<summary>Example output — first alignment / data lines (syngmap.1map)</summary>

```bash
ONEview "$ANALYSIS_DIR/syngmap.1map" | awk '/^S /{p=1} p' | head -n 16
```

After the schema (`~` lines) and counts, each query appears as an `S`/`I` block with syncmer MEMs (`M`), missing syncmers (`X`), and unique mappings (`U`):

```text
S 1 13523
I 27 m64079_211006_084742/18/ccs
M 1 429 13
X 416 63 aacaagtatattaagaaacgattttgtttcattgataaatcagaaaaaaacgaaatagaaatc
M 460 724 1
M 506 965 2
X 936 63 aagaatattgcacaaactaattaaaatatccaaattatcacaattttatctcttaatcgattt
X 939 63 aatattgcacaaactaattaaaatatccaaattatcacaattttatctcttaatcgatttatt
X 940 63 atattgcacaaactaattaaaatatccaaattatcacaattttatctcttaatcgatttatta
M 984 1529 2
X 1521 63 atagctctgatttgagcaataagcgtgtctctgatggaattatctctgggtttgggagtgaga
M 1568 2235 1
M 1611 2301 6
X 2293 63 cgatccgctcattcttttcagcttgagactgtttcagagaaaccaatagttggtaagatgcag
M 2305 3192 5
X 3168 63 cccttcagtagcatatgcacaatgccgagaactaagtggattggtggatactgttaacataca
```

</details>

---


### impg

We can project the same FLC gene region through the all-vs-all alignments, extract the homologous sequences across the pangenome, and build a local variation graph with `impg`.

**Expect:** FLC `impg query` ~30 s; ~0.5 GB RAM. (Rebuilding the `.impg` index, if needed, is a few seconds.)

```bash
IMPG_PAF_DIR=$DATA/athaliana/alignments
IMPG_PAF_LIST=$OUT/ath_pafs.txt
IMPG_IDX=$OUT/alignments.paf.impg
```

This step needs the *A. thaliana* pair PAFs (`alignments_athaliana.tar.gz` → `$IMPG_PAF_DIR`) and the `.impg` index built in Practical 1 (`$IMPG_IDX`). Confirm both are present:

```bash
ls -lh "$IMPG_PAF_DIR"/*.paf | head
ls -lh "$IMPG_IDX"
```

If the alignments directory is missing, unpack that archive into your course directory (outside the container) and restart / remount so `/course/data/athaliana/alignments/` is visible.

If the index is missing (you did not finish the `impg index` step in Practical 1), rebuild it:

```bash
find "$IMPG_PAF_DIR" -name '*.paf' | sort > "$IMPG_PAF_LIST"
impg index --alignment-list "$IMPG_PAF_LIST" -i "$IMPG_IDX"
```

Task: Query the FLC region (`CP138175.1:3173000-3179000`) from the *A. thaliana* alignments and build a GFA graph of the homologous sequences.

<details>
<summary>Show impg command</summary>

```bash
find "$IMPG_PAF_DIR" -name '*.paf' | sort > "$IMPG_PAF_LIST"
impg query -i "$IMPG_IDX" --alignment-list "$IMPG_PAF_LIST" -r "CP138175.1:3173000-3179000" -d 1000 -x -o gfa --sequence-files "$ATHALIANA_FASTA_DIR"/*.fa -O "$ANALYSIS_DIR/flc"
```

</details>

<details>
<summary>Example output (flc.gfa)</summary>

```bash
head -n 10 "$ANALYSIS_DIR/flc.gfa"
```

```text
H	VN:Z:1.0
S	1	CTTAAAAGAAAAAAATAGGTACAATTT
S	2	A
S	3	T
S	4	CAAAAACATTCCTATATTATGCGCAAAGTGACAAGAAAAAGAAACAGCCCTATCAAAAGTTTTAGTTAG
S	5	AT
S	6	ATATATTAATTAAGTCGGTGCAATCAACTCTCGTTATCCATCGAGTATGATGCAACCGAAGAACAAG
S	7	T
S	8	G
S	9	T
```

</details>


Inputs: all-vs-all pairwise alignments (PAF), the corresponding FASTA sequences, and a query interval

Outputs: extracted homologous FASTA sequences aligned into a local GFA graph (`flc.gfa`)

---


### BandageNG

[BandageNG](https://github.com/asl/BandageNG) visualises assembly / variation graphs. Here we render the FLC graph built with impg above.

**Expect:** a few seconds; <0.1 GB RAM.

Task: Produce an image of the FLC variation graph.

<details>
<summary>Show BandageNG command</summary>

```bash
# thicker node/edge settings (same idea as Practical 4) so the SVG stays readable
BANDAGE_OPTS=(
  --nodelen 50
  --minnodlen 80
  --edgewidth 40
  --nodewidth 200
  --depwidth 0
  --outline 8
  --colour uniform
  --unicolpos '#08306b'
  --unicolneg '#08306b'
  --edgecol '#000000'
  --outcol '#000000'
)
BandageNG image "${BANDAGE_OPTS[@]}" "$ANALYSIS_DIR/flc.gfa" "$ANALYSIS_DIR/flc.svg"
```

</details>

Optionally, thicken strokes further for readability. The Bandage settings above already draw thicker nodes/edges into the SVG; if you also want a thickened raster (e.g. for slides), use the same helper script from Practical 4 (`p4_thicken_png.py` — Pillow is already in the course image):

<details>
<summary>Optional: thicken a Bandage PNG</summary>

```bash
wget -q -O "$ANALYSIS_DIR/p4_thicken_png.py" \
  https://raw.githubusercontent.com/WCSCourses/Scalable_Genomics_and_Pangenomics_2026/main/course_modules_2026/module_5_pangenome_practicals/scripts/p4_thicken_png.py
BandageNG image "${BANDAGE_OPTS[@]}" "$ANALYSIS_DIR/flc.gfa" "$ANALYSIS_DIR/flc.png"
python3 "$ANALYSIS_DIR/p4_thicken_png.py" --rounds 1 \
  -i "$ANALYSIS_DIR/flc.png" -o "$ANALYSIS_DIR/flc.thick.png"
```

</details>

<details>
<summary>Example output (FLC variation graph layout, thickened SVG)</summary>

![BandageNG layout of the FLC impg GFA (thickened)](figures/p2_flc_bandage.svg)

</details>


Inputs: a variation graph in GFA format

Outputs: an SVG image of the graph layout (`flc.svg`; optionally a thickened PNG too)

---


### panacus

[panacus](https://github.com/codialab/panacus) computes coverage and growth statistics over a pangenome graph (GFA). This is useful for determining how much of the pangenome is core vs accessory based on the multiple alignment encoded in the graph topology.

**Expect:** a few seconds; <0.1 GB RAM.

Task: Build an HTML report with coverage histogram and growth curves for the FLC graph from impg.

<details>
<summary>Show panacus example command</summary>

```bash
cat > "$ANALYSIS_DIR/flc_report.yaml" <<EOF
- !Gfa
  graph: $ANALYSIS_DIR/flc.gfa
  count_type: Bp
  analyses:
    - !Hist
    - !Growth
      coverage: 1,1,2
      quorum: 0,0.9,0
EOF

panacus report "$ANALYSIS_DIR/flc_report.yaml" > "$ANALYSIS_DIR/flc_report.html"
```

</details>

<details>
<summary>Example output (panacus growth plot)</summary>

![Panacus growth curves for the FLC GFA (course dry-run)](figures/p2_flc_panacus_report.png)

</details>


Inputs: a variation graph in GFA format

Outputs: an interactive HTML report (`flc_report.html`) with histogram and growth analyses

---


### vg giraffe

In Practical 1 we built two Giraffe indexes for the MHC region (prefer your own builds under `$OUT/vg_giraffe/`; fall back to the shipped indexes under `$DATA/human/vg_giraffe/`):

**Expect:** MHC graph mapping ~15–30 s and ~0.5–1 GB RAM; CHM13 linear mapping ~5–15 s and ~0.5 GB.

```bash
VG_DIR=$OUT/vg_giraffe
# VG_DIR=$DATA/human/vg_giraffe   # fallback: pre-built
```
- `$VG_DIR/mhc` — variation graph over the 5 human haplotypes
- `$VG_DIR/mhc_chm13` — linear CHM13 sequence for the same region

Task: Align the provided human HiFi long reads with [`vg giraffe`](https://github.com/vgteam/vg/wiki/Mapping-long-reads-with-Giraffe) against **both** indexes and compare the resulting alignments.

<details>
<summary>Show solution</summary>

```bash
# map against the 5-haplotype MHC graph (long-read indexes + hifi preset)
vg giraffe -b hifi \
  -Z "$VG_DIR/mhc/mhc.giraffe.gbz" \
  -d "$VG_DIR/mhc/mhc.dist" \
  -m "$VG_DIR/mhc/mhc.longread.withzip.min" \
  -z "$VG_DIR/mhc/mhc.longread.zipcodes" \
  -f "$DATA/human/reads/i002c_mhc_1k.fa" > "$ANALYSIS_DIR/reads.mhc.gam"

# map against the linear CHM13 MHC index
vg giraffe -b hifi \
  -Z "$VG_DIR/mhc_chm13/mhc_chm13.giraffe.gbz" \
  -d "$VG_DIR/mhc_chm13/mhc_chm13.dist" \
  -m "$VG_DIR/mhc_chm13/mhc_chm13.longread.withzip.min" \
  -z "$VG_DIR/mhc_chm13/mhc_chm13.longread.zipcodes" \
  -f "$DATA/human/reads/i002c_mhc_1k.fa" > "$ANALYSIS_DIR/reads.mhc_chm13.gam"
```

</details>

Inputs: long-read Giraffe indexes from Practical 1 (`lr-giraffe`) and long HiFi reads

Outputs: GAM alignments for each index

#### compare the two alignments (example of reference bias)

The HiFi query reads are **I002C paternal**, not CHM13. Aligning them to a **single linear CHM13** path may contain *reference bias* setup: where I002C containes variation absent from CHM13, alignments to a linear index may be poor quality (or not possible). Including multiple human haplotypes increases the variation in the graph, so the aligner could choose a path that fits the read better, yielding a better alignment.

**Expect:** a few seconds per `vg stats -a` call; negligible RAM.

Task: Summarize both GAMs with `vg stats -a`, then skim the differences in aligned fraction, alignment score, and soft-clip burden.

```bash
# VG_DIR as set above for vg giraffe
echo "=== MHC pangenome graph ==="
vg stats -a "$ANALYSIS_DIR/reads.mhc.gam" "$VG_DIR/mhc/mhc.giraffe.gbz"

echo "=== linear CHM13 MHC ==="
vg stats -a "$ANALYSIS_DIR/reads.mhc_chm13.gam" "$VG_DIR/mhc_chm13/mhc_chm13.giraffe.gbz"
```

<details>
<summary>Example output (vg stats -a)</summary>

```text
=== MHC pangenome graph ===
Total alignments: 1000
Total primary: 1000
Total secondary: 0
Total aligned: 995
Total perfect: 5
Total gapless (softclips allowed): 10
Total paired: 0
Total properly paired: 0
Alignment score: mean 3342.85, median 2651, stdev 3044.03, max 22159 (1 reads)
Mapping quality: mean 2.90553, median 2, stdev 7.20156, max 60 (13 reads)
Insertions: 52432 bp in 9628 read events
Deletions: 62076 bp in 15304 read events
Substitutions: 256484 bp in 230556 read events
Matches: 5914461 bp (5944.18 bp/aligned)
Softclips: 9132670 bp (59.4728% of bases, 9178.56 bp/aligned) in 1700 read events
Total time: 64.1198 seconds
Speed: 15.5958 reads/second
Unvisited nodes: 505720/555723 (10084387 bp)
Single-visited nodes: 33685/555723 (736546 bp)
Significantly biased heterozygous sites: 0/0

=== linear CHM13 MHC ===
Total alignments: 1000
Total primary: 1000
Total secondary: 0
Total aligned: 995
Total perfect: 0
Total gapless (softclips allowed): 3
Total paired: 0
Total properly paired: 0
Alignment score: mean 2643.13, median 2203, stdev 2258.73, max 18033 (1 reads)
Mapping quality: mean 5.90151, median 2, stdev 13.1715, max 60 (33 reads)
Insertions: 57256 bp in 10302 read events
Deletions: 80161 bp in 15740 read events
Substitutions: 293367 bp in 263003 read events
Matches: 5846585 bp (5875.96 bp/aligned)
Softclips: 9158839 bp (59.6432% of bases, 9204.86 bp/aligned) in 1709 read events
Total time: 27.2562 seconds
Speed: 36.6888 reads/second
Unvisited nodes: 123508/153766 (3952237 bp)
Single-visited nodes: 21931/153766 (701792 bp)
Significantly biased heterozygous sites: 0/0
```

</details>

**Discuss:**

- Total **aligned** counts are often similar. What does that imply about *whether* reads map, versus *how well* they map?
- Mean **alignment score** and **perfect / gapless** counts: which index looks more confident for these reads, and why might that track donor haplotype rather than “CHM13 correctness”?

---


## 2) Optional: build an MHC graph with impg and compare to [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)

This optional exercise rebuilds an MHC variation graph with `impg` and compares it to the shipped [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) MHC graph using BandageNG and graph stats.

[impg](https://github.com/pangenome/impg) indexes and enables querying of genomic intervals across a pangenome using all-vs-all pairwise genome alignments. Running all-pairs alignments is not practical on a laptop, so we ship pre-computed human alignments (`alignments_human.tar.gz` — optional large download from the [course data bucket](https://kj11.cog.sanger.ac.uk/#kmer_course_data/kmer_course_data/); unpacks to `data/human/alignments.paf`). Do **not** run `impg align` here. A companion `.impg` index for the shipped human PAF is already in `course_data` as `$DATA/human/alignments.paf.impg`. Skip this section if you did not download `alignments_human.tar.gz`.

Task: Using the shipped index and alignments, build a graph of the MHC region across the set of human genomes ([CHM13 coords](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/CHM13@all/OtherDifficult/CHM13_MHC.bed.gz) -> chr6:28381448-33301940). Use `-o gfa:seqwish`, which skips the slow smoothing step for the purposes of this tutorial (graph induction via [seqwish](https://github.com/pangenome/seqwish)). Then visualise the impg graph and the shipped [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) graph side by side with BandageNG, and compare basic graph statistics.

This step needs the human assemblies. Decompress them if you have not already, and confirm the optional human PAF is present:

```bash
AGC=$DATASETS/human.agc
FASTA_DIR=$OUT/human_fastas
mkdir -p "$FASTA_DIR"

# [wfmash](https://github.com/waveygang/wfmash) ships with impg under libexec — symlink into conda bin for --aligner wfmash
ln -sfn /opt/conda/envs/scalable_course/libexec/impg/wfmash /opt/conda/envs/scalable_course/bin/wfmash

ls -lh "$DATA/human/alignments.paf" "$DATA/human/alignments.paf.impg"
agc getcol -o "$FASTA_DIR" "$AGC"
```

<details>
<summary>Show impg query command</summary>

```bash
impg query -i "$DATA/human/alignments.paf.impg" -a "$DATA/human/alignments.paf" -r "chr6:28381448-33301940" -d 1000 -x -o gfa:seqwish --sequence-files "$FASTA_DIR"/*.fa -O "$ANALYSIS_DIR/mhc_impg" --aligner wfmash
```

</details>

Inputs: all-vs-all pairwise alignments (PAF) and the corresponding FASTA sequences

Outputs: a GFA graph for the queried MHC region (e.g. `$ANALYSIS_DIR/mhc_impg.gfa`)

### BandageNG and stats: impg vs [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)

The course ships a smoothed minigraph-cactus MHC graph as `$DATA/human/mhc/mhc.full.gfa.gz`. After the impg query finishes, render both graphs and print summary stats so you can compare topology / size at a glance.

<details>
<summary>Show BandageNG and stats commands</summary>

```bash
# unpack the shipped MC graph next to the impg output
gunzip -c "$DATA/human/mhc/mhc.full.gfa.gz" > "$ANALYSIS_DIR/mhc_mc.full.gfa"

# visualise both
BandageNG image "$ANALYSIS_DIR/mhc_impg.gfa" "$ANALYSIS_DIR/mhc_impg.svg"
BandageNG image "$ANALYSIS_DIR/mhc_mc.full.gfa" "$ANALYSIS_DIR/mhc_mc.svg"

# basic graph stats (nodes / edges / length)
vg stats -z "$ANALYSIS_DIR/mhc_impg.gfa"
vg stats -z "$ANALYSIS_DIR/mhc_mc.full.gfa"
```

</details>

Inputs: the impg MHC GFA and the shipped minigraph-cactus MHC GFA

Outputs: SVG layouts for both graphs plus printed `vg stats` summaries to compare

---
