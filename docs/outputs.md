The pipeline outputs a number of files for further analysis and exploration, as well as to provide an overview of the 
read preprocessing and distribution.

## Report

`report.html`

This file is generated by Snakemake and outlines a lot of information relating to the Hecatomb run.
Under the Results tabs are summary files for things like reads for each sample following different preprocessing steps
as well as some summary plots.

## SeqTable

`hecatomb_out/RESULTS/seqtable.fasta`

The SeqTable is the primary output of the read preprocessing and serves as the input for Taxonomic assignment.
It is composed of all the representative sequences from the clustered reads for all samples.
Samples are clustered individually, and the seq IDs for this fasta file follows the format `>sampleID:count:seqNumber`.
Here, `count` is the number of reads in that cluster which is important for statistical exploration.
Sequences are numbered sequentially (`seqNumber`) to ensure unique IDs.

## BigTable

`hecatomb_out/RESULTS/bigtable.tsv`

The BigTable is the main output of Taxonomic assignment and can be directly imported into R or Python.
The BigTable combines the seqtable IDs with their sampleID, counts, normalised counts, alignment information, taxonomic assignments and Baltimore classification.
This file is big, hence the name, but is designed to make merging with sample metadata, plotting, and statistical interrogation as easy as possible.

The header looks like this:

```text
seqID  sampleID  count  percent  alnType  targetID  evalue  pident  fident  nident  mismatches  qcov  tcov  qstart  qend  qlen  tstart  tend  tlen  alnlen  bits  targetName  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

<details>
    <summary>Column header definitions</summary>
    <b>seqID:</b> Sequence ID (format = sampleID:count:uniqueInt)<br>
    <b>sampleID:</b> The sample IDs derived from the read files<br>
    <b>count:</b> The number of reads represented by the sequence<br>
    <b>percent:</b> Percent of the host-removed reads (normalised count)<br>
    <b>alnType:</b> Type of alignment (aa = amino acid, nt = nucleotide)<br>
    <b>targetID:</b> The UniProt or NCBI ID of the database target sequence<br>
    <b>evalue:</b> expect value of the alignment (less is better)<br>
    <b>pident:</b> percent identity (of the alignment)<br>
    <b>fident:</b> fraction identity<br>
    <b>nident:</b> number of identical bases/residues<br>
    <b>mismatches:</b> number of mismatched bases/residues<br>
    <b>qcov:</b> coverage of query sequence (query = seqtable sequence)<br>
    <b>tcov:</b> coverage of target sequence (target = database sequence)<br>
    <b>qstart:</b> query start position<br>
    <b>qend:</b> query end position<br>
    <b>qlen:</b> query sequence length<br>
    <b>tstart:</b> target start position<br>
    <b>tend:</b> target end position<br>
    <b>tlen:</b> target sequence length<br>
    <b>alnlen:</b> alignment length<br>
    <b>bits:</b> bit score (more is better)<br>
    <b>targetName:</b> target sequence name<br>
    <b>taxMethod:</b> Method used to assign taxonomy (either Lowest Common Ancestor, "LCA"; or top hit sequence, "topHit")<br>
    <b>kingdom/phylum/class/order/family/genus/species:</b> Taxonomy annotations<br>
    <b>baltimoreType:</b> Baltimore classification (double/single strand, DNA/RNA, +/-)<br>
    <b>baltimoreGroup:</b> Baltimore classification group<br>
</details>

<br>

## TaxonLevelCounts

`hecatomb_report/taxonLevelCounts.tsv`

This file is derived from the BigTable and summarises the total sequence counts, for each sample, at all taxonomic levels.
The TaxonLevelCounts combines the sampleID with the taxonomic level for which the counts refer, the full taxonomic path, 
the taxon name, and the total and normalised read counts.
The purpose of this file is to expedite statistical interrogation of your data.
For instance, if you wanted to compare the numbers of say Flaviviridae reads between two groups of samples, 
those counts have already been collected, and you can simply run your analysis and plotting on the relevant slice of the table.  

The file looks something like this:

```text
sampleID    taxonLevel  taxonPath                                   taxonName       count   CPM
sample1     Kingdom     k_Bacteria                                  Bacteria        3162    3178.818
sample1     phylum      K_Viruses,p_Phixviricota                    Phixviricota    1216    1222.467
sample1     class       K_Viruses,p_Uroviricota,c_Caudoviricetes    Caudoviricetes  1234    1240.564
etc.
```

## Assembly

`hecatomb_out/RESULTS/assembly.fasta`

These are the contigs generated, unless you run Hecatomb with the `--skipAssembly` flag.
The assembly is used for producing the ContigSeqTable and ContigKrona plots, as well as the direct contig annotations.

`hecatomb_out/RESULTS/contig_count_table.tsv`

The contig count table contains the coverage information of all contigs for each sample.

```text
Sample  Contig  Length  Reads  RPKM  FPKM  SPM  AverageFold  ReferenceGC  CoveragePercentage  CoverageBases  MedianFold
```

<details>
    <summary>Column header definitions</summary>
    <b>Sample:</b> The sample IDs derived from the read files <br>
    <b>Contig:</b> The contig ID in assembly.fasta <br>
    <b>RPKM:</b> Reads Per Kilobase Million - see https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ <br>
    <b>FPKM:</b> Fragments Per Kilobase Million - see https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ <br>
    <b>SPM:</b> Sequences Per Million (counts normalized by library size) <br>
    <b>AverageFold:</b> Average read coverage of contig <br>
    <b>ReferenceGC:</b> Contig GC content <br>
    <b>CoveragePercentage:</b> Percent of contig covered by alignments <br>
    <b>CoverageBases:</b> Total bases of contig covered by alignments <br>
    <b>MedianFold:</b> Median read coverage of contig <br>
</details>

<br>

## Contig Annotations

`hecatomb_out/RESULTS/contigAnnotations.tsv`

The contig annotations follows a similar format to the bigtable.
Refer to bigtable for column definitions.

```text
contigID  evalue  pident  fident  nident  mismatch  qcov  tcov  qstart  qend  qlen  tstart  tend  tlen  alnlen  bits  target  kingdom  phylum  class  order  family  genus  species
```

## ContigSeqTable

`hecatomb_out/RESULTS/contigSeqTable.tsv`

The ContigSeqTable combines the read mapping information for the assembly with the read-based taxonomic assignments.
This file is intended to assist the user in identifying and binning assembly contigs by applying a consensus approach to contig taxonomic assignment.
The file includes the positional mapping information and can also enable investigation of more complex features such as 
chimeric contigs, recombination or horizontal transfer events.

The header looks like this:

```text
contigID  seqID  start  stop  len  qual  count  percent  alnType  taxMethod  kingdom  phylum  class  order  family  genus  species  baltimoreType  baltimoreGroup
```

<details>
    <summary>Column header definitions</summary>
    <b>contigID:</b> Contid ID from assembly.fasta<br>
    <b>seqID:</b> Sequence ID (format = sampleID:count:uniqueInt)<br>
    <b>start:</b> Alignment start position on contig<br>
    <b>stop:</b> Alignment end position on contig<br>
    <b>len:</b> Alignment length<br>
    <b>qual:</b> Alignment quality<br>
    <b>count:</b> The number of reads represented by the sequence<br>
    <b>percent:</b> Percent of the host-removed reads (normalised count)<br>
    <b>alnType:</b> Type of alignment (aa = amino acid, nt = nucleotide)<br>
    <b>taxMethod:</b> Method used to assign taxonomy (either Lowest Common Ancestor, "LCA"; or top hit sequence, "topHit")<br>
    <b>kingdom/phylum/class/order/family/genus/species:</b> Taxonomy annotations<br>
    <b>baltimoreType:</b> Baltimore classification (double/single strand, DNA/RNA, +/-)<br>
    <b>baltimoreGroup:</b> Baltimore classification group<br>
</details>

<br>


## krona.html and contigKrona.html

`hecatomb_report/krona.html`

`hecatomb_report/contigKrona.html`

The Krona plots are to assist in visual exploration of the read annotations.
krona.html is derived from the bigtable and shows the raw distribution of taxon assignments.
contigKrona.html is derived from the contigSeqTable and includes the taxon assignment method (either tophit or LCA).
The contigKrona plot helps to visualise the distributions of topHit versus LCA assigned reads as well as the 
distributions over contigs of the identified species, and the distribution of taxonomic assignments for each contig.
