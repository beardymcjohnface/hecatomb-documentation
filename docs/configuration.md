## Recommended customisation

If you're running Hecatomb on a HPC cluster, we absolutely recommend setting up a 
[Snakemake profile](profiles.md).

## Changing the Hecatomb configuration

The Hecatomb configuration file contains settings related to resources and cutoffs for various stages of the pipeline. 
The different config settings are outlined further on.
You can permanently change the behaviour of your Hecatomb installation by modifying the values in the system config file.

Alternatively, you can copy and modify the default config file.
Before Hecatomb runs, it will copy the system default config file to the output directory and use it for your analysis.
To customise your run, you can copy the system default config file before running Hecatomb like so:

```shell
hecatomb config
```

You can then edit your new `hecatomb.out/hecatomb.config.yaml` file to suit your needs.
It will automatically be used in your Hecatomb run, or if you rename it you can specify the file with `--configfile`:

```shell
hecatomb run --configfile myRenamedHecatomb.config.yaml
```

## Database location

The databases are large (~55 GB) and if your Hecatomb installation is on a partition with limited on space,
you might want to specify a new location to house the database files. 
By default, this config setting is blank and the pipeline will use the install location.
You can specify the directory in the Hecatomb system config file under `args: databases: `, 
e.g:

```yaml
args:
    databases: /scratch/HecatombDatabases
```

and rerun the installation 

```shell
hecatomb install
```

## Default resources

The Hecatomb config file contains some sensible defaults for resources.
While these should work for most datasets, they may fail for larger ones.
You may also have more CPUs etc at your disposal and want to minimise runtime of the pipeline.
Currently, the slowest steps are the MMSeqs searches (under `resources: big:`); 
increasing the CPUs and RAM can significantly improve runtime.
The other settings (under `resources: med:`) will yield more modest improvement.

```yaml
resources:
    big:
        mem: 64000     # Memory for big (mostly mmseqs) jobs in megabytes (e.g 64GB = 64000, recommend >= 64000)
        cpu: 24        # Threads (recommend >= 16)
        time: 1440     # Max runtime in minutes (allows to set max time for the scheduler via snakemake profiles)
    med:
        mem: 32000      # Memory for most jobs in megabytes (recommend >= 32000)
        cpu: 16         # CPUs for most jobs in megabytes (recommend >= 16)
    ram:
        mem: 16000    # Memory for slightly RAM-hungry jobs in megabytes (recommend >= 16000)
        cpu: 2        # CPUs for slightly RAM-hungry jobs (recommend >= 2)
```

## Preprocessing settings

Fastp settings can be modified in the config under `qc: fastp:`. 
Refer to Fastp's documentation for more details.

```yaml
qc:
    compression: 1      # Compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest. Default is 1, based on assumption of large scratch space
    fastp:
        --qualified_quality_phred 15
        --length_required 90
        --cut_tail 
        --cut_tail_window_size 25
        --cut_tail_mean_quality 15
        --dedup
        --dup_calc_accuracy 4
        --trim_poly_x
```

## Alignment filtering

Hecatomb has settings for the various MMSeqs2 steps.
Read-clustering parameters using mmseqs' linclust can be modified under `mmseqs: linclustParams:`.
For annotation with MMSeqs2, alignment filtering can be modified for each of the primary (viral) and secondary
(multi-kingdom) searches under `filtAAprimary`, `filtAAsecondary`, `filtNTprimary`, and `filtNTsecondary`.
Search sensitivity parameters when running Hecatomb with `--search fast` are specified under `perfAAfast` and `perfNTfast`,
and with the default `--search sensitive` under `perfAA` and `perfNT`.
Refer to the MMSeqs2 documentation for more information on these settings.
Finally, you can specify taxIDs to ignore (i.e. defer to the best hit) under `mmseqs: taxIdIgnore:`.

```yaml
mmseqs:
    linclustParams:
        --kmer-per-seq-scale 0.3
        -c 0.8
        --cov-mode 1
        --min-seq-id 0.97
        --alignment-mode 3
    filtAAprimary:
        --min-length 30
        -e 1e-3
    filtAAsecondary:
        --min-length 30
        -e 1e-5
    filtNTprimary:
        --min-length 90
        -e 1e-10
    filtNTsecondary:
        --min-length 90
        -e 1e-20
    perfAA:
        --start-sens 1
        --sens-steps 3
        -s 7
        --lca-mode 2
        --shuffle 0
    perfAAfast:
        -s 4.0
        --lca-mode 2
        --shuffle 0
    perfNT:
        --start-sens 2
        -s 7
        --sens-steps 3
    perfNTfast:
        -s 4.0
    taxIdIgnore: 0 1 2 10239 131567 12429 2759
```

## Assembly settings

All assembly settings are specified under `assembly:`.
Canu assembler is used for `--preprocess longreads` and megahit is used for shortreads.
Flye is used to collapse individual sample contigs into a population assembly.
Refer to Canu, Flye, and Megahit's documentation for more information on these settings.

```yaml
assembly:
    canu:
        correctedErrorRate=0.16
        maxInputCoverage=10000
        minInputCoverage=0
        corOutCoverage=10000
        corMhapSensitivity=high
        corMinCoverage=0
        useGrid=False
        stopOnLowCoverage=False
        genomeSize=10M
        -nanopore
        # -pacbio
        # -pacbio-hifi
    megahit:
        --k-min 45
        --k-max 225
        --k-step 26
        --min-count 2
        --min-contig-len 1000
    flye:
        -g 1g
```


## Additional Snakemake commands

As mentioned, Hecatomb is powered by Snakemake but runs via a launcher for your convenience.
Snakemake itself has many command line options, and the launcher can pass additional commands on to Snakemake.

One such example is if you're not production ready you might wish to do a 'dry-run', where the run is simulated but no 
jobs are submitted, just to see if everything is configured correctly.
To do that, Snakemake needs the dry run flag (`--dry-run`, `--dryrun`, or `-n`).
In Hecatomb, simply tack it on to the end of your command:

```shell
hecatomb run --reads fasq/ --profile slurm --dry-run
```

Hecatomb prints the Snakemake command to the terminal window before running and you should see these additional options 
added to the Snakemake command. Have a look at the full list of available Snakemake options with `snakemake --help`. 
__Any unrecognised command will be passed on to Snakemake verbatim__, so use with caution :p

