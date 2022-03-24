# Mulberry
Mulberry is (currently) an R script used to analyze and vizualize AMR gene data.

### Prerequisites

What dependencies you need to run Mulberry:


- [NextFlow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
- data!

### Using the pipeline
The pipeline is designed to start from short read fastq files. All files must be in the same directory. Start the pipeline using:

```
$ nextflow run ~/mulberry/mulberry.nf --reads /path/to/reads
```

#### To rename the output pdf

```
$ nextflow run ~/mulberry/mulberry.nf  --title "Title name in quotes" --reads /path/to/reads
```

### Output files
The default directory for Mulberry output data is mulberry_results, unless changed by using the ```--outdir``` flag:
```
$ nextflow run ~/mulberry/mulberry.nf --outdir mulberry_results_2 --reads /path/to/reads
```


### Author

* **Rachael St. Jacques** - *Bioinformatics Principal Scientist* - [stjacqrm](https://github.com/stjacqrm)
