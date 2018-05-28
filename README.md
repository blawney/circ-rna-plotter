# circ-rna-plotter
Visualizations for circRNA transcripts compatible with KNIFE algorithm output files

Example output from a potential circRNA transcript featuring 3 exons of the MARK1 gene:

![alt text](https://github.com/blawney/circ-rna-plotter/raw/master/chr1%7CMARK1:220792097%7CMARK1:220771672%7Crev%7C%2B.png)


### Usage:

Note that you will need to use python3 and have the appropriate dependencies installed.  There is a `requirements.txt` file provided for pip3.

Use the `-h` arg to print the arguments:

```
usage: circular_transcript_plotter.py [-h] -j JUNCTION_FILE -r READS_FILE -g
                                      GTF_FILE -n READ_COUNT_THRESHOLD -p
                                      JUNCTION_CDF_THRESHOLD

optional arguments:
  -h, --help            show this help message and exit
  -j JUNCTION_FILE, --junctions JUNCTION_FILE
                        Path to the junction output file. This file ends with
                        "__circJuncProbs.txt"
  -r READS_FILE, --reads READS_FILE
                        Path to the reads output file, which has information
                        about the individual reads and their junction info.
                        This file ends with "__output.txt"
  -g GTF_FILE, --gtf GTF_FILE
                        Path to a GTF file (e.g. from UCSC or Ensembl)
  -n READ_COUNT_THRESHOLD, --num_reads READ_COUNT_THRESHOLD
                        The required minimum number of supporting reads
  -p JUNCTION_CDF_THRESHOLD, --cdf JUNCTION_CDF_THRESHOLD
                        The required minimum threshold for the confidence as
                        reported via the CDF. Think of this like a p-value,
                        except you want HIGH values. The publication considers
                        values greater than 0.9.

```
Using the provided files in this repository:
```
python3 circular_transcript_plotter.py \
  -j demo__circJuncProbs.txt \
  -r demo__output.txt \
  -g demo.gtf \
  -n 5 \
  -p 0.95
```
