# KIQ

## About

Author: Peter Menzel, pmenzel@gmail.com

KIQ is a program for counting the occurrences of a pre-defined static set of k-mers in a large collection of sequence data, for example RNA-Seq.
The counts are stored in a database file, which can be queried for retrieving the counts of a particular k-mer.

KIQ uses 32-mers comprised of the characters A, C, G, and T, which are represented as a 64bit unsigned integer.


## Installation

### Compilation from source

```
git clone https://github.com/pmenzel/kiq
cd kiq/src
make
```
After compilation, the executable `kiq` is located in the `kiq/bin/` folder.

### Support for SRA files

KIQ can be built with support for directly reading SRA files by using the libraries [ncbi-vdb](https://github.com/ncbi/ncbi-vdb) and [ngs](https://github.com/ncbi/ngs) from NCBI.

To this end, download, compile, and install both of these libraries first:
```
git clone https://github.com/ncbi/ngs
cd ngs
./configure
make
sudo make install
```
and
```
git clone https://github.com/ncbi/ncbi-vdb
cd ncbi-vdb
./configure
make
sudo make install
```
Afterwards, compile KIQ with
```
make sra
```

If the installation of `ncbi-vdb` and `ngs` cannot be done with superuser privileges, then specify an existing custom folder as installation target for `ncbi-vdb` and `ngs` when running their `./configure` scripts, for example:
```
./configure --prefix /home/username/software
```
Afterwards, compile KIQ with:
```
make sra NCBI_DIR=/home/username/software
```

## Usage

### Create k-mer index

First, KIQ needs to make an index from a fixed set of k-mers, which are read line by line from the text file specified with the option `-l` to the command `kiq index`.  
Additionally, the mandatory options `-i` and `-k` are required for specifying the file names for the k-mer index and the (initially empty) k-mer count database, which are later used for counting and querying.
```
kiq index -i kmerindex.bin -k kmercounts.bin -l kmers.txt
```

### Count indexed k-mers in sequence data

Second, the indexed k-mers are counted from sequencing data, either from FASTA/Q or SRA files.

For **FASTA/Q** files:
```
kiq db -i kmerindex.bin -k kmercounts.bin -m experiments.bin -l input.tsv -z 5
```
The file `input.tsv` is a two-column tab-separated file, in which the first column denotes the experiment or sample name and the second column contains the path to a FASTQ/A file containing the sequence data, for example:
```
Sample1	/path/to/reads/sample1/reads_1.fastq
Sample1	/path/to/reads/sample1/reads_2.fastq
Sample2	/path/to/reads/sample2/reads.fasta
```
The file `experiments.bin` will be created by KIQ and contains the experiment names and read counts.

For **SRA** files:
```
kiq sra -i kmerindex.bin -k kmercounts.bin -m experiments.bin -l input.tsv -z 5
```
Here, the second column in the file `input.tsv` contains the path to a SRA file
associated with the sample.  Since KIQ uses the NCBI SRA API, the second
column can also contain just the SRA identifier (starting with SRR), in which
case the file will be downloaded on the fly and stored in the NCBI local cache
directory. This directory can be modified using the
[vdb-config](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) tool
from the [sra-tools package](https://github.com/ncbi/sra-tools/).


### Query database by a k-mer

After counting k-mers from the FASTQ/SRA files, the index can be queried with `kiq query`, either by providing a list of query k-mers in a text
file using option `-q` or by providing a single query k-mer using option `-Q`.

For example:
```
kiq query -i kmerindex.bin -k kmercounts.bin -m experiments.bin -q query.txt
```
where query.txt contains the k-mers to be searched, one per line.

```
kiq query -i kmerindex.bin -k kmercounts.bin -m experiments.bin -Q ACGTACGTACGTACGTACGTACGTACGTACGT

```


### Acknowledgments

KIQ uses the [BBHash](https://github.com/rizkg/BBHash) library for indexing of the k-mers
and [zstr](https://github.com/mateidavid/zstr) library for reading compressed FASTA/Q files.


### License
See the file LICENSE.

