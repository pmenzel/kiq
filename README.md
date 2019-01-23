# KIQ

## About

Copyright (c) 2018,2019 Peter Menzel <pmenzel@gmail.com>

KIQ is a program for counting the occurrences of a pre-defined static set of
k-mers in a large collection of sequencing experiments, for example RNA-Seq.

Occurrence counts for each k-mer per experiment are stored in a database file,
which can then be queried for retrieving those experiments that contain a
particular k-mer.

KIQ uses 32-mers comprised of the characters A, C, G, and T, which are
represented as a 64bit unsigned integer.


## Installation

### Compilation from source

```
git clone https://github.com/pmenzel/kiq
cd kiq/src
make
```
After compilation, the executable `kiq` is located in the `kiq/bin/` folder.

### Support for SRA files

KIQ can be built with support for directly reading SRA files by using the
libraries [ncbi-vdb](https://github.com/ncbi/ncbi-vdb) and
[ngs](https://github.com/ncbi/ngs) from NCBI.

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

If the installation of `ncbi-vdb` and `ngs` cannot be done with superuser
privileges, then you can specify an existing custom folder as installation
target for `ncbi-vdb` and `ngs` when running their `./configure` scripts, for
example:
```
./configure --prefix /home/username/software
```
In this case, compile KIQ with:
```
make sra NCBI_DIR=/home/username/software
```

## Usage

### Create k-mer index

First, KIQ needs to make an index from a fixed set of k-mers, which are read
line by line from the text file specified with the option `-l` to the command
`kiq index`. Additionally, the mandatory options `-i` and `-k` are required
for specifying the file names for the k-mer index and the (initially empty)
k-mer count database, which are later used for counting and querying.
```
kiq index -i kmer_index.bin -k kiq_database.bin -l kmers.txt
```

### Count indexed k-mers in sequence data

Second, the indexed k-mers are counted from sequencing data, either from FASTA/Q or SRA files.

For **FASTA/Q** files:
```
kiq db -i kmer_index.bin -k kiq_database.bin -l input.tsv -z 5
```
The file `input.tsv` is a tab-separated file, in which the first column denotes
the experiment or sample name/ID and the second column contains the path to a
FASTQ/A file containing the sequence data.  Optionally, a third column can
contain a description of the sample, which will also be stored in the database.

For **SRA** files:
```
kiq sra -i kmer_index.bin -k kiq_database.bin -l input.tsv -z 5
```
Here, the second column in the file `input.tsv` contains the path to a SRA file
associated with the sample.  Since KIQ uses the NCBI SRA API, the second column
can also contain just the SRA identifier (starting with `SRR`).  In this case,
the program first checks if the SRA file is already in the local cache and then
reads it from there, or otherwise it will be downloaded on the fly and stored
in the local cache.  The cache directory can be set using the
[vdb-config](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) tool
from the [sra-tools package](https://github.com/ncbi/sra-tools/).

The option `-z` specifies the number threads that are used for k-mer counting.

Additional datasets can be added to an existing database by using the option `-a`.


### Query database by a k-mer

After counting the k-mers, the database can be queried
with `kiq query`, either by providing a list of query k-mers in a text file
using option `-q` or by providing query k-mer(s) using option `-Q`.

For example:
```
kiq query -i kmer_index.bin -k kiq_database.bin -q query.txt

kiq query -i kmer_index.bin -k kiq_database.bin -Q ACGTACGTACGTACGTACGTACGTACGTACGT

kiq query -i kmer_index.bin -k kiq_database.bin -Q ACGTACGTACGTACGTACGTACGTACGTACGT,GCGTACGTACGTACGTACGTACGTACGTACGG
```

### Export k-mer database

The k-mer database and metadata can be exported using `kiq dump`:
```
kiq dump -i kmer_index.bin -k kiq_database.bin -p stats

kiq dump -i kmer_index.bin -k kiq_database.bin -p metadata

kiq dump -i kmer_index.bin -k kiq_database.bin -p long

kiq dump -i kmer_index.bin -k kiq_database.bin -p stats
```

### Modify database
KIQ's k-mer database can be modified using `kiq modify`, which reads a
tab-separated file containing instructions for modifying the database. The
first column contains the experiment ID, the second column contains a command
and the third column contains an optional argument.  Two different commands are
available: `delete` for deleting an experiment (and its associated k-mer
counts) from the database, and `update_desc` for updating the description for a
given experiment ID.

For example:
```
SRRZZZ	delete
SRRXXX	update_desc	Total RNA-Seq from Drosophila heads.
```
Run `kiq modify`:
```
kiq modify -i kmer_index.bin -k kiq_database.bin -c modifications.tsv
```


### Acknowledgments

KIQ uses the [BBHash](https://github.com/rizkg/BBHash) library for indexing of the k-mers
and [zstr](https://github.com/mateidavid/zstr) library for reading compressed FASTA/Q files.


### License
See the file LICENSE.

