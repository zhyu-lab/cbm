# CBM

CBM is a tool for detecting tumor subclones from single-cell mutation data.

## Requirements

* Linux systems.
* CMake3.0+.
* g++.

## Installation

To build binary, do as follows:

```
tar -zxvf CBM.tar.gz
cd CBM
cmake .
make
```

After the installation, the main program of CBM is generated in “bin” directory. Type following command if you want to add CBM to system path:
```
make install
```

## Usage

CBM uses single-cell binary genotype matrix to infer tumor subclones and their genotypes.

Example:

```
cbm -i testdata/example.txt -o testdata/example
```

## Input Files

### Genotype Matrix

The SNVs of single cells are denoted as a genotype matrix. Each row defines the mutation profiles of a single cell, and each column represents a mutation. Columns are separated by tabs. The genotype matrix is binary.

The entry at position [i,j] should be

* 0 if mutation j is not observed in cell i,
* 1 if mutation j is observed in cell i, or
* 3 if the genotype information is missing

## Output Files

The base name of the output files is provided by users.

### Genotypes of subclones

The genotypes of subclones are written to a file with suffix "clones.txt".

### Cell assignments

The cell-to-subclone assignments are written to a file with suffix "assignment.txt".

## Arguments

* `-i, --input <filename>` Replace \<filename\> with the file containing the genotype matrix.

* `-o, --output <string>` Replace \<string\> with the base name of the output file.

## Optional arguments

* `-K, --maxc <INT>` Set \<INT\> to a positive integer. This specifies the maximum number of subclones to consider.

* `-a, --alpha <Double>` Set \<Double\> to estimated false positive rate of the single-cell sequencing experiment.

* `-b, --beta <Double>` Set \<Double\> to estimated false negative rate of the single-cell sequencing experiment.

## Contact

If you have any questions, please contact zhyu@nxu.edu.cn.