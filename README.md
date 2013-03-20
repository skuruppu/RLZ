Description
===========

RLZ is a program that compresses sets of very similar DNA sequences where one
of the sequences in the dataset acts as a reference sequence, and the
remaining sequences are compressed by an LZ77 parsing with respect to the
chosen reference sequence.

This package includes the source code for compression and decompression, as
well as the RLZ self-index.

Distribution
============

Refer to the `COPYING` file in the directory for licensing information of
this program.

Citations
=========

If this code is used, please cite the following papers and thesis:

S. Kuruppu, S. J. Puglisi and J. Zobel,
[*Relative Lempel-Ziv Compression of Genomes for Large-Scale Storage and
Retrieval*](http://dx.doi.org/10.1007/978-3-642-16321-0_20).
Proc. 17th International Symposium on String Processing and Information
Retrieval (SPIRE 2010)
Lecture Notes in Computer Science, Volume 6393, (2010) pp. 201-206.

S. Kuruppu, S. J. Puglisi and J. Zobel,
[*Optimized Relative Lempel-Ziv Compression of
Genomes*](http://crpit.com/confpapers/CRPITV113Kuruppu.pdf).
Proc. 34th Australian Computer Science Conference
(ACSC2011) CRPIT, Volume 11, (2011) pp. 91-98.

S. Kuruppu,
[*Compression of Large DNA
Databases*](http://dtl.unimelb.edu.au//exlibris/dtl/d3_1/apache_media/L2V4bGlicmlzL2R0bC9kM18xL2FwYWNoZV9tZWRpYS8yODQ3Nzg=.pdf).
Ph.D. Thesis, 2012
The University of Melbourne.

Requirements
============

1. Version 1.0.8 of the [libcds library](https://github.com/fclaude/libcds).
   The software is not tested for earlier versions of libcds.

2. [libdivsufsort](https://code.google.com/p/libdivsufsort)

Install the above libraries and add the paths to the include files to the
`CPLUS_INCLUDE_PATH` environment variable, and the paths to the library files
(`.so` and `.a` files) to `LIBRARY_PATH` and `LD_LIBRARY_PATH` environment
variables.

Building
========

Run `make`.

Usage
=====

Sequence files
--------------

The sequence files must only contain the nucleotides '`a`', '`c`', '`g`',
'`t`' and '`n`' (in lower case). If the sequence contains any other
nucleotides, it's best to replace them with the '`n`' nucleotides prior to
using RLZ. The program assumes that the sequence files being read uses a byte
per symbol. In other words, if the sequence is of length n, then the sequence
file should be n bytes. The file should not have any terminating symbols such
as '`\0`' or '`\n`'.

RLZ
---

    Usage: rlz [OPTIONS] REF FILE1 FILE2 ...
        -d: Decompress (all other options ignored)
        -e: Type of encoding (t: text, b: binary) (default: b)
        -i: Output a self-index with given name (all options except -r ignored)
        -l: Enable LISS encoding
        -r: Only enable random access in the index (used with -i)
        -s: Output short factors as substring and length pairs
        REF: Name of reference sequence
        FILE1 ...: Names of files to be compressed

### Examples ###

Basic RLZ compression can be done by invoking the following command:

    ./rlz <ref seq file name> <seq1 file name> <seq2 file name> ...

This produces output files with the same name as the input sequence files and
a '`.fac`' suffix. For example, '`<seq1 file name>.fac`'.

RLZ decompression can be done by invoking the following command:

    ./rlz -d <ref seq file name> <seq1 file name> <seq2 file name> ...

The input file names to the decompression command should not have the '`.fac`'
suffix and should be the same as the input to the compression invocation.

An RLZ self-index can be produced with the following command:

    ./rlz -i <index file name> <ref seq file name> <seq1 file name> <seq2 file name> ...

Calling with the `-r` option will produce an index that just supports the
`display()` query.

RLZ Index
---------

    Usage: rlz_index MODE FILE
         MODE: Mode for using the index (d: display, c: count, l: locate)
         FILE: Compressed index file

### Examples ###

To use the RLZ self-index for `display()` queries, invoke the following
command:

    ./rlzindex d <index file name>

The program will continue to prompt for query input until `Ctrl+D` is pressed.
A query should be in one line in the following format:

    <sequence number> <start position> <end position>

Sequences are numbered starting from 0, where sequence 0 is the reference and
the remaining sequences are ordered in the same order as the arguments were
passed in to the '`rlz`' invocation.

Positions also start from 0.

Then end position is one greater than the last position to retrieve.
Therefore, (<end position> - <start position>) = length of subsequence to
retrieve.

To use the RLZ self-index for `count()` queries, invoke the following command:

    ./rlzindex c <index file name>

The program will continue to prompt for query input until `Ctrl+D` is pressed.
A query should be in one line in the following format:

    <pattern to search>

To use the RLZ self-index for `locate()` queries, invoke the following
command:

    ./rlzindex l <index file name>

The program will continue to prompt for query input until `Ctrl+D` is pressed.
A query should be in one line in the following format:

    <pattern to search>

Tests
=====

Please refer to the runTests.sh script for further examples on how to run the
`rlz` and `rlzindex` programs. Test input are also provided in the test
directory.

To run the tests, invoke the following command:

    ./runTests.sh

A script, `queryTestGenerator.py`, is also included for your convenience to
generate test cases for queries.
