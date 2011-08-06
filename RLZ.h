#include <iostream>
#include <SuffixTree.h>
#include <Array.h>

#ifdef _cplusplus
#define _STDC_CONSTANT_MACROS
#ifdef _STDINT_H
#undef _STDINT_H
#endif
#include <cstdint>
#endif

class RLZ
{
    public:
        static const uint64_t BITSPERBASE = 3;
        static const char *NUCLALPHA;
        static const uint64_t NUCLALPHASIZE = 5;

        /** Constructor for the RLZ class.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         */
        RLZ(char **filenames, uint64_t numfiles);

        /** Temporary constructor that implements the suffix tree
         * instead of a suffix array.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         * @param state Random parameter to overload the constructor
         */
        RLZ(char **filenames, uint64_t numfiles, bool state);

        ~RLZ();

        void compress();

        void decompress();

    private:
        // Reference sequence as a bit vector with 3bpb encoding
        // {a,c,g,t,n}
        cds_utils::Array *refseq;
        uint64_t refseqlen;
        uint64_t logrefseqlen;

        // Suffix tree of the reference sequence
        cds_static::SuffixTree *st;

        // Suffix array of the reference sequence
        cds_utils::Array *sa;

        // File names of sequences to be compressed
        char **filenames;
        uint64_t numfiles;

        /** Store a sequence containing nucleotides from alphabet
         * NUCLALPHA using BITSPERBASE bits each.
         * @param sequence Character array of the sequence
         * @param filename Name of input file
         * @param dest Place to store the sequence to
         * @param length Number of symbols to store
         */
        void store_sequence(char *sequence, char *filename,
                            Array *dest, uint64_t length);

        /** Store a sequence containing nucleotides from alphabet
         * NUCLALPHA using BITSPERBASE bits each.
         * @param infile Input stream containing sequence
         * @param filename Name of input file
         * @param dest Place to store the sequence to
         * @param length Number of symbols to store
         */
        void store_sequence(ifstream& infile, char *filename,
                            Array *dest, uint64_t length);

        /** Conducts the relative Lempel-Ziv compression of the sequence
         * inside the infile and writes the output to outfile.
         * @param infile Input file stream
         * @param filename Name of the input file
         * @param outfile Output file stream
         */
        void relative_LZ_factorise(ifstream& infile, char *filename,
                                   ofstream& outfile);

        /** Conducts the relative Lempel-Ziv compression of the sequence
         * inside the infile and writes the output to outfile.
         * @param infile Input file stream
         * @param filename Name of the input file
         * @param outfile Output file stream
         * @param random parameter to overload the method
         */
        void relative_LZ_factorise(ifstream& infile, char *filename,
                                   ofstream& outfile, bool state);
};
