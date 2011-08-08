#include <iostream>
#include <SuffixTree.h>
#include <Array.h>
#include "Bits.h"

#ifdef _cplusplus
#define _STDC_CONSTANT_MACROS
#ifdef _STDINT_H
#undef _STDINT_H
#endif
#include <cstdint>
#endif

class FactorWriter
{
    public:

        /** Destructor for the class. */
        virtual ~FactorWriter() {}

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        virtual void output_factor(uint64_t pos, uint64_t len) = 0;
};

class FactorWriterText : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         */
        FactorWriterText(ofstream& outfile);
        
        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void output_factor(uint64_t pos, uint64_t len);

    private:
        
        // Output stream to write to
        ofstream& outfile;
};

class FactorWriterBinary : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param maxposbits Number of bits for encoding positions
         */
        FactorWriterBinary(ofstream& outfile, uint64_t maxposbits);

        /** Destructor for the class. */
        ~FactorWriterBinary();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void output_factor(uint64_t pos, uint64_t len);

    private:

        // To write bits and integers
        BitWriter *bwriter;

        // To Golomb encode numbers
        GolombCoder *gcoder;

        // Maximum number of bits to use to encode a position
        uint64_t maxposbits;
};

// A base class for RLZ compression and decompression
class RLZ
{
    public:
        static const uint64_t BITSPERBASE = 3;
        static const char *NUCLALPHA;
        static const uint64_t NUCLALPHASIZE = 5;

    protected:
        
        // Reference sequence as a bit vector with 3bpb encoding
        // {a,c,g,t,n}
        cds_utils::Array *refseq;
        uint64_t refseqlen;
        uint64_t logrefseqlen;

        // File names of sequences to be compressed or decompressed
        char **filenames;
        uint64_t numfiles;

        /** Store a sequence containing nucleotides from alphabet
         * NUCLALPHA using BITSPERBASE bits each.
         * @param sequence Character array of the sequence
         * @param filename Name of input file
         * @param dest Place to store the sequence to
         * @param length Number of symbols to store
         */
        virtual void store_sequence(char *sequence, char *filename,
                                    Array *dest, uint64_t length);

        /** Store a sequence containing nucleotides from alphabet
         * NUCLALPHA using BITSPERBASE bits each.
         * @param infile Input stream containing sequence
         * @param filename Name of input file
         * @param dest Place to store the sequence to
         * @param length Number of symbols to store
         */
        virtual void store_sequence(ifstream& infile, char *filename,
                                    Array *dest, uint64_t length);

};

class RLZCompress : RLZ
{
    public:

        /** Constructor for the RLZ compress class.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         * @param encoding Type of encoding to be used
         */
        RLZCompress(char **filenames, uint64_t numfiles, char encoding='b');

        /** Temporary constructor that implements the suffix tree
         * instead of a suffix array.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         * @param state Random parameter to overload the constructor
         */
        RLZCompress(char **filenames, uint64_t numfiles, bool state);

        /** Destructor for the class. */
        ~RLZCompress();

        /** Method for compression. */
        void compress();

    private:

        // Suffix tree of the reference sequence
        cds_static::SuffixTree *st;

        // Suffix array of the reference sequence
        cds_utils::Array *sa;

        // Type of encoding
        char encoding;

        /** Conducts the relative Lempel-Ziv compression of the sequence
         * inside the infile and writes the output to outfile.
         * @param infile Input file stream
         * @param filename Name of the input file
         * @param facwriter FactorWriter object to output factors
         */
        void relative_LZ_factorise(ifstream& infile, char *filename,
                                   FactorWriter& facwriter);

        /** Conducts the relative Lempel-Ziv compression of the sequence
         * inside the infile and writes the output to outfile.
         * @param infile Input file stream
         * @param filename Name of the input file
         * @param outfile Output file stream
         * @param state random parameter to overload the method
         */
        void relative_LZ_factorise(ifstream& infile, char *filename,
                                   ofstream& outfile, bool state);


        /** Conducts a binary search in the suffix array for a symbol at
         * a particular offset of the suffixes.
         * @param pl Left boundary to begin with
         * @param pr Right boundary to begin with
         * @param c Symbol to find at offset
         * @param offset Offset from the beginning of the suffix
         * @param cl Output left boundary
         * @param cr Output right boundary
         */
        void sa_binary_search(uint64_t pl, uint64_t pr, int c, 
                              uint64_t offset, uint64_t *cl, uint64_t *cr);

};

class RLZDecompress
{

    public:

        /** Constructor for the RLZ decompressor class.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         * @param encoding Type of encoding to be used
         */
        RLZDecompress(char **filenames, uint64_t numfiles);

        /** Destructor for the class. */
        ~RLZDecompress();

        /** Method for decompression. */
        void decompress();

    private:
};
