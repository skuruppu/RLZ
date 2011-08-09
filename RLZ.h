#include <iostream>
#include <vector>
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

        /** Constructor for the class. */
        FactorWriter();

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param encoding Type of encoding to be used
         * @param isshort Whether to short factor encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         * @param logrefseqlen Number of bits for encoding positions
         */
        FactorWriter(ofstream& outfile, char encoding, bool isshort,
                     cds_utils::Array *refseq, uint64_t refseqlen,
                     uint64_t logrefseqlen);

        /** Destructor for the class. */
        virtual ~FactorWriter();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        virtual void write_factor(uint64_t pos, uint64_t len);

    protected:

        // Constants for Golomb coding
        static const unsigned int GOLOMBDIV = 64;
        static const unsigned int GOLOMBDIVSHORT = 8;

        // 2*len+len/GOLOMBDIVSHORT+(LOG2GOLOMBDIVSHORT+1) <
        // logrefseqlen+len/GOLOMBDIV+(LOG2GOLOMBDIV+1)
        uint64_t SHORTFACTHRESH;
        
        // Reference sequence as a bit vector with 3bpb encoding
        // {a,c,g,t,n}
        cds_utils::Array *refseq;
        uint64_t refseqlen;
        uint64_t logrefseqlen;

    private:
        
        FactorWriter *facwriter;
};

class FactorWriterText : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param isshort Whether to short factor encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         */
        FactorWriterText(ofstream& outfile, bool isshort,
                         cds_utils::Array *refseq, uint64_t refseqlen,
                         uint64_t logrefseqlen);
        
        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void write_factor(uint64_t pos, uint64_t len);

    private:
        
        // Output stream to write to
        ofstream& outfile;

        // Whether to short factor encode or not
        bool isshort;
};

class FactorWriterBinary : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param isshort Whether to short factor encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         * @param logrefseqlen Number of bits for encoding positions
         */
        FactorWriterBinary(ofstream& outfile, bool isshort,
                           cds_utils::Array *refseq, uint64_t refseqlen,
                           uint64_t logrefseqlen);

        /** Destructor for the class. */
        ~FactorWriterBinary();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void write_factor(uint64_t pos, uint64_t len);

    private:

        // To write bits and integers
        BitWriter *bwriter;

        // To Golomb encode numbers
        GolombCoder *gcoder;

        // To Golomb encode short numbers
        GolombCoder *gcodershort;

        // Whether to short factor encode or not
        bool isshort;
};

class FactorReader
{
    public:

        FactorReader ();

        /** Constructor for the class.
         * @param infile Input file stream
         * @param logrefseqlen Number of bits for encoding positions
         */
        FactorReader(ifstream& infile, uint64_t logrefseqlen);

        /** Destructor for the class. */
        virtual ~FactorReader();

        /** Read an RLZ factor.
         * @param pos Output position component of factor
         * @param len Output length component of factor
         * @param substr Vector to store short factor if needed
         * @return Status to indicate if a factor was read successfully
         */
        virtual bool read_factor(uint64_t *pos, uint64_t *len, 
                                 std::vector<char>& substr);

    private:
        
        FactorReader *facreader;
};

class FactorReaderText : public FactorReader
{
    public:

        /** Constructor for the class.
         * @param infile Input file stream
         * @param isshort Whether some factors will be short factors
         */
        FactorReaderText(ifstream& infile, bool isshort);
        
        /** Read an RLZ factor.
         * @param pos Output position component of factor
         * @param len Output length component of factor
         * @param substr Vector to store short factor if needed
         * @return Status to indicate if a factor was read successfully
         */
        bool read_factor(uint64_t *pos, uint64_t *len,
                         std::vector<char>& substr);

    private:
        
        // Input stream to read from
        ifstream& infile;

        // Whether some factors will be short factors
        bool isshort;
};

class FactorReaderBinary : public FactorReader
{
    public:

        /** Constructor for the class.
         * @param infile Input file stream
         * @param logrefseqlen Number of bits for encoding positions
         * @param isshort Whether some factors will be short factors
         */
        FactorReaderBinary(ifstream& infile, uint64_t logrefseqlen, 
                           bool isshort);

        /** Destructor for the class. */
        ~FactorReaderBinary();

        /** Read an RLZ factor.
         * @param pos Output position component of factor
         * @param len Output length component of factor
         * @param substr Vector to store short factor if needed
         * @return Status to indicate if a factor was read successfully
         */
        bool read_factor(uint64_t *pos, uint64_t *len,
                         std::vector<char>& substr);

    private:

        // To read bits and integers
        BitReader *breader;

        // To Golomb decode numbers
        GolombCoder *gdecoder;

        // To Golomb decode short numbers
        GolombCoder *gdecodershort;

        // Whether some factors will be short factors
        bool isshort;

        // Maximum number of bits to use to encode a position
        uint64_t logrefseqlen;
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
         * @param isshort Encode shorter factors as substr,len pairs
         */
        RLZCompress(char **filenames, uint64_t numfiles, 
                    char encoding='b', bool isshort=false);

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

        // Short factor encoding
        bool isshort;

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

class RLZDecompress : RLZ
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

        /** Conducts the relative Lempel-Ziv decompression of the
         * sequence inside the infile and writes the output to outfile.
         * @param facreader FactorReader object to read factors from
         * @param filename Name of the input file
         * @param outfile Output file stream
         */
        void relative_LZ_defactorise(FactorReader& facreader,
                                     char *filename, ofstream& outfile);
};
