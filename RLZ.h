/* RLZ compress
 * Copyright (C) 2011 Shanika Kuruppu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * RLZ - Relative Lempel Ziv
 * Implements the RLZ compression algorithm.
 * Authors: Shanika Kuruppu (kuruppu@csse.unimelb.edu.au)
 *          Simon J. Puglisi (simon.puglisi@rmit.edu.au)
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <SuffixTree.h>

#include "Bits.h"
#include "lib_wrapper/wrapper.h"

#ifdef _cplusplus
#define _STDC_CONSTANT_MACROS
#ifdef _STDINT_H
#undef _STDINT_H
#endif
#include <cstdint>
#endif

#ifdef CDS
    #define Array lib_wrapper::CDSArray
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
         * @param isliss Whether to LISS encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         * @param logrefseqlen Number of bits for encoding positions
         */
        FactorWriter(ofstream& outfile, char encoding, bool isshort,
                     bool isliss, Array *refseq, uint64_t refseqlen,
                     uint64_t logrefseqlen);

        /** Destructor for the class. */
        virtual ~FactorWriter();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        virtual void write_factor(uint64_t pos, uint64_t len);

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         * @param lissfac Factor to be encoded is part of the LISS
         */
        virtual void write_factor(uint64_t pos, uint64_t len, 
                                  bool lissfac) {}

        /** Finalise any writing that hasn't completed yet. */
        virtual void end_of_sequence();

    protected:

        // Constants for Golomb coding
        static const unsigned int GOLOMBDIV = 64;
        static const unsigned int GOLOMBDIVSHORT = 8;

        // 2*len+len/GOLOMBDIVSHORT+(LOG2GOLOMBDIVSHORT+1) <
        // logrefseqlen+len/GOLOMBDIV+(LOG2GOLOMBDIV+1)
        uint64_t SHORTFACTHRESH;
        
        // Reference sequence as a bit vector with 3bpb encoding
        // {a,c,g,t,n}
        Array *refseq;
        uint64_t refseqlen;
        uint64_t logrefseqlen;

        // Whether to short factor encode or not
        bool isshort;

        // LISS factor encoding
        bool isliss;
        std::vector<uint64_t> positions;
        std::vector<uint64_t> lengths;

    private:
        
        FactorWriter *facwriter;

        /** Finds longest strictly increasing subsequence (LISS).
         *  O(n log k) algorithm, k is the length of the LISS.
         *  @param a Input set of integers
         *  @param b Indexes of integers in a that belong to the LISS
         */
        void find_LISS(std::vector<uint64_t>& a,
                       std::vector<uint64_t>& b);
};

class FactorWriterText : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param isshort Whether to short factor encode or not
         * @param isliss Whether to LISS encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         */
        FactorWriterText(ofstream& outfile, bool isshort, bool isliss,
                         Array *refseq, uint64_t refseqlen,
                         uint64_t logrefseqlen);
        
        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void write_factor(uint64_t pos, uint64_t len);

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         * @param lissfac Factor to be encoded is part of the LISS
         */
        void write_factor(uint64_t pos, uint64_t len, bool lissfac);

        /** Finalise any writing that hasn't completed yet. */
        void end_of_sequence();

    private:
        
        // Output stream to write to
        ofstream& outfile;
        
        // Keeps track of if the first LISS factor has been encoded or
        // not
        bool firstliss;
};

class FactorWriterBinary : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param isshort Whether to short factor encode or not
         * @param isliss Whether to LISS encode or not
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         * @param logrefseqlen Number of bits for encoding positions
         */
        FactorWriterBinary(ofstream& outfile, bool isshort, bool isliss,
                           Array *refseq, uint64_t refseqlen,
                           uint64_t logrefseqlen);

        /** Destructor for the class. */
        ~FactorWriterBinary();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void write_factor(uint64_t pos, uint64_t len);

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         * @param lissfac Factor to be encoded is part of the LISS
         */
        void write_factor(uint64_t pos, uint64_t len, bool lissfac);

        /** Finalise any writing that hasn't completed yet. */
        void end_of_sequence();

    private:

        // To write bits and integers
        BitWriter *bwriter;

        // To Golomb encode numbers
        GolombCoder *gcoder;

        // To Golomb encode short numbers
        GolombCoder *gcodershort;

        // Keeps track of if the first LISS factor has been encoded or
        // not
        bool firstliss;
};

class FactorWriterIndex : public FactorWriter
{
    public:

        /** Constructor for the class.
         * @param outfile Output file stream
         * @param refseq Reference sequence
         * @param refseqlen Length of reference sequence
         * @param logrefseqlen Number of bits for encoding positions
         * @param displayonly Output only display() query structures
         */
        FactorWriterIndex(ofstream& outfile, Array *refseq,
                          Array *sa, uint64_t refseqlen,
                          uint64_t logrefseqlen, bool displayonly);


        /** Destructor for the class. */
        ~FactorWriterIndex();

        /** Output an RLZ factor.
         * @param pos Position component of factor
         * @param len Length component of factor
         */
        void write_factor(uint64_t pos, uint64_t len);

        /** Finalise any writing that hasn't completed yet. */
        void end_of_sequence();

        /** Write the index out to disk. */
        void write_index();

    private:

        // Output stream to write to
        ofstream& outfile;

        // To write bits and integers
        BitWriter *bwriter;

        // Factor start positions
        std::vector<bool> facstarts;

        // Cumulative sequence lengths
        std::vector<uint64_t> cumseqlens;

        // List of factor positions
        unsigned int *positions;
        size_t posarraylen;

        // Total number of factors
        uint64_t numfacs;

        // Accumulate the length of a sequence
        uint64_t cumlen;

        // If this is true, only write the data structures required to
        // implement display() query
        bool displayonly;

        // Suffix array of the reference sequence
        Array *sa;

        // Nested level lists of the sorted factors
        Array *nll;
        Array *levelidx;
        uint32_t numlevels;

        // Sequence start positions in factors
        std::vector<bool> seqstarts;
        uint32_t currseqfacnum;

        // Bit vectors to store which positions factors start and end at
        cds_utils::BitString *isstart;
        cds_utils::BitString *isend;

        void construct_nested_level_list(cds_static::BitSequenceSDArray&
                                         facstartssdarray);
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
         * @param isliss Whether LISS encode was used
         */
        FactorReaderText(ifstream& infile, bool isshort, bool isliss);
        
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

        // LISS factor encoding
        bool isliss;

        // Variables need for LISS encoding
        bool firstliss;
        uint64_t prevpos, cumlen;
};

class FactorReaderBinary : public FactorReader
{
    public:

        /** Constructor for the class.
         * @param infile Input file stream
         * @param logrefseqlen Number of bits for encoding positions
         * @param isshort Whether some factors will be short factors
         * @param isliss Whether LISS encode was used
         */
        FactorReaderBinary(ifstream& infile, uint64_t logrefseqlen, 
                           bool isshort, bool isliss);

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

        // LISS factor encoding
        bool isliss;

        // Maximum number of bits to use to encode a position
        uint64_t logrefseqlen;

        // Variables need for LISS encoding
        bool firstliss;
        uint64_t prevpos, cumlen;
};

// A base class for RLZ compression and decompression
class RLZ
{
    protected:
        
        // Reference sequence as a bit vector with 3bpb encoding
        // {a,c,g,t,n}
        Array *refseq;
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
         * @param isliss Enable LISS factor encoding
         */
        RLZCompress(char **filenames, uint64_t numfiles, 
                    char encoding='b', bool isshort=false, 
                    bool isliss=false);

        /** Constructor for the RLZ compress class.
         * @param filenames Filenames for sequences to be compressed
         * @param numfiles Number of files in the dataset
         * @param idxname Name to give to the index
         * @param displayonly Only implement display() query
         */
        RLZCompress(char **filenames, uint64_t numfiles, char *idxname,
                    bool displayonly);

        /** Destructor for the class. */
        ~RLZCompress();

        /** Method for compression. */
        void compress();

    private:

        // Suffix array of the reference sequence
        Array *sa;

        // Type of encoding
        char encoding;

        // Short factor encoding
        bool isshort;

        // LISS encoding
        bool isliss;

        // Name of the output index file
        char *idxname;
        
        // If this is true, only write the data structures required to
        // implement display() query
        bool displayonly;

        /** Read the reference sequence and construct the suffix array
         */
        void read_refseq_and_construct_sa();

        /** Read the reference sequence and construct the suffix array
         */
        void read_refseq_and_sa();

        /** Conducts the relative Lempel-Ziv compression of the sequence
         * inside the infile and writes the output to outfile.
         * @param infile Input file stream
         * @param filename Name of the input file
         * @param facwriter FactorWriter object to output factors
         */
        void relative_LZ_factorise(ifstream& infile, char *filename,
                                   FactorWriter& facwriter);

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
