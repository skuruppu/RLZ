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
         */
        RLZ(char **filenames);

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

        // File names of sequences to be compressed
        char **filenames;

        /** Store a sequence containing nucleotides from alphabet
         * NUCLALPHA using BITSPERBASE bits each.
         * @param infile Input file
         * @param filename Name of input file
         * @param sequence Place to store the sequence to
         * @param length Number of symbols to store
         */
        void store_sequence(ifstream& infile, char *filename,
                            Array *sequence, uint64_t length);
};
