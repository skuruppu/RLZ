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

        RLZ(std::ifstream& infile);

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

        void store_sequence(ifstream& infile, Array *sequence, uint64_t length);
};
