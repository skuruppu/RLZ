#include "RLZ.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

unsigned char nucl_to_int[256] = {0};
unsigned char int_to_nucl[6] = {0};

void initialise_nucl_converters()
{
    nucl_to_int['a'] = 1;
    nucl_to_int['A'] = 1;
    nucl_to_int['c'] = 2;
    nucl_to_int['C'] = 2;
    nucl_to_int['g'] = 3;
    nucl_to_int['G'] = 3;
    nucl_to_int['n'] = 4;
    nucl_to_int['N'] = 4;
    nucl_to_int['t'] = 5;
    nucl_to_int['T'] = 5;

    int_to_nucl[1] = 'a';
    int_to_nucl[2] = 'c';
    int_to_nucl[3] = 'g';
    int_to_nucl[4] = 'n';
    int_to_nucl[5] = 't';
}

RLZ::RLZ(ifstream& infile)
{
    // Get the length of the file
    infile.seekg(0, ios::end);
    refseqlen = infile.tellg();
    infile.seekg(0, ios::beg);

    // Get the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((1<<i) != refseqlen) ? i+1 : i;

    initialise_nucl_converters();

    // Read the reference sequence
    refseq = new Array(refseqlen, (1<<BITSPERBASE)-1);
    store_sequence(infile, refseq, refseqlen);

}

RLZ::~RLZ()
{

}

void RLZ::store_sequence(ifstream& infile, Array *sequence,
                         uint64_t length)
{
    uint64_t i;
    int c;
    unsigned int v;

    for (i=0; i<length; i++)
    {
        c = infile.get();
        v = nucl_to_int[c];
        if (v > 0)
        {
            refseq->setField(i, v);
        }
        else
        {
            cerr << "Invalid symbol " << c << " at position " << i 
                 << ".\n";
            exit(1);
        }
    }

    return;
}
