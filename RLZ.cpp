#include "RLZ.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

const char *RLZ::NUCLALPHA = "acgnt";

unsigned char nucl_to_int[256] = {0};
unsigned char int_to_nucl[RLZ::NUCLALPHASIZE+1];

void initialise_nucl_converters()
{
    uint64_t i;

    for (i=0; i<RLZ::NUCLALPHASIZE; i++)
    {
        nucl_to_int[(unsigned int)RLZ::NUCLALPHA[i]] = i+1;
        nucl_to_int[toupper(RLZ::NUCLALPHA[i])] = i+1;
    }

    for (i=0; i<RLZ::NUCLALPHASIZE; i++)
    {
        int_to_nucl[i+1] = RLZ::NUCLALPHA[i];
    }
}

RLZ::RLZ(char **filenames)
{
    this->filenames = filenames;

    // Open reference sequence file
    ifstream infile;
	infile.open(filenames[0], ifstream::in);
    if (infile.bad())
    {
        cerr << "Could not open file " << filenames[0] << ".\n";
        exit(1);
    }

    // Get the length of the file
    infile.seekg(0, ios::end);
    refseqlen = infile.tellg();
    infile.seekg(0, ios::beg);

    // Get the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = (((unsigned)1<<i) != refseqlen) ? i+1 : i;

    initialise_nucl_converters();

    // Read the reference sequence
    refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
    store_sequence(infile, filenames[0], refseq, refseqlen);

}

RLZ::~RLZ()
{
    delete refseq;
    //delete st;
}

void RLZ::store_sequence(ifstream& infile, char *filename, 
                         Array *sequence, uint64_t length)
{
    uint64_t i;
    int c;
    unsigned int v;

    for (i=0; i<length; i++)
    {
        c = infile.get();
        v = nucl_to_int[c];
        // Valid nucleotide
        if (v > 0)
        {
            refseq->setField(i, v);
        }
        else
        {
            cerr << "Invalid symbol " << c << " at position " << i;
            cerr << " of sequence in " << filename << ".\n";
            exit(1);
        }
    }

    return;
}
