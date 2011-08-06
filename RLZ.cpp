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

    char *sequence = NULL;
    if (loadText(filenames[0], &sequence, &refseqlen))
    {
        cerr << "Couldn't read reference sequence.\n";
        exit(1);
    }
    // loadText places an extra byte at the end
    refseqlen--;

    // Get the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = (((unsigned)1<<i) != refseqlen) ? i+1 : i;

    initialise_nucl_converters();

    // Read the reference sequence
    refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
    store_sequence(sequence, filenames[0], refseq, refseqlen);

    char stfilename[1024];
    sprintf(stfilename, "%s.st", filenames[0]);
    ifstream infile;
    infile.open(stfilename, ifstream::in);
    // Construct the suffix tree if it doesn't exist
    if (!infile.good())
    {
        infile.close();
        SuffixTreeY sty(sequence, refseqlen+1, DAC, CN_NPR, 32);
        ofstream outfile(stfilename);
        sty.save(outfile);
        outfile.close();
        infile.open(stfilename, ifstream::in);
    }
    // Load from saved suffix tree file
    st = SuffixTree::load(infile);
    infile.close();

    delete [] sequence;
}

RLZ::~RLZ()
{
    delete refseq;
    delete st;
}

void RLZ::store_sequence(char *sequence, char *filename,
                         Array *dest, uint64_t length)
{
    uint64_t i;
    unsigned int v;

    for (i=0; i<length; i++)
    {
        v = nucl_to_int[(unsigned int)sequence[i]];
        // Valid nucleotide
        if (v > 0)
        {
            dest->setField(i, v);
        }
        else
        {
            cerr << "Invalid symbol " << sequence[i] << " at position " << i;
            cerr << " of sequence in " << filename << ".\n";
            exit(1);
        }
    }

    return;
}
