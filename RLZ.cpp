#include <divsufsort64.h>
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


RLZ::RLZ(char **filenames, uint64_t numfiles)
{
    this->filenames = filenames;
    this->numfiles = numfiles;

    initialise_nucl_converters();

    uint64_t i;
    char safilename[1024];
    sprintf(safilename, "%s.sa", filenames[0]);
    ifstream infile;
    infile.open(safilename, ifstream::in);
    // Need to construct suffix tree
    if (!infile.good())
    {
        infile.close();
        // Read reference sequence into memory since its needed by
        // suffix tree constructor
        char *sequence = NULL;
        if (loadText(filenames[0], &sequence, &refseqlen))
        {
            cerr << "Couldn't read reference sequence.\n";
            exit(1);
        }
        // loadText places an extra byte at the end
        refseqlen--;

        // Read the reference sequence
        refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(sequence, filenames[0], refseq, refseqlen);

        // Construct suffix array
        uint64_t *sufarray = new uint64_t[refseqlen+1];
        if (divsufsort64((sauchar_t*)sequence, (saidx64_t*)sufarray,
            refseqlen+1) < 0)
        {
            cerr << "Error in constructing suffix array.\n";
            exit(1);
        }

        sa = new Array(refseqlen+1, logrefseqlen);
        for (i=0; i<=refseqlen; i++)
        {
            sa->setField(i, sufarray[i]);
        }

        // Write out suffix array to disk for later use
        ofstream outfile(safilename);
        sa->save(outfile);
        outfile.close();

        delete [] sequence;
        delete [] sufarray;
    }
    else
    {
        // Load suffix array from saved suffix array file
        sa = new Array(infile);
        infile.close();
        
        // Open reference sequence file
        infile.open(filenames[0], ifstream::in);
        if (!infile.good())
        {
            cerr << "Couldn't open file " << filenames[0] << ".\n";
            exit(1);
        }

        // Get the reference sequence length
        infile.seekg(0, ios::end);
        refseqlen = infile.tellg();
        infile.seekg(0, ios::beg);

        // Read the reference sequence
        refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(infile, filenames[0], refseq, refseqlen);
        infile.close();
    }
}

RLZ::RLZ(char **filenames, uint64_t numfiles, bool state)
{
    this->filenames = filenames;
    this->numfiles = numfiles;

    initialise_nucl_converters();

    char stfilename[1024];
    sprintf(stfilename, "%s.st", filenames[0]);
    ifstream infile;
    infile.open(stfilename, ifstream::in);
    // Need to construct suffix tree
    if (!infile.good())
    {
        infile.close();
        // Read reference sequence into memory since its needed by
        // suffix tree constructor
        char *sequence = NULL;
        if (loadText(filenames[0], &sequence, &refseqlen))
        {
            cerr << "Couldn't read reference sequence.\n";
            exit(1);
        }
        // loadText places an extra byte at the end
        refseqlen--;

        // Read the reference sequence
        refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(sequence, filenames[0], refseq, refseqlen);

        // Construct suffix tree
        SuffixTreeY *sty = new SuffixTreeY(sequence, refseqlen+1, PHI,
                                           CN_NPR, 32);

        // Write out suffix tree to disk for later use
        ofstream outfile(stfilename);
        sty->save(outfile);
        outfile.close();

        st = sty;

        delete [] sequence;
    }
    else
    {
        // Load suffix tree from saved suffix tree file
        st = SuffixTree::load(infile);
        infile.close();
        
        // Open reference sequence file
        infile.open(filenames[0], ifstream::in);
        if (!infile.good())
        {
            cerr << "Couldn't open file " << filenames[0] << ".\n";
            exit(1);
        }

        // Get the reference sequence length
        infile.seekg(0, ios::end);
        refseqlen = infile.tellg();
        infile.seekg(0, ios::beg);

        // Read the reference sequence
        refseq = new Array(refseqlen, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(infile, filenames[0], refseq, refseqlen);
        infile.close();
    }
}

RLZ::~RLZ()
{
    delete refseq;
    delete sa;
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

void RLZ::store_sequence(ifstream &infile, char *filename,
                         Array *dest, uint64_t length)
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
            dest->setField(i, v);
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

void RLZ::compress()
{
    uint64_t i;
    char outfilename[1024];
    ifstream infile;
    ofstream outfile;

    for (i=1; i<numfiles; i++)
    {
        // Open sequence file to be compressed
        infile.open(filenames[i], ifstream::in);
        if (!infile.good())
        {
            cerr << "Couldn't open file " << filenames[i] << ".\n";
            exit(1);
        }

        // Open output file to write compressed sequence to
        sprintf(outfilename, "%s.fac", filenames[i]);
        outfile.open(outfilename, ofstream::out);
        if (!outfile.good())
        {
            cerr << "Couldn't open file " << outfilename << ".\n";
            exit(1);
        }

        relative_LZ_factorise(infile, filenames[i], outfile);

        infile.close();
        outfile.close();
    }
}


void RLZ::relative_LZ_factorise(ifstream& infile, char *filename,
                                ofstream& outfile)
{
}

void RLZ::relative_LZ_factorise(ifstream& infile, char *filename,
                                ofstream& outfile, bool state)
{
    int c;
    uint64_t i, len;
    size_t pl, pr, cl, cr;
    bool runofns;

    i = 0;
    st->Root(&pl, &pr);
    len = 0;
    runofns = false;
    while (1)
    {
        c = infile.get();
        if (infile.eof())
        {
            break;
        }
        if (nucl_to_int[c] == 0)
        {
            cerr << "Invalid symbol " << c << " at position " << i;
            cerr << " of sequence in " << filename << ".\n";
            exit(1);
        }

        if (c == 'n')
        {
            if (!runofns && len > 0)
            {
                cout << st->Locate(pl,pl) << ' ' << len << endl;
                st->Root(&pl, &pr);
                len = 0; 
            }
            runofns = true;
            len++;
        }
        else
        {
            if (runofns)
            {
                cout << refseqlen << ' ' << len << endl;
                runofns = false;
                len = 0;
            }

            st->Child(pl, pr, (unsigned char)c, &cl, &cr);
            cout << (char)c << ' ' << cl << ' ' << cr << endl;
            
            if (len == 12)
            {
                st->FChild(pl, pr, &cl, &cr);
                cout << cl << ' ' << cr << endl;
                break;
            }

            if (cl == (uint64_t)(-1) || cr == (uint64_t)(-1))
            {
                for (i=pl; i<=pr; i++)
                {
                    cout << st->Locate(i,i) << endl;
                }
                //cout << st->Locate(pl,pl) << ' ' << len << endl;
                st->Root(&pl, &pr);
                len = 0;
                break;
            }
            else
            {
                len++;
                pl = cl;
                pr = cr;
            }
        }

        i++;
    }
}
