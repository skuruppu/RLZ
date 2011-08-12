#include <divsufsort64.h>
#include "RLZ.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

const char *RLZ::NUCLALPHA = "acgnt";

unsigned char nucl_to_int[256] = {0};
unsigned char int_to_nucl[RLZ::NUCLALPHASIZE+1];
unsigned char int_to_2bpb[RLZ::NUCLALPHASIZE+1];
unsigned char bpb_to_char[RLZ::NUCLALPHASIZE+1];

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

    int_to_2bpb[1] = 0; // a
    int_to_2bpb[2] = 1; // c
    int_to_2bpb[3] = 2; // g
    int_to_2bpb[4] = -1; // n (invalid)
    int_to_2bpb[5] = 3; // t

    bpb_to_char[0] = 'a';
    bpb_to_char[1] = 'c';
    bpb_to_char[2] = 'g';
    bpb_to_char[3] = 't';
}


RLZCompress::RLZCompress(char **filenames, uint64_t numfiles, 
                         char encoding, bool isshort, bool isliss)
{
    this->filenames = filenames;
    this->numfiles = numfiles;
    this->encoding = encoding;
    this->isshort = isshort;
    this->isliss = isliss;
    this->st = NULL;

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
        refseq = new Array(refseqlen+1, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(sequence, filenames[0], refseq, refseqlen);

        // Construct suffix array
        uint64_t *sufarray = new uint64_t[refseqlen+1];
        if (divsufsort64((sauchar_t*)sequence, (saidx64_t*)sufarray,
            refseqlen+1) < 0)
        {
            cerr << "Error in constructing suffix array.\n";
            exit(1);
        }

        sa = new Array(refseqlen+1, refseqlen);
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
        refseq = new Array(refseqlen+1, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(infile, filenames[0], refseq, refseqlen);
        infile.close();

    }

    // Calculate the log of the reference sequence length
    i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;
}

RLZCompress::RLZCompress(char **filenames, uint64_t numfiles, 
                         bool state)
{
    this->filenames = filenames;
    this->numfiles = numfiles;
    this->sa = NULL;

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
        refseq = new Array(refseqlen+1, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(sequence, filenames[0], refseq, refseqlen);

        // Construct suffix tree
        SuffixTreeY *sty = new SuffixTreeY(sequence, refseqlen+1, DAC,
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
        refseq = new Array(refseqlen+1, ((unsigned)1<<BITSPERBASE)-1);
        store_sequence(infile, filenames[0], refseq, refseqlen);
        infile.close();
    }

    // Calculate the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;
}

RLZCompress::~RLZCompress()
{
    delete refseq;
    if (sa != NULL) delete sa;
    if (st != NULL) delete st;
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
    dest->setField(i, 0);

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
    dest->setField(i, 0);

    return;
}

void RLZCompress::compress()
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

        // Intialise the factor writer
        FactorWriter *facwriter = new FactorWriter(outfile, encoding,
                                                   isshort, isliss,
                                                   this->refseq,
                                                   refseqlen,
                                                   logrefseqlen);

        relative_LZ_factorise(infile, filenames[i], *facwriter);
        //relative_LZ_factorise(infile, filenames[i], outfile, true);

        delete facwriter;

        infile.close();
        outfile.close();
    }
}


void RLZCompress::relative_LZ_factorise(ifstream& infile, 
                                        char *filename,
                                        FactorWriter &facwriter)
{
    int c;
    uint64_t i, len;
    uint64_t pl, pr, cl, cr;
    bool runofns;

    i = 0;
    runofns = false;
    pl = 0; pr = refseqlen; len = 0;
    while (1)
    {
        // EOF reached
        if ((c = infile.get()) == EOF) break;

        // Invalid symbol in file
        if (nucl_to_int[c] == 0)
        {
            cerr << "Invalid symbol " << c << " at position " << i;
            cerr << " of sequence in " << filename << ".\n";
            exit(1);
        }

        // Part of a run of Ns
        if (c == 'n')
        {
            // A new run of Ns so print earlier factor and reset suffix
            // array boundaries
            if (!runofns && len > 0)
            {
                facwriter.write_factor(sa->getField(pl), len);
                pl = 0; pr = refseqlen; len = 0; 
            }
            runofns = true;
            len++;
        }
        else
        {
            // A run of Ns just ended so print the factor
            if (runofns)
            {
                facwriter.write_factor(refseqlen, len);
                runofns = false;
                len = 0;
            }

            // Try to extend the current match
            sa_binary_search(pl, pr, nucl_to_int[c], len, &cl, &cr);
            
            // Couldn't extend current match so print factor
            if (cl == (uint64_t)(-1) || cr == (uint64_t)(-1))
            {
                facwriter.write_factor(sa->getField(pl), len);
                infile.unget();
                pl = 0; pr = refseqlen; len = 0;
            }
            // Set the suffix array boundaries to narrow the search for
            // the next symbol
            else
            {
                pl = cl; pr = cr; len++;
            }
        }
        i++;
    }
    
    // Print last factor
    if (len > 0)
    {
        if (runofns)
            facwriter.write_factor(refseqlen, len);
        else
            facwriter.write_factor(sa->getField(pl), len);
    }

    // Call this to indicate the end of input to the factor writer
    facwriter.finalise();
}

void RLZCompress::relative_LZ_factorise(ifstream& infile, 
                                        char *filename,
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
        // EOF reached
        if ((c = infile.get()) == EOF) break;

        // Invalid symbol in file
        if (nucl_to_int[c] == 0)
        {
            cerr << "Invalid symbol " << c << " at position " << i;
            cerr << " of sequence in " << filename << ".\n";
            exit(1);
        }

        // Part of a run of Ns
        if (c == 'n')
        {
            // A new run of Ns so print earlier factor and go back to
            // the root of the suffix tree
            if (!runofns && len > 0)
            {
                cout << st->Locate(pl,pl) << ' ' << len << endl;
                st->Root(&pl, &pr); len = 0; 
            }
            runofns = true;
            len++;
        }
        else
        {
            // A run of Ns just ended so print the factor
            if (runofns)
            {
                cout << refseqlen << ' ' << len << endl;
                runofns = false;
                len = 0;
            }

            // The previous suffix tree branch taken covers more than
            // one symbol
            if (len < st->SDepth(pl, pr))
            {
                // Couldn't extend current match so print factor
                if ((unsigned char)c != st->Letter(pl, pr, len+1))
                {
                    cout << st->Locate(pl,pl) << ' ' << len << endl;
                    infile.unget(); i--;
                    // Reset the search range to the entire suffix array
                    st->Root(&pl, &pr); len = 0;
                }
                else
                {
                    len++;
                }
            }
            // Need to traverse to a new child
            else
            {
                st->Child(pl, pr, (unsigned char)c, &cl, &cr);
            
                // Couldn't extend current match so print factor
                if (cl == (uint64_t)(-1) || cr == (uint64_t)(-1))
                {
                    cout << st->Locate(pl,pl) << ' ' << len << endl;
                    infile.unget(); i--;
                    // Reset the search range to the entire suffix array
                    st->Root(&pl, &pr); len = 0;
                }
                // Set the suffix array boundaries to narrow the search
                // for the next symbol
                else
                {
                    pl = cl; pr = cr; len++;
                }
            }
        }
        i++;
    }

    // Print last factor
    if (len > 0)
    {
        if (runofns)
            cout << refseqlen << ' ' << len << endl;
        else
            cout << st->Locate(pl,pl) << ' ' << len << endl;
    }
}

void RLZCompress::sa_binary_search(uint64_t pl, uint64_t pr, int c,
                                   uint64_t offset, uint64_t *cl,
                                   uint64_t *cr)
{
    uint64_t low, high, mid; 
    int midval, midvalleft, midvalright;

    // Binary search left
    low = pl; high = pr;
    while (low <= high)
    {
        mid = (low + high) >> 1;

        midval = refseq->getField(sa->getField(mid)+offset);
        // Move left boundary to the middle
        if (midval < c)
            low = mid + 1;
        // Move right boundary to the middle
        else if (midval > c)
            high = mid - 1;
        else
        {
            // Mid is at the left boundary
            if(mid == pl)
            {
                *cl = mid;
                break;
            }
            midvalleft = refseq->getField(sa->getField(mid-1)+offset);
            // Discard mid and values to the right of mid
            if(midvalleft == midval)
                high = mid - 1;
            // Left-most occurrence found
            else
            {
                *cl = mid;
                break;
            }
        }
    }

    // Key not found so return not found symbols
    if (low > high)
    {
        *cl = (uint64_t)(-1);
        *cr = (uint64_t)(-1);
        return;
    }

    // Binary search right
    low = *cl; high = pr;
    while (low <= high)
    {
        mid = (low + high) >> 1;


        midval = refseq->getField(sa->getField(mid)+offset);
        // Move left bounary to the middle
        if (midval < c)
            low = mid + 1;
        // Move right boundary to the middle
        else if (midval > c)
            high = mid - 1;
        else
        { 
            // Rightmost occurrence of key found
            if(mid == pr)
            {
                *cr = mid;
                break;
            }
            midvalright = refseq->getField(sa->getField(mid+1)+offset);
            // Discard mid and the ones to the left of mid
            if(midvalright == midval)
                low = mid + 1; 
            // Rightmost occurrence of key found
            else 
            {
                *cr = mid;
                break;
            }
        }
    }

    // Key not found so return not found symbols
    if (low > high)
    {
        *cl = (uint64_t)(-1);
        *cr = (uint64_t)(-1);
        return;
    }

    return;
}

RLZDecompress::RLZDecompress(char **filenames, uint64_t numfiles)
{

    this->filenames = filenames;
    this->numfiles = numfiles;

    initialise_nucl_converters();

    // Open reference sequence file
    ifstream infile(filenames[0], ifstream::in);
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
    refseq = new Array(refseqlen+1, ((unsigned)1<<BITSPERBASE)-1);
    store_sequence(infile, filenames[0], refseq, refseqlen);
    infile.close();

    // Calculate the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;
}

RLZDecompress::~RLZDecompress()
{
    delete refseq;
}

void RLZDecompress::decompress()
{
    uint64_t i;
    char infilename[1024], outfilename[1024];
    ifstream infile;
    ofstream outfile;

    for (i=1; i<numfiles; i++)
    {
        // Open sequence file to be decompressed
        sprintf(infilename, "%s.fac", filenames[i]);
        infile.open(infilename, ifstream::in);
        if (!infile.good())
        {
            cerr << "Couldn't open file " << infilename << ".\n";
            exit(1);
        }

        // Open output file to write decompressed sequence to
        sprintf(outfilename, "%s.dec", filenames[i]);
        outfile.open(outfilename, ofstream::out);
        if (!outfile.good())
        {
            cerr << "Couldn't open file " << outfilename << ".\n";
            exit(1);
        }

        // Intialise the factor reader
        FactorReader *facreader = new FactorReader(infile, logrefseqlen);

        relative_LZ_defactorise(*facreader, infilename, outfile);

        delete facreader;

        infile.close();
        outfile.close();
    }
}

void RLZDecompress::relative_LZ_defactorise(FactorReader& facreader,
                                            char *filename,
                                            ofstream& outfile)
{
    uint64_t pos, len, i;
    vector<char> substr;

    try
    {
        // Read factors until EOF
        while (facreader.read_factor(&pos, &len, substr))
        {
            // A short factor
            if (substr.size() == len)
            {
                for (i=0; i<len; i++)
                {
                    outfile << substr.at(i);
                }
            }
            // Run length encoded Ns
            else if (pos == refseqlen)
            {
                for (i=0; i<len; i++)
                {
                    outfile << 'n';
                }
            }
            // Standard factor
            else
            {
                for (i=pos; i<pos+len; i++)
                {
                    outfile << (char)int_to_nucl[refseq->getField(i)];
                }
            }
            substr.clear();
        }
    }
    catch (exception e)
    {
        cerr << e.what() << endl;
        cerr << "Could not read from file " << filename << ".\n";
        exit(1);
    }
}

FactorWriter::FactorWriter() :
    facwriter(NULL) {}

FactorWriter::FactorWriter(ofstream& outfile, char encoding,
                           bool isshort, bool isliss, Array *refseq, 
                           uint64_t refseqlen, uint64_t logrefseqlen)
{
    unsigned char encbyte = 0;

    this->isshort = isshort;
    this->isliss = isliss;

    // The encoding format is tbsl---- where t=text, b=binary,
    // s=shortfac encoding, l=liss encoding

    if (isshort)
        encbyte |= (unsigned)1<<5;

    if (isliss)
        encbyte |= (unsigned)1<<4;

    if (encoding == 'b')
    {
        encbyte |= (unsigned)1<<6;

        // Output the type of encoding
        outfile << encbyte;

        facwriter = new FactorWriterBinary(outfile, isshort, isliss,
                                           refseq, refseqlen,
                                           logrefseqlen);
    }
    else if (encoding == 't')
    {
        encbyte |= (unsigned)1<<7;

        // Output the type of encoding
        outfile << encbyte << endl;

        facwriter = new FactorWriterText(outfile, isshort, isliss,
                                         refseq, refseqlen,
                                         logrefseqlen);
    }
    else
    {
        cerr << "Unknown encoding type.\n";
        exit(1);
    }
    
}

void FactorWriter::write_factor(uint64_t pos, uint64_t len)
{
    // If LISS encoding then can't output factors yet so just store them
    if (isliss)
    {
        positions.push_back(pos);
        lengths.push_back(len);
        return;
    }

    facwriter->write_factor(pos, len);
}

void FactorWriter::finalise()
{
    vector<uint64_t> liss;
    uint64_t i, j;
    uint64_t prevpos, cumlen;

    if (isliss)
    {
        find_LISS(positions, liss);
        i = j = 0;
        // Write standard factors until the first LISS factor
        while (j < liss[i])
        {
            facwriter->write_factor(positions[j], lengths[j], false);
            j++;
        }
        facwriter->write_factor(positions[j], lengths[j], true);
        // Initialise the position and length accumulators to keep track
        // of the last LISS factor encountered
        prevpos = positions[j];
        cumlen = lengths[j];
        i++; j++;
        while (i < liss.size())
        {
            // Write standard factors
            while (j < liss[i])
            {
                facwriter->write_factor(positions[j], lengths[j],
                                        false);
                cumlen += lengths[j];
                j++;
            }
            // Write LISS factor
            facwriter->write_factor(positions[j]-(prevpos+cumlen),
                                    lengths[j], true);
            prevpos = positions[j];
            cumlen = lengths[j];
            i++; j++;
        }
        // Write the remaining standard factors
        while (j<positions.size())
        {
            facwriter->write_factor(positions[j], lengths[j], false);
            j++;
        }
    }

    facwriter->finalise();
}

FactorWriter::~FactorWriter()
{
    delete facwriter;
}

// Finds longest strictly increasing subsequence (LISS).
// O(n log k) algorithm, k is the length of the LISS.
void FactorWriter::find_LISS(vector<uint64_t>& a, vector<uint64_t>& b)
{
    uint64_t asize = a.size();
    vector<uint64_t> p(asize);
    uint64_t u, v, c;
    uint64_t i;
   
    i = 0;
    // Ignore factors that represent Ns
    while (a[i] == refseqlen)
        i++;
    // Add the first factor
    b.push_back(i);
   
    for (; i < asize; i++) 
    {
        // Ignore factors that represent Ns
        if (a[i] == refseqlen)
            continue;

        // If the last inserted index in b has a position less than the
        // position at the current index then no need to binary search so
        // just store the values 
        if (a[b.back()] < a[i]) 
        {
            p[i] = b.back();
            b.push_back(i);
            continue;
        }
      
        // Binary search till you find the place where the current
        // position will fit in b
        for (u = 0, v = b.size()-1; u < v;) 
        {
            c = (u + v) / 2;
            if (a[b[c]] < a[i]) 
                u=c+1; 
            else 
                v=c;
        }
      
        // Found a spot to insert the current position 
        if (a[i] < a[b[u]]) 
        {
            if (u > 0) 
                p[i] = b[u-1];
            b[u] = i;
        }	
    }

    // Back-track to figure out which indices are part of the LISS
    for (u = b.size(), v = b.back(); u>0; v = p[v], u--) 
    {
        b[u-1] = v;
    }
}

FactorWriterText::FactorWriterText(ofstream& outfile, bool isshort, 
                                   bool isliss, Array *refseq, uint64_t
                                   refseqlen, uint64_t logrefseqlen) :
    outfile(outfile)
{
    this->refseq = refseq;
    this->refseqlen = refseqlen; 
    this->logrefseqlen = logrefseqlen;
    this->isshort = isshort;
    this->isliss = isliss;

    if (isshort)
    {
        // 2*len+len/GOLOMBDIVSHORT+(LOG2GOLOMBDIVSHORT+1) <
        // logrefseqlen+len/GOLOMBDIV+(LOG2GOLOMBDIV+1)
        SHORTFACTHRESH = (64.0/135)*(logrefseqlen+3);
    
    }

    // Initialise a flag needed to start LISS encoding
    if (isliss)
    {
        firstliss = true;
    }
}

void FactorWriterText::write_factor(uint64_t pos, uint64_t len)
{
    uint64_t i;

    // 2*len+len/8+4 < logrefseqlen+len/64+7
    if (isshort && pos!=refseqlen && len <= SHORTFACTHRESH)
    {
        for (i=pos; i<pos+len; i++)
        {
            // Don't short factor encode if there are 'n's in the factor
            if (int_to_nucl[refseq->getField(i)] == 'n')
            {
                outfile << pos << ' ' << len << endl;
                return;
            }
        }

        // Output a short factor, substring followed by the length
        for (i=pos; i<pos+len; i++)
        {
            outfile << int_to_nucl[refseq->getField(i)];
        }
        outfile << ' ' << len << endl;
        return;
    }

    // Just a standard factor
    outfile << pos << ' ' << len << endl;
}

void FactorWriterText::write_factor(uint64_t pos, uint64_t len,
                                    bool lissfac)
{
    if (lissfac)
    {
        // Write the first LISS factor as a standard factor
        if (firstliss)
        {
            outfile << "* " << pos << ' ' << len << endl;
            firstliss = false;
            return;
        }

        // Output the diff and the length of the factor
        outfile << "* " << (int64_t)pos << ' ' << len << endl;
        return;
    }

    // Standard factor
    write_factor(pos, len);
}

void FactorWriterText::finalise()
{

}

FactorWriterBinary::FactorWriterBinary(ofstream& outfile, bool isshort,
                                       bool isliss, Array *refseq,
                                       uint64_t refseqlen,
                                       uint64_t logrefseqlen)
{
    bwriter = new BitWriter(outfile);
    gcoder = new GolombCoder(*bwriter, GOLOMBDIV);
    gcodershort = new GolombCoder(*bwriter, GOLOMBDIVSHORT);

    this->isshort = isshort;
    this->isliss = isliss;
    this->refseq = refseq;
    this->refseqlen = refseqlen; 
    this->logrefseqlen = logrefseqlen;

    // Output the Golomb coding parameter
    bwriter->int_to_binary(GOLOMBDIV, 8);

    // Ouptut the Golomb coding parameter for short ints
    if (isshort)
    {
        bwriter->int_to_binary(GOLOMBDIVSHORT, 8);

        // 2*len+len/GOLOMBDIVSHORT+(LOG2GOLOMBDIVSHORT+1) <
        // logrefseqlen+len/GOLOMBDIV+(LOG2GOLOMBDIV+1)
        SHORTFACTHRESH = (64.0/135)*(logrefseqlen+3);
    }

    // Initialise a flag needed to start LISS encoding
    if (isliss)
    {
        firstliss = true;
    }
}

FactorWriterBinary::~FactorWriterBinary()
{
    delete bwriter;
    delete gcoder;
    delete gcodershort;
}

void FactorWriterBinary::write_factor(uint64_t pos, uint64_t len)
{
    uint64_t i;

    // 2*len+len/8+4 < logrefseqlen+len/64+7
    if (isshort && pos!=refseqlen && len <= SHORTFACTHRESH)
    {
        for (i=pos; i<pos+len; i++)
        {
            // Don't short factor encode if there are 'n's in the factor
            if (int_to_nucl[refseq->getField(i)] == 'n')
            {
                // Indicate this is not a short factor
                bwriter->write_bit(1);
                // Encode pos and len pair
                bwriter->int_to_binary(pos, logrefseqlen);
                gcoder->golomb_encode(len);
                return;
            }
        }

        // Valid short factor so first output a 0 bit to indicate this
        bwriter->write_bit(0);
        // Encode length followed by 2bpb substring
        gcodershort->golomb_encode(len);
        for (i=pos; i<pos+len; i++)
        {
            bwriter->int_to_binary(int_to_2bpb[refseq->getField(i)], 2);
        }
        return;
    }

    // Not a short factor so output 1 bit to indicate this is not a
    // short factor
    if (isshort) bwriter->write_bit(1);

    // Output factor as a standard factor
    bwriter->int_to_binary(pos, logrefseqlen);
    gcoder->golomb_encode(len);
}

void FactorWriterBinary::write_factor(uint64_t pos, uint64_t len,
                                      bool lissfac)
{
    if (lissfac)
    {
        // A bit to indicate that it's an LISS factor
        bwriter->write_bit(1);

        // Write the first LISS factor as a standard factor
        if (firstliss)
        {
            if (isshort) bwriter->write_bit(1);
            // Output factor as a standard factor
            bwriter->int_to_binary(pos, logrefseqlen);
            gcoder->golomb_encode(len);
            firstliss = false;
            return;
        }

        // The next LISS position can be inferred from the previous LISS
        // position and the cumulative length of factors inbetween
        if ((int64_t)pos == 0)
        {
            bwriter->write_bit(0);
            // Only need the length
            gcoder->golomb_encode(len);
        }
        // The next LISS position is to the left of the inferred LISS
        // position
        else if ((int64_t)pos < 0)
        {
            // 1 bit to say position is a diff
            bwriter->write_bit(1);
            // 0 bit to say the diff is negative
            bwriter->write_bit(0);
            gcoder->golomb_encode(abs((int64_t)pos));
            gcoder->golomb_encode(len);
        }
        // The next LISS position is to the right of the inferred LISS
        // position
        else
        {
            // 1 bit to say position is a diff
            bwriter->write_bit(1);
            // 1 bit to say the diff is positive
            bwriter->write_bit(1);
            gcoder->golomb_encode(pos);
            gcoder->golomb_encode(len);
        }

        return;
    }

    // A bit to indicate that it's not encoded as an LISS factor
    bwriter->write_bit(0);
    // Standard factor
    write_factor(pos, len);
}

void FactorWriterBinary::finalise()
{
    bwriter->flush();
}

FactorReader::FactorReader() :
    facreader(NULL) {}

FactorReader::FactorReader(ifstream& infile, uint64_t logrefseqlen)
{
    infile.exceptions(ifstream::failbit | ifstream::badbit |
                      ifstream::eofbit);

    // Read a byte to figure out what encoding was used
    unsigned char encodings;
    infile >> encodings;

    bool isshort = false;
    bool isliss = false;

    // Short factor encoding enabled
    if (encodings & ((unsigned)1<<5))
        isshort = true;

    // LISS encoding enabled
    if (encodings & ((unsigned)1<<4))
        isliss = true;

    // Plain text 10000000
    if (encodings & ((unsigned)1<<7))
    {
        facreader = new FactorReaderText(infile, isshort, isliss);
        // Read the new line
        infile.get();
    }
    // Binary output
    else if (encodings & ((unsigned)1<<6))
    {
        facreader = new FactorReaderBinary(infile, logrefseqlen,
                                           isshort, isliss);
    }
    else
    {
        cerr << "Invalid encoding in file.\n";
        exit(1);
    }
}

FactorReader::~FactorReader()
{
    delete facreader;
}

bool FactorReader::read_factor(uint64_t *pos, uint64_t *len,
                               vector<char>& substr)
{
    return facreader->read_factor(pos, len, substr);
}

FactorReaderText::FactorReaderText(ifstream& infile, bool isshort,
                                   bool isliss) :
    infile(infile), isshort(isshort), isliss(isliss)
{
    // Initialise variables needed for LISS factor encoding
    if (isliss)
    {
        firstliss = true;
        prevpos = cumlen = 0;
    }
}
        
bool FactorReaderText::read_factor(uint64_t *pos, uint64_t *len,
                                   vector<char>& substr)
{
    int c;
    int64_t diff;

    try
    {
        // LISS factor
        if (isliss && infile.peek() == '*')
        {
            // Remove the * and the space after it
            assert(infile.get() == '*');
            assert(infile.get() == ' ');
            // First LISS factor so it's a standard factor
            if (firstliss)
            {
                infile >> *pos >> *len;
                firstliss = false;
            }
            else
            {
                // Retrieve position diff and length
                infile >> diff >> *len;
                *pos = prevpos + cumlen + diff; 
            }
            prevpos = *pos;
            cumlen = 0;
        }
        // Short factor
        else if (isshort && nucl_to_int[infile.peek()] > 0)
        {
            while ((c = infile.get()) != ' ')
                substr.push_back(c);
            infile >> *len;
        }
        // Standard factor
        else
        {
            infile >> *pos >> *len;
        }
        // Accumulate in between factor lengths
        if (isliss) cumlen += *len;
        // Get rid of the newline
        assert(infile.get() == '\n');
    }
    catch (ifstream::failure e)
    {
        // EOF detected
        if (infile.eof())
        {
            *pos = NULL; *len = NULL; substr.clear();
            return false;
        }
        // Some other IO error
        throw e;
    }

    return true;
}

FactorReaderBinary::FactorReaderBinary(ifstream& infile, 
                                       uint64_t logrefseqlen,
                                       bool isshort, bool isliss)
{
    uint64_t divisor;

    breader = new BitReader(infile);
    divisor = breader->binary_to_int(8);
    gdecoder = new GolombCoder(*breader, (unsigned int)divisor);

    this->logrefseqlen = logrefseqlen;
    this->isshort = isshort;
    this->isliss = isliss;

    // Get the divisor for Golomb decoding lengths for short factors
    if (isshort)
    {
        divisor = breader->binary_to_int(8);
        gdecodershort = new GolombCoder(*breader, (unsigned int)divisor);
    }

    // Initialise variables needed for LISS factor encoding
    if (isliss)
    {
        firstliss = true;
        prevpos = cumlen = 0;
    }
}

FactorReaderBinary::~FactorReaderBinary()
{
    delete breader;
    delete gdecoder;

    if (isshort)
        delete gdecodershort;
}

bool FactorReaderBinary::read_factor(uint64_t *pos, uint64_t *len,
                                     vector<char>& substr)
{
    uint64_t i;
    int64_t diff;

    try
    {
        // LISS factor
        if (isliss && breader->read_bit() == 1)
        {
            // First LISS factor so it's a standard factor
            if (firstliss)
            {
                // Get rid of the short factor encoding bit
                if (isshort) assert(breader->read_bit() == 1);
                *pos = breader->binary_to_int(logrefseqlen);
                *len = gdecoder->golomb_decode();
                firstliss = false;
            }
            else
            {
                // No diff between predicted and actual position
                if (breader->read_bit() == 0)
                {
                    diff = 0;
                }
                else
                {
                    // Retrieve position diff
                    diff = (breader->read_bit() == 0) ? -1 : 1;
                    diff *= gdecoder->golomb_decode();
                }
                *pos = prevpos + cumlen + diff;
                *len = gdecoder->golomb_decode();
            }
            prevpos = *pos;
            cumlen = 0;
        }
        // Short factor
        else if (isshort && breader->read_bit() == 0)
        {
            *len = gdecodershort->golomb_decode();
            for (i=0; i<*len; i++)
            {
                substr.push_back(bpb_to_char[breader->binary_to_int(2)]);
            }
        }
        else
        {
            // Standard factor
            *pos = breader->binary_to_int(logrefseqlen);
            *len = gdecoder->golomb_decode();
        }
        // Accumulate in between factor lengths
        if (isliss) cumlen += *len;
    }
    // EOF reached
    catch (BitsEOFException& eofexp)
    {
        *pos = NULL; *len = NULL;
        return false;
    }
    return true;
}

FactorWriterIndex::FactorWriterIndex(ofstream& outfile, 
                                     cds_utils::Array *refseq, 
                                     uint64_t refseqlen, 
                                     uint64_t logrefseqlen) :
    outfile(outfile)
{
    this->refseq = refseq;
    this->refseqlen = refseqlen;
    this->logrefseqlen = logrefseqlen;

    bwriter = new BitWriter(outfile);

    numfacs = 0;
    cumlen = 0;
}

FactorWriterIndex::~FactorWriterIndex()
{
    delete bwriter;
}

void FactorWriterIndex::write_factor(uint64_t pos, uint64_t len)
{
    uint64_t i;

    facstarts.push_back(true);
    for (i=1; i<len; i++)
    {
        facstarts.push_back(false);
    }

    bwriter->int_to_binary(pos, logrefseqlen);

    numfacs++;
}

void FactorWriterIndex::finalise()
{
    cumseqlens.push_back(cumlen);

    cumlen = 0;
}
