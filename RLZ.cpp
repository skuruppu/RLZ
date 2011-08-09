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


RLZCompress::RLZCompress(char **filenames, uint64_t numfiles, 
                         char encoding, bool isshort)
{
    this->filenames = filenames;
    this->numfiles = numfiles;
    this->encoding = encoding;
    this->isshort = isshort;

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

    // Calculate the log of the reference sequence length
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;
}

RLZCompress::~RLZCompress()
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
                                                   isshort,
                                                   logrefseqlen);

        relative_LZ_factorise(infile, filenames[i], *facwriter);

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

    try
    {
        // Read factors until EOF
        while (facreader.read_factor(&pos, &len))
        {
            // Run length encoded Ns
            if (pos == refseqlen)
            {
                for (i=0; i<len; i++)
                {
                    outfile << 'n';
                }
                continue;
            }
            // Standard factor
            for (i=pos; i<pos+len; i++)
            {
                outfile << (char)int_to_nucl[refseq->getField(i)];
            }
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
                           bool isshort, uint64_t maxposbits)
{
    // The encoding format is tbsl---- where t=text, b=binary,
    // s=shortfac encoding, l=liss encoding
    if (encoding == 'b')
    {
        // Output the type of encoding
        // 01100000
        if (isshort)
            outfile << (unsigned char)96;
        // 01000000
        else
            outfile << (unsigned char)64;

        facwriter = new FactorWriterBinary(outfile, isshort,
                                           maxposbits);
    }
    else if (encoding == 't')
    {
        // Output the type of encoding
        // 10100000
        if (isshort )
            outfile << (unsigned char)160 << endl;
        // 10000000
        else
            outfile << (unsigned char)128 << endl;

        facwriter = new FactorWriterText(outfile, isshort);
    }
    else
    {
        cerr << "Unknown encoding type.\n";
        exit(1);
    }
}

void FactorWriter::write_factor(uint64_t pos, uint64_t len)
{
    facwriter->write_factor(pos, len);
}

FactorWriter::~FactorWriter()
{
    delete facwriter;
}

FactorWriterText::FactorWriterText(ofstream& outfile, bool isshort) :
    outfile(outfile), isshort(isshort) {}

FactorWriterBinary::FactorWriterBinary(ofstream& outfile, bool isshort,
                                       uint64_t maxposbits)
{
    bwriter = new BitWriter(outfile);
    gcoder = new GolombCoder(*bwriter, 64);

    this->isshort = isshort;
    this->maxposbits = maxposbits;

    // Output the Golomb coding parameter
    bwriter->int_to_binary(64, 8);
}

FactorWriterBinary::~FactorWriterBinary()
{
    bwriter->flush();

    delete bwriter;
    delete gcoder;
}

void FactorWriterText::write_factor(uint64_t pos, uint64_t len)
{
    outfile << pos << ' ' << len << endl;
}

void FactorWriterBinary::write_factor(uint64_t pos, uint64_t len)
{
    bwriter->int_to_binary(pos, maxposbits);
    gcoder->golomb_encode(len);
}

FactorReader::FactorReader() :
    facreader(NULL) {}

FactorReader::FactorReader(ifstream& infile, uint64_t maxposbits)
{
    infile.exceptions(ifstream::failbit | ifstream::badbit |
                      ifstream::eofbit);

    // Read a byte to figure out what encoding was used
    unsigned char encodings;
    infile >> encodings;

    // Plain text 10000000
    if (encodings & ((unsigned)1<<7))
    {
        facreader = new FactorReaderText(infile);
        // Read the new line
        infile.get();
    }
    // Binary output
    else if (encodings & ((unsigned)1<<6))
    {
        facreader = new FactorReaderBinary(infile, maxposbits);
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

bool FactorReader::read_factor(uint64_t *pos, uint64_t *len)
{
    return facreader->read_factor(pos, len);
}

FactorReaderText::FactorReaderText(ifstream& infile) :
    infile(infile) {}
        
bool FactorReaderText::read_factor(uint64_t *pos, uint64_t *len)
{
    try
    {
        infile >> *pos >> *len;
    }
    catch (ifstream::failure e)
    {
        // EOF detected
        if (infile.eof())
        {
            *pos = NULL; *len = NULL;
            return false;
        }
        // Some other IO error
        throw e;
    }

    return true;
}

FactorReaderBinary::FactorReaderBinary(ifstream& infile, 
                                       uint64_t maxposbits)
{
    breader = new BitReader(infile);
    uint64_t divisor = breader->binary_to_int(8);
    gdecoder = new GolombCoder(*breader, (unsigned int)divisor);

    this->maxposbits = maxposbits;
}

FactorReaderBinary::~FactorReaderBinary()
{
    delete breader;
    delete gdecoder;
}

bool FactorReaderBinary::read_factor(uint64_t *pos, uint64_t *len)
{
    try
    {
        *pos = breader->binary_to_int(maxposbits);
        *len = gdecoder->golomb_decode();
    }
    // EOF reached
    catch (BitsEOFException& eofexp)
    {
        *pos = NULL; *len = NULL;
        return false;
    }
    return true;
}
