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

#include <divsufsort64.h>
#include <BitSequenceSDArray.h>
#include <BitSequenceRRR.h>
#include "RLZ.h"
#include "alphabet.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

RLZCompress::RLZCompress(char **filenames, uint64_t numfiles, 
                         char encoding, bool isshort, bool isliss)
{
    this->filenames = filenames;
    this->numfiles = numfiles;
    this->encoding = encoding;
    this->isshort = isshort;
    this->isliss = isliss;
    this->idxname = NULL;
    this->displayonly = false;
    this->usecsa = false;
    this->st = NULL;

    read_refseq_and_sa();
}

RLZCompress::RLZCompress(char **filenames, uint64_t numfiles,
                         char *idxname, bool displayonly, bool usecsa)
{
    this->filenames = filenames;
    this->numfiles = numfiles;
    this->encoding = 'i';
    this->isshort = false;
    this->isliss = false;
    this->idxname = idxname;
    this->displayonly = displayonly;
    this->usecsa = usecsa;
    this->st = NULL;

    read_refseq_and_construct_sa();
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

void RLZCompress::read_refseq_and_construct_sa()
{
    uint64_t i;

    initialise_nucl_converters();

    // Read reference sequence into memory since its needed by
    // suffix tree constructor
    char *sequence = NULL;
    size_t seqlen;
    if (loadText(filenames[0], &sequence, &seqlen))
    {
        cerr << "Couldn't read reference sequence.\n";
        exit(1);
    }

    // loadText places an extra byte at the end
    refseqlen = seqlen-1;

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

    // Calculate the log of the reference sequence length
    i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;

    //delete [] sequence;
    delete [] sufarray;
}

void RLZCompress::read_refseq_and_sa()
{
    initialise_nucl_converters();

    uint64_t i;
    char safilename[1024];
    sprintf(safilename, "%s.sa", filenames[0]);
    ifstream infile;
    infile.open(safilename, ifstream::in);
    // Need to construct suffix array
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
    FactorWriter *facwriter = NULL;

    if (encoding == 'i')
    {
        // Open output file to write compressed sequence to
        outfile.open(idxname, ofstream::out);
        if (!outfile.good())
        {
            cerr << "Couldn't open file " << outfilename << ".\n";
            exit(1);
        }

        facwriter = new FactorWriterIndex(outfile, refseq, sa,
                                          refseqlen, logrefseqlen,
                                          displayonly, usecsa);
    }

    for (i=1; i<numfiles; i++)
    {
        // Open sequence file to be compressed
        infile.open(filenames[i], ifstream::in);
        if (!infile.good())
        {
            cerr << "Couldn't open file " << filenames[i] << ".\n";
            exit(1);
        }

        // Intialise the factor writer
        if (encoding != 'i')
        {

            // Open output file to write compressed sequence to
            sprintf(outfilename, "%s.fac", filenames[i]);
            outfile.open(outfilename, ofstream::out);
            if (!outfile.good())
            {
                cerr << "Couldn't open file " << outfilename << ".\n";
                exit(1);
            }

            facwriter = new FactorWriter(outfile, encoding, isshort,
                                         isliss, this->refseq,
                                         refseqlen, logrefseqlen);
        }

        relative_LZ_factorise(infile, filenames[i], *facwriter);
        //relative_LZ_factorise(infile, filenames[i], outfile, true);

        if (encoding != 'i')
        {
            delete facwriter;
            outfile.close();
        }

        infile.close();
    }

    if (encoding == 'i')
    {
        ((FactorWriterIndex*)facwriter)->write_index();
        delete facwriter;
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
    facwriter.end_of_sequence();
}

void RLZCompress::relative_LZ_factorise(ifstream& infile, 
                                        char *filename,
                                        ofstream& outfile, bool state)
{
    int c;
    uint64_t i, j, len, depth=0, saval=0;
    size_t pl, pr, cl, cr;
    bool runofns;
    vector<unsigned char> substr;

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
                st->Root(&pl, &pr); len = 0; depth = 0; saval = 0;
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
            if (len < depth)
            {
                // Couldn't extend current match so print factor
                if ((unsigned char)c != int_to_nucl[refseq->getField(saval+len)])
                {
                    cout << saval << ' ' << len << endl;
                    infile.unget(); i--;
                    // Reset the search range to the entire suffix array
                    st->Root(&pl, &pr); len = 0; depth = 0; saval = 0;
                }
                else
                {
                    len++; j++;
                }
            }
            // Need to traverse to a new child
            else
            {
                st->Child(pl, pr, (unsigned char)c, &cl, &cr);
            
                // Couldn't extend current match so print factor
                if (cl == (uint64_t)(-1))
                {
                    cout << st->Locate(pl,pl) << ' ' << len << endl;
                    infile.unget(); i--;
                    // Reset the search range to the entire suffix array
                    st->Root(&pl, &pr); len = 0; depth = 0; saval = 0;
                }
                // Set the suffix array boundaries to narrow the search
                // for the next symbol
                else
                {
                    depth = st->SDepth(cl, cr);
                    if (len+1 < depth)
                        saval = st->Locate(cl,cl);
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
    while (low <= high && high != (uint64_t)-1)
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
    if (low > high || high == (uint64_t)-1)
    {
        *cl = (uint64_t)(-1);
        *cr = (uint64_t)(-1);
        return;
    }

    // Binary search right
    low = *cl; high = pr;
    while (low <= high && high != (uint64_t)-1)
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
    if (low > high || high == (uint64_t)-1)
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

void FactorWriter::end_of_sequence()
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

    facwriter->end_of_sequence();
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

void FactorWriterText::end_of_sequence()
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

void FactorWriterBinary::end_of_sequence()
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
            substr.clear();
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
        return false;
    }
    return true;
}

// Class that supports sorting of factors according to their start
// and end positions.
class SortClass
{
private:
    uint64_t cumseqlen;
    BitSequenceSDArray& facstarts;
    uint64_t numfacs;
public:
    SortClass(uint64_t cumseqlen, BitSequenceSDArray& facstarts,
              uint64_t numfacs) : cumseqlen(cumseqlen),
              facstarts(facstarts), numfacs(numfacs) {}
    bool operator() (const uint32_t i, const uint32_t j)
    {
        return sort_function(i, j);
    }

    // Factor sorting function
    bool sort_function(const uint32_t i, const uint32_t j)
    {
        uint32_t len1, len2;

        // Get the lengths of the two factors 
        if (i+1 == numfacs)
            len1 = cumseqlen - facstarts.select1(i+1);
        else
            len1 = facstarts.select1(i+2) - facstarts.select1(i+1);

        if (j+1 == numfacs)
            len2 = cumseqlen - facstarts.select1(j+1);
        else
            len2 = facstarts.select1(j+2) - facstarts.select1(j+1);

        // If positions are equal then order it on decreasing length
        return (len1 > len2);
    }

};

FactorWriterIndex::FactorWriterIndex(ofstream& outfile, 
                                     cds_utils::Array *refseq, 
                                     cds_utils::Array *sa, 
                                     uint64_t refseqlen, 
                                     uint64_t logrefseqlen,
                                     bool displayonly,
                                     bool usecsa) :
    outfile(outfile)
{
    this->refseq = refseq;
    this->sa = sa;
    this->refseqlen = refseqlen;
    this->logrefseqlen = logrefseqlen;
    this->displayonly = displayonly;
    this->usecsa = usecsa;

    bwriter = new BitWriter(outfile);

    numfacs = 0;
    cumlen = 0;

    // Fill the first two slots with zeros
    cumseqlens.push_back(0);
    cumseqlens.push_back(0);

    posarraylen = 1000;
    positions = new unsigned int[posarraylen];
    isstart = new BitString(refseqlen+1);

    if (!displayonly)
    {
        isend = new BitString(refseqlen);

        currseqfacnum = 0;

        nll = NULL;
        levelidx = NULL;
    }
}

void FactorWriterIndex::write_index()
{
    uint64_t i;

    // First write the data structure that are just necessary to answer
    // display() queries
    
    // Write the reference sequence
    refseq->save(outfile);
    cout << "refseq: " << refseq->getSize() << endl;

    // Create the compressed bit vector and write it
    BitSequenceSDArray compfacstarts(facstarts);
    compfacstarts.save(outfile);
    cout << "facstarts: " << compfacstarts.getSize() << endl;

    // Create the compressed bit vector for isstart
    BitSequenceRRR compisstart(*isstart);
    compisstart.save(outfile);
    cout << "isstart: " << compisstart.getSize() << endl;

    // Write out the positions but map them to consecutive integers
    // using isstart to find the positions that it should map to
    Array posarray(numfacs, refseqlen);
    for (i=0; i<numfacs; i++)
    {
        posarray.setField(i, get_field_64(positions, logrefseqlen, i));
    }
    posarray.save(outfile);
    cout << "positions: " << posarray.getSize() << endl;

    // Write out the cumulative sequence lengths
    Array cumseqlensarray(cumseqlens.size(), cumseqlens.back());
    for (i=0; i<cumseqlens.size(); i++)
        cumseqlensarray.setField(i, cumseqlens.at(i));
    cumseqlensarray.save(outfile);
    cout << "cumseqlens: " << cumseqlensarray.getSize() << endl;

    // Write out the data structures necessary to implement count() and
    // locate() queries
    if (!displayonly)
    {
        if (usecsa)
        {
            size_t seqlen = refseqlen+1;
            char *sequence1 = new char[seqlen];
            for (i=0; i<seqlen-1; i++)
                sequence1[i] = int_to_nucl[refseq->getField(i)];
            sequence1[i] = '\0';
            TextIndexCSA *csa = new TextIndexCSA((uchar*)sequence1, seqlen,
                                                 NULL);
            csa->save(outfile);
            cout << "csa: " << csa->getSize() << endl;
        }
        else
        {
            // Write out the suffix array
            sa->save(outfile);
            cout << "sa: " << sa->getSize() << endl;
        }

        /*
        // Construct suffix tree
        SuffixTreeY sty(sequence, refseqlen+1);
        sty.save(outfile);
        cout << "st: " << sty.getSize() << endl;
        delete [] sequence;
        */

        // Construct and write the nested level list
        construct_nested_level_list(compfacstarts);
        nll->save(outfile);
        levelidx->save(outfile);
        cout << "nll: " << nll->getSize()+(numlevels+1)*sizeof(uint32_t);
        cout << endl;

        // Create the compressed bit vector for isend
        BitSequenceRRR compisend(*isend);
        compisend.save(outfile);
        cout << "isend: " << compisend.getSize() << endl;

        // Create the compressed bit vector for factor start positions
        BitSequenceSDArray compseqstarts(seqstarts);
        compseqstarts.save(outfile);
        cout << "compseqstarts: " << compseqstarts.getSize() << endl;
    }
}

FactorWriterIndex::~FactorWriterIndex()
{
    delete bwriter;
    delete [] positions;
    delete isstart;
    if (!displayonly)
    {
        if (nll) delete nll;
        if (levelidx) delete levelidx;
        delete isend;
    }
}

void FactorWriterIndex::write_factor(uint64_t pos, uint64_t len)
{
    uint64_t i;

    // Set a bit at the factor start position and set the rest of the
    // bits to zero
    facstarts.push_back(true);
    for (i=1; i<len; i++)
        facstarts.push_back(false);

    // Check if there's enough memory to store the position and if not
    // allocate more memory
    if (numfacs*logrefseqlen/(sizeof(unsigned int)*8)+1 >= posarraylen)
    {
        unsigned int *newarray = new unsigned int[posarraylen*2];
        memcpy(newarray, positions, posarraylen*sizeof(unsigned int));
        delete [] positions;
        positions = newarray;
        posarraylen *= 2;
    }
    set_field_64(positions, logrefseqlen, numfacs, pos);

    isstart->setBit(pos);

    cumlen += len;
    numfacs++;

    if (!displayonly)
    {
        if (pos != refseqlen)
            isend->setBit(pos+len);
        currseqfacnum++;
    }
}

void FactorWriterIndex::end_of_sequence()
{
    // Store the cumulative length of this sequence
    cumseqlens.push_back(cumseqlens.back()+cumlen);

    // Reset the cumulative length in preparation for the next sequence
    cumlen = 0;

    if (!displayonly)
    {
        // Keep track of the factors at which new sequences start
        seqstarts.push_back(true);
        for (uint32_t i=1; i<currseqfacnum; i++)
            seqstarts.push_back(false);
        currseqfacnum = 0;
    }
}


void FactorWriterIndex::construct_nested_level_list
     (cds_static::BitSequenceSDArray& compfacstarts)
{
    vector<uint32_t> *posindices = new vector<uint32_t>[refseqlen];
    vector< vector<uint32_t> > nestedlevels;
    uint64_t pos, pos1, pos2, len1, len2, idx, cumlen, i, j, k;

    // Count the number of factors that have position i as the starting
    // position 
    for (i=0; i<numfacs; i++)
    {
        pos = get_field_64(positions, logrefseqlen, i); 
        // Only count positions that are not part of run length encoded
        // Ns 
        if (pos < refseqlen)
        {
            posindices[pos].push_back(i);
        }
    }

    vector<uint32_t> sortedpositions;
    // The object that has the sorting function
    SortClass sortfunc(cumseqlens.back(), compfacstarts, numfacs);
    // Sort the factors that belong to each position and put it in a
    // temporary vector
    for (i=0; i<refseqlen; i++)
    {
        if (posindices[i].size() > 0)
        {
            sort(posindices[i].begin(), posindices[i].end(), sortfunc);
            for (j=0; j<posindices[i].size(); j++)
            {
                sortedpositions.push_back(posindices[i][j]);
            }
            posindices[i].clear();
        }
    }
    delete [] posindices;

    // Create an empty vector for the first level
    nestedlevels.push_back(vector<uint32_t>());
    // Add the first factor onto the vector
    nestedlevels.back().push_back(sortedpositions.at(0));
    // Store the pos and len of the first factor
    pos1 = get_field_64(positions, logrefseqlen, sortedpositions.at(0));
    if (sortedpositions.at(0)+1 == numfacs)
        len1 = cumseqlens.back() -
               compfacstarts.select1(sortedpositions.at(0)+1);
    else
        len1 = compfacstarts.select1(sortedpositions.at(0)+2) - 
               compfacstarts.select1(sortedpositions.at(0)+1);
    i = 1;
    // Go through all factors and put them in the appropriate level
    // depending on how they are ordered
    while (i < sortedpositions.size())
    {
        pos2 = get_field_64(positions, logrefseqlen, sortedpositions.at(i));
        if (sortedpositions.at(i)+1 == numfacs)
            len2 = cumseqlens.back() - 
                   compfacstarts.select1(sortedpositions.at(i)+1);
        else
            len2 = compfacstarts.select1(sortedpositions.at(i)+2) - 
                   compfacstarts.select1(sortedpositions.at(i)+1);
        // Current factor is not within the last factor that was
        // inserted at level 0
        if (pos1 <= pos2 && (pos1+len1) <= (pos2+len2))
        {
            nestedlevels.at(0).push_back(sortedpositions.at(i));
            pos1 = pos2;
            len1 = len2;
        }
        // Current factor is within the last factor inserted at level 0
        else
        {
            j = 1;
            // Keep going up the levels until the current factor is not
            // within the last factor inserted into that level
            while (1)
            {
                // Not enough levels so make a new level
                if (nestedlevels.size() == j)
                {
                    nestedlevels.push_back(vector<uint32_t>());
                    nestedlevels.back().push_back(sortedpositions.at(i));
                    break;
                }
                else
                {
                    idx = nestedlevels.at(j).back();
                    pos1 = get_field_64(positions, logrefseqlen, idx);
                    if (idx+1 == numfacs)
                        len1 = cumseqlens.back() -
                               compfacstarts.select1(idx+1);
                    else
                        len1 = compfacstarts.select1(idx+2) - 
                               compfacstarts.select1(idx+1);
                    // Current factor is not contained by previous
                    // factor in the level so insert current factor into
                    // the level
                    if (pos1 <= pos2 && (pos1+len1) <= (pos2+len2))
                    {
                        nestedlevels.at(j).push_back(sortedpositions.at(i));
                        break;
                    }
                    // Try the next level
                    j++;
                }
            }
            // Make the previous factor be the last factor inserted into
            // level 0
            idx = nestedlevels.at(0).back();
            pos1 = get_field_64(positions, logrefseqlen, idx);
            if (idx+1 == numfacs)
                len1 = cumseqlens.back() -
                       compfacstarts.select1(idx+1);
            else
                len1 = compfacstarts.select1(idx+2) - 
                       compfacstarts.select1(idx+1);
        }
        i++;
    }

    // Allocate memory for storing the indices of the positions for each
    // nested level 
    nll = new Array(sortedpositions.size(), numfacs-1);
    levelidx = new Array(nestedlevels.size()+1, sortedpositions.size());
    numlevels = nestedlevels.size();

    k = 0;
    cumlen = 0;
    // Store the position indices and the index into the start of a
    // level
    for (i=0; i<numlevels; i++)
    {
        levelidx->setField(i,cumlen);
        cumlen += nestedlevels[i].size();
        for (j=0; j<nestedlevels[i].size(); j++,k++)
        {
            nll->setField(k, nestedlevels[i][j]);
        }
        nestedlevels[i].clear();
    }
    // Have one more than the length of the reference sequence such that
    // levelidx[j+1] - levelidx[j] = number of factors that are in that
    // level
    levelidx->setField(i,cumlen);

}
