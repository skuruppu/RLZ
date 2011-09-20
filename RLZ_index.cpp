/* RLZ display()
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
 * Implements search for a compressed RLZ index. So far only supports
 * random access (display()).
 * Authors: Simon Puglisi (simon.puglisi@rmit.edu.au)
 *          Shanika Kuruppu (kuruppu@csse.unimelb.edu.au)
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>

#include <libcdsBasics.h>
#include <cppUtils.h>
#include <BitString.h>

#include "RLZ_index.h"
#include "Bits.h"

#define BPB 3
#define INVALID 0xffffffff
#define BUFSIZE 1000

using namespace std;
using namespace cds_utils;
using namespace cds_static;

uint64_t cumfaclen = 0;

int main (int argc, char **argv)
{
    char usage[] = "Usage: rlz_index FILE\n \
    FILE: Compressed index file\n";

    // Check for the correct number of arguments
    if (argc < 2)
    {
        cerr << usage;
        exit(1);
    }

    RLZ_index *rlzidx = new RLZ_index(argv[1]);

    rlzidx->size();
    rlzidx->display();

    return 0; 
}

RLZ_index::RLZ_index(char *filename) :
    // Don't want to include the reference sequence 
    numseqs(0),
    refseq(NULL),
    positions(NULL),
    facstarts(NULL),
    cumseqlens(NULL)
{
    // Open the index file
    ifstream idxfile(filename, ifstream::in);
    if (!idxfile.good())
    {
        cerr << "Error opening index file.\n";
        exit(1);
    }

    // Read the meta data
    BitReader breader = BitReader(idxfile);
    numfacs = breader.binary_to_int(sizeof(uint64_t)*8);
    numseqs = breader.binary_to_int(sizeof(uint64_t)*8);

    // Read the reference sequence
    refseq = new Array(idxfile);

    // Read the factor start positions
    facstarts = BitSequenceSDArray::load(idxfile);

    // Create a compact array to store the positions
    positions = new Array(idxfile);

    // Read the cumulative sequence lengths
    cumseqlens = new uint64_t[numseqs];
    for (uint64_t i=0; i<numseqs; i++)
        idxfile.read((char*)&cumseqlens[i], sizeof(uint64_t));

    // Calculate the log of the reference sequence length
    refseqlen = refseq->getLength()-1; // length includes null byte
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;
}

RLZ_index::~RLZ_index()
{
    delete refseq;
    delete positions;
    delete facstarts;
    delete [] cumseqlens;
}

/*
void RLZ_index::read_reference_sequence(char *filename)
{
    char *refseqfilename = new char[1024];
    ifstream refseqfile;
    uint64_t i;
    int c, val;
    
    // Initialize the reference sequence reading 
    strcpy(refseqfilename, filename);
    refseqfile.open(refseqfilename, ifstream::in);
    if (!refseqfile.good())
    {
       cerr << "Error opening file of base sequence " << filename << endl;
       exit(1);
    }

    // Get the reference sequence and its length
    refseqfile.seekg(0,ios::end);
    refseqlen = (uint64_t)refseqfile.tellg();
    int j = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<j) != refseqlen) ? j+1 : j;

    // Reset file pointer to the start of the file 
    refseqfile.seekg(0,ios::beg);

    // Initialise the reference sequence storage 
	refseq = new Array(refseqlen, 7);

    // Read each char and encode to 3bpb 
    for (i=0; i<refseqlen; i++)
    {
        c = refseqfile.get();
        switch (c)
        {
            case 'a':
                val = 0;
                break;
            case 'c':
                val = 1;
                break;
            case 'g':
                val = 2;
                break;
            case 'n':
                val = 3;
                break;
            case 't':
                val = 4;
                break;
            default:
                cerr << "RLZ_index: Unsupported character in reference sequence.\n";
                exit(1);
                break;
        }
		refseq->setField(i, val);
    }

    delete [] refseqfilename;

    return;
}

void RLZ_index::read_compressed_factors()
{
    uint64_t i, totseqlen, cumlen;
    vector<uint64_t> poss;
    vector<uint64_t> lens;
    char *facfilename = new char[1024];
    ifstream facfile;

    // Allocate memory to store the cumseqlens 
    cumseqlens = new uint64_t[numseqs+2];
    cumseqlens[0] = 0;

    // Read the factors into memory for each factor file 
    totseqlen = 0;
    numfacs = 0;
    for (i=1; i<=numseqs; i++)
    {
        strcpy(facfilename, filenames[i]);
        facfile.open(facfilename, ifstream::in);
        if (!facfile.good())
        {
            cerr << "Error opening file " << filenames[i] << "." << endl;
            exit(1);
        }
        
        cumseqlens[i] = totseqlen;
        decompress_factors(facfile, poss, lens, totseqlen);
        facfile.close();
    }
    cumseqlens[numseqs+1] = totseqlen; // The end of all the sequences

    // Create a compact array to store the positions
    positions = new Array(poss);
    
    // Store the facstarts as a compressed bit vector
    BitString lengtharray(totseqlen);
    cumlen = 0;
    for (i=0; i<numfacs; i++)
    {
        lengtharray.setBit(cumlen);
        cumlen += lens.at(i);
    }
    facstarts = new BitSequenceSDArray(lengtharray);

    delete [] facfilename;
}

void RLZ_index::decompress_factors(ifstream &facfile, vector<uint64_t> &poss, 
                           vector<uint64_t> &lens, uint64_t &totseqlen)
{
    uint64_t pos, len, cumlen;

    BitReader breader = BitReader(facfile);
    GolombCoder gdecoder = GolombCoder(breader, 64);

    cumlen = 0;
    while (1)
    {
        // Decode standard factor
        try
        {
            pos = breader.binary_to_int(logrefseqlen);
            len = gdecoder.golomb_decode();
            cumfaclen += len;
        }
        catch (BitsEOFException& eofexp)
        {
            break;
        }
        catch (BitsUnexpectedException& unexp)
        {
            cerr << unexp.what();
            exit(1);
        }

        poss.push_back(pos);
        lens.push_back(len);
        totseqlen += len;
        numfacs ++;
    }

}
*/

void RLZ_index::decode()
{
    uint64_t i, j, maxfacnum, len;
    ofstream outfile;
    char *filename = new char[1024];

    sprintf(filename, "para.dec");
    outfile.open(filename, ofstream::out);

    j = 0;
    for (i=0; i<numseqs-2; i++)
    {
        // Compressed sequence filenames start at index 1
        //sprintf(filename, "%s.dec", filenames[i+1]);
        //outfile.open(filename, ofstream::out);

        // Get the cumseqlens for the set of factors that belong to the
        // current sequence being decoded
        // Special case since facstarts->rank1(cumseqlens[i+1]-1) would result 
        // in a seg fault
        if (cumseqlens[i+1] == 0)
            maxfacnum = j + facstarts->rank1(cumseqlens[i+2]-1);
        else
            maxfacnum = j + facstarts->rank1(cumseqlens[i+2]-1) - 
                        facstarts->rank1(cumseqlens[i+1]-1);
        for (; j<maxfacnum; j++)
        {
            // Get the length of the factor
            if (j+1 == numfacs)
                len = cumseqlens[numseqs+1] - facstarts->select1(j+1);
            else
                len = facstarts->select1(j+2) - facstarts->select1(j+1);
            // Standard factor decoding
            outfile << positions->getField(j) << ' ';
            // Print the length
            outfile << len << endl;
        }
        //outfile.close();
    }
    
    outfile.close();

    delete [] filename;
}

void RLZ_index::display()
{
    uint64_t i, seq, start, end, v;
    vector <uint> substring;
    long totaltime = 0, totalchars = 0, totalqueries = 0; 

    // Read each query at a time 
    while (scanf("%lu %lu %lu", &seq, &start, &end) == 3)
    {
        // Reset the counters since we want to measure the time for
        // displaying after everything is in cache 
        if (totalqueries == 1000)
        {
            totaltime = 0;
            totalchars = 0;
            totalqueries = 0;
        }
        totaltime += display(seq, start, end, substring); 
        totalchars += (end-start);
        totalqueries ++;
        
        for (i=0; i<substring.size(); i++)
        {
            v = substring.at(i);
            switch(v)
            {
                case 1:
                    cout << 'a';
                    break;
                case 2:
                    cout << 'c';
                    break;
                case 3:
                    cout << 'g';
                    break;
                case 4:
                    cout << 'n';
                    break;
                case 5:
                    cout << 't';
                    break;
                default:
                    cerr << "RLZ_index: Unknown nucleotide bit representation." << v << "\n";
                    exit(1);
                    break;
            }
        }
        cout << '\n';
        substring.clear();
        
    }

    cerr << '\n';
    cerr << "Total factors: " << numfacs << '\n';
    cerr << "Total queries: " << totalqueries << '\n';
    cerr << "Total chars: " << totalchars << '\n';
    cerr << "Total time: " << (float)totaltime/totalchars << " microsecs/char\n";
    return;
}

long RLZ_index::display(uint64_t seq, uint64_t start, uint64_t end, vector <uint> &substring)
{
    uint64_t i=0, rk, p, l, b, s;
    long msec1, msec2;
    timeval tv;

    // Start timer 
    gettimeofday(&tv, NULL);
    msec1 = tv.tv_sec*1000*1000 + tv.tv_usec;

    // Substring to be retrieved from reference sequence
    if (seq == 0)
    {
        for (i=start; i<end; i++)
            substring.push_back(refseq->getField(i));
        
        // End timer 
        gettimeofday(&tv, NULL);
        msec2 = tv.tv_sec*1000*1000 + tv.tv_usec;

        return (msec2 - msec1);
    }

    rk = facstarts->rank1(cumseqlens[seq]+start);
    b = facstarts->select1(rk);
    if (rk == numfacs)
        l = cumseqlens[numseqs-1] - b;
    else
        l = facstarts->select1(rk+1) - b;

    // Just a standard factor
    p = positions->getField(rk-1);
    // A substring of Ns
    if (p == refseqlen)
    {
        s = p+cumseqlens[seq]+start-b;
        for (i=s; i<(p+l) && i<(s+end-start); i++)
        {
            substring.push_back(3);
        }
    }
    else
    {
        s = p+cumseqlens[seq]+start-b;
        for (i=s; i<(p+l) && i<(s+end-start) && i<refseqlen; i++)
        {
            substring.push_back(refseq->getField(i));
        }
    }

    // If more factors need to be decoded
    if (i == p+l)
    {
        // Get the chars from subsequent factors  
        start = b - cumseqlens[seq] + l;
        while (start < end)
        {
            rk ++;
            b = facstarts->select1(rk);
            if (rk == numfacs)
                l = cumseqlens[numseqs-1] - b;
            else
                l = facstarts->select1(rk+1) - b;

            // Just a standard factor
            p = positions->getField(rk-1);

            // A substring of Ns
            if (p == refseqlen)
            {
                s = p+cumseqlens[seq]+start-b;
                for (i=s; i<(p+l) && i<(s+end-start); i++)
                {
                    substring.push_back(3);
                }
            }
            else
            {
                s = p+cumseqlens[seq]+start-b;
                for (i=s; i<(p+l) && i<(s+end-start) && i<refseqlen; i++)
                {
                    substring.push_back(refseq->getField(i));
                }
            }

            start = b - cumseqlens[seq] + l;
        }
    }

    // End timer 
    gettimeofday(&tv, NULL);
    msec2 = tv.tv_sec*1000*1000 + tv.tv_usec;

    return (msec2 - msec1);
}

int RLZ_index::size()
{
    unsigned int totalsize = 0, size;

    cerr << "Reference sequence length: " << refseqlen << endl;
    cerr << "Number of factors: " << numfacs << endl;
    cerr << endl;

    // Size of refseq and refseqlen 

    size = 0;
    // Contents of the refseq array
	size += (unsigned int)refseq->getSize();
    // Size of the refseqlen variables
    size += (sizeof(refseqlen) + sizeof(logrefseqlen)); 
    totalsize += size;

    cerr << "Space consumption:\n";
    cerr << "Reference sequence: " << size << " bytes\n";

    // Size of cumseqlens, positions and facstarts arrays 

    size = 0;
    // Contents of cumseqlens array
    size += (sizeof(uint64_t)*numseqs);
    // Contents of positions array
    size += (unsigned int)positions->getSize();
    cerr << "positions: " << (unsigned int)positions->getSize() << " bytes\n";
    // Contents of facstarts data structure
    size += (facstarts->getSize());
    cerr << "facstarts: " << (unsigned int)facstarts->getSize() << " bytes\n";
    // Size of the pointers
    size += (sizeof(cumseqlens)+sizeof(facstarts));
    // Size of the numseqs and numfacs variables
    size += (sizeof(numseqs)+sizeof(numfacs));
    totalsize += size;

    cerr << "Factors: " << size << " bytes\n";
    cerr << endl;

    cerr << "Total: " << totalsize << " bytes\n";

    return totalsize;
}
