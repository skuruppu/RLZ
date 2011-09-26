/* RLZ index 
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
 * Implements a self-index using RLZ.
 * Authors: Shanika Kuruppu (kuruppu@csse.unimelb.edu.au)
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
#include "alphabet.h"

#define BPB 3
#define INVALID 0xffffffff
#define BUFSIZE 1000

using namespace std;
using namespace cds_utils;
using namespace cds_static;

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

    initialise_nucl_converters();

    RLZ_index *rlzidx = new RLZ_index(argv[1]);

    rlzidx->size();
    rlzidx->locate();

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

    // Read in data structures necessary to just implement display()

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
    //idxfile.read((char*)&cumseqlens, numseqs*sizeof(uint64_t));
    for (uint64_t i=0; i<numseqs; i++)
        idxfile.read((char*)&cumseqlens[i], sizeof(uint64_t));

    // Calculate the log of the reference sequence length
    refseqlen = refseq->getLength()-1; // length includes null byte
    uint64_t i = floor(log2(refseqlen));
    logrefseqlen = ((unsigned)(1<<i) != refseqlen) ? i+1 : i;

	idxfile.get();
	if (idxfile.eof())
	{
		displayonly = true;
		idxfile.close();
		return;
	}
	idxfile.unget();

    // Read in data stuctures needed to implement locate() and count()
    // queries

    // Read in the suffix array
    sa = new Array(idxfile);

    // Read in the suffix tree
    //st = SuffixTree::load(idxfile);

    // Read the nested level list and level index
    nll = new Array(idxfile);
    idxfile.read((char*)&numlevels, sizeof(uint32_t));
    levelidx = new uint32_t[numlevels+1];
    idxfile.read((char*)levelidx, (numlevels+1)*sizeof(uint32_t));

    // Read the isstart and isend bit vectors
    isstart = BitSequenceRRR::load(idxfile);
    isend = BitSequenceRRR::load(idxfile);

    // Read the seqfacstart bit vector
    seqfacstart = BitSequenceSDArray::load(idxfile);

    idxfile.close();
}

RLZ_index::~RLZ_index()
{
    delete refseq;
    delete positions;
    delete facstarts;
    delete [] cumseqlens;
    delete sa;
    delete nll;
    delete [] levelidx;
    delete isstart;
    delete isend;
}

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
    uint64_t i, seq, start, end;
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
            cout << int_to_nucl[substring.at(i)];
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
    uint64_t i=0;
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

    return display(cumseqlens[seq]+start, cumseqlens[seq]+end, substring);
}

long RLZ_index::display(uint64_t start, uint64_t end, vector <uint> &substring)
{
    uint64_t i=0, rk, p, l, b, s;
    long msec1, msec2;
    timeval tv;

    // Start timer 
    gettimeofday(&tv, NULL);
    msec1 = tv.tv_sec*1000*1000 + tv.tv_usec;

    rk = facstarts->rank1(start);
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
        s = p+start-b;
        for (i=s; i<(p+l) && i<(s+end-start); i++)
        {
            substring.push_back(nucl_to_int['n']);
        }
    }
    else
    {
        s = p+start-b;
        for (i=s; i<(p+l) && i<(s+end-start) && i<refseqlen; i++)
        {
            substring.push_back(refseq->getField(i));
        }
    }

    // If more factors need to be decoded
    if (i == p+l)
    {
        // Get the chars from subsequent factors  
        start = b + l;
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
                s = p+start-b;
                for (i=s; i<(p+l) && i<(s+end-start); i++)
                {
                    substring.push_back(nucl_to_int['n']);
                }
            }
            else
            {
                s = p+start-b;
                for (i=s; i<(p+l) && i<(s+end-start) && i<refseqlen; i++)
                {
                    substring.push_back(refseq->getField(i));
                }
            }

            start = b + l;
        }
    }

    // End timer 
    gettimeofday(&tv, NULL);
    msec2 = tv.tv_sec*1000*1000 + tv.tv_usec;

    return (msec2 - msec1);
}

void RLZ_index::count()
{
    char *pattern = new char[1000];
    unsigned int ptnlen;
    uint64_t occurrences, totptns, totlen, totocc;
    vector<occ_t> occs(0);
    long msec1, msec2, tottime;
    timeval tv;

    tottime = totptns = totlen = totocc = 0;
    while (true)
    {
        cin.getline(pattern, 1000);
        ptnlen = strlen(pattern);

        if (ptnlen == 0) break;

        // Start timer 
        gettimeofday(&tv, NULL);
        msec1 = tv.tv_sec*1000*1000 + tv.tv_usec;

        occurrences = search(pattern, ptnlen, occs, true);

        // End timer 
        gettimeofday(&tv, NULL);
        msec2 = tv.tv_sec*1000*1000 + tv.tv_usec;

        // Update statistics
        tottime += (msec2-msec1);
        totptns ++;
        totlen += ptnlen;
        totocc += occurrences;

        cout << pattern << " : " << occurrences << endl;
    }

    cerr << (float)tottime/totocc << " microseconds/occurrences\n";

    delete [] pattern;
}

void RLZ_index::locate()
{
    char *pattern = new char[1000];
    unsigned int ptnlen;
    uint64_t occurrences, totptns, totlen, totocc;
    vector<occ_t> occs;
    long msec1, msec2, tottime;
    timeval tv;

    tottime = totptns = totlen = totocc = 0;

    while (true)
    {
        cin.getline(pattern, 1000);
        ptnlen = strlen(pattern);

        if (ptnlen == 0) break;

        // Start timer 
        gettimeofday(&tv, NULL);
        msec1 = tv.tv_sec*1000*1000 + tv.tv_usec;

        occurrences = search(pattern, ptnlen, occs);

        // End timer 
        gettimeofday(&tv, NULL);
        msec2 = tv.tv_sec*1000*1000 + tv.tv_usec;

        // Update statistics
        tottime += (msec2-msec1);
        totptns ++;
        totlen += ptnlen;
        totocc += occurrences;

        for (uint64_t i=0; i<occurrences; i++)
            cout << occs.at(i).seq << ' ' << occs.at(i).pos << endl;

        occs.clear();
    }

    cerr << (float)tottime/totocc << " microseconds/occurrences\n";

    delete [] pattern;
}

uint64_t RLZ_index::search(const char *pattern, unsigned int ptnlen,
                           vector<occ_t>& occs, bool iscount)
{
    uint64_t lb, rb, pfxlen, suflen, pos, occurrences, seq;
    uint64_t i, j, k, l;
    uint32_t poslb, posrb, facidx;
    uint64_t prevpos, prevlen, nextpos, nextlen, abspos, occpos;
    vector <uint> substr;
    occ_t occ;

    // Convert the pattern to use 3bpb
    Array intpattern(ptnlen, NUCLALPHASIZE);
    for (i=0; i<ptnlen; i++)
        intpattern.setField(i, nucl_to_int[(int)pattern[i]]);

    occurrences = 0;
    // Search for patterns occurring within factors
    // First get the positions at which the pattern occurs in the
    // reference sequence
    sa_binary_search(intpattern, &lb, &rb);
    //st_search((unsigned char*)pattern, ptnlen, lb, rb);
    if (lb != (uint64_t)-1 && rb != (uint64_t)-1)
    {
        // Add the reference sequence occurrences to the number of
        // occurrences
        occurrences += (rb - lb + 1);

        // Store the seq,pos pairs at which the pattern occurs
        if (!iscount)
        {
            for (i=lb; i<=rb; i++)
            {
                occ.seq = 0;
                occ.pos = sa->getField(i);
                //occ.pos = st->Locate(i,i);
                occs.push_back(occ);
            }
        }

        // For each position of occurrence in the 
        for (i=lb; i<=rb; i++)
        {
            // Look for factors that contain this interval in all levels
            // of the nll
            pos = sa->getField(i);
            //pos = st->Locate(i,i);
            for (j=0; j<numlevels; j++)
            {
                poslb = levelidx[j]; posrb = levelidx[j+1] - 1;
                facs_binary_search(pos, pos+ptnlen, &poslb, &posrb);
                if (poslb == (uint32_t)-1 || posrb == (uint32_t)-1)
                    continue;

                occurrences += (posrb - poslb + 1);

                // Store the seq,pos pairs at which the pattern occurs
                if (!iscount)
                {
                    for (k=poslb; k<=posrb; k++)
                    {
                        facidx = nll->getField(k);
                        seq = seqfacstart->rank1(facidx);
                        abspos = facstarts->select1(facidx+1)
                                 + (pos-positions->getField(facidx));
                        occpos = abspos - cumseqlens[seq];
                        occ.seq = seq; occ.pos = occpos;
                        occs.push_back(occ);
                    }
                }
            }
        }
    }

    // Split the pattern into two and binary search for factors starting
    // with the suffix then look for the complete occurrences
    for (i=1; i<=ptnlen/2; i++)
    {
        pfxlen = i;
        suflen = ptnlen-i;

        // Copy the 3bpb version of the pattern suffix
        Array intsufptn(suflen, NUCLALPHASIZE);
        for (j=i; j<ptnlen; j++)
            intsufptn.setField(j-i, nucl_to_int[(int)pattern[j]]);

        // Copy the 3bpb version of the pattern prefix
        Array intpfxptn(pfxlen, NUCLALPHASIZE);
        for (j=0; j<pfxlen; j++)
            intpfxptn.setField(j, nucl_to_int[(int)pattern[j]]);

        // Search for the positions in the reference sequence at which
        // the current suffix occurs
        sa_binary_search(intsufptn, &lb, &rb);
        //st->SLink(lb, rb, &lb, &rb);

        // The suffix doesn't occur in the reference sequence
        if (lb == (uint64_t)-1 || rb == (uint64_t)-1)
            continue;

        // Search for pairs of factors where the first factor ends with
        // the prefix and the second factor starts with the suffix
        for (l=lb; l<=rb; l++)
        {
            pos = sa->getField(l);
            //pos = st->Locate(l, l);
            // Ignore start positions at which factors don't start
            if (!isstart->access(pos))
                continue;
            for (j=0; j<numlevels; j++)
            {
                poslb = levelidx[j]; posrb = levelidx[j+1] - 1;
                factor_start_binary_search(pos, &poslb, &posrb);
                if (poslb == (uint32_t)-1 || posrb == (uint32_t)-1)
                    continue;

                for (k=poslb; k<=posrb; k++)
                {
                    // Get the position at which the previous factor
                    // occurs in the positions array
                    facidx = nll->getField(k)-1;
                    // Get the position component of the previous factor
                    prevpos = positions->getField(facidx);
                    // Ignore factors that are all Ns
                    if (prevpos == refseqlen) continue;

                    prevlen = factor_length(facidx);
                    // Ignore previous factor since it's too short to
                    // contain the prefix we're looking for
                    // TODO: Modify algorithm to look for patterns
                    // occurring across multiple factors
                    if (prevlen < pfxlen) continue;

                    // Compare the prefix with the equivalent length
                    // suffix of the previous factor and if they are
                    // equal then we have a match
                    if (compare_substr_to_refseq(intpfxptn,
                        prevpos+prevlen-pfxlen, pfxlen))
                    {
                        occurrences ++; 

                        // Store the seq,pos pairs at which the pattern
                        // occurs
                        if (!iscount)
                        {
                            seq = seqfacstart->rank1(facidx);
                            occpos = facstarts->select1(facidx+2)
                                     - pfxlen - cumseqlens[seq];
                            occ.seq = seq; occ.pos = occpos;
                            occs.push_back(occ);
                        }
                    }
                }
            }
        }
    }

    //st_search((unsigned char*)pattern, i-1, lb, rb);
    //uint64_t depth = st->SDepth(lb, rb);
    //uint64_t saval = st->Locate(lb, lb);

    // Split the pattern into two and binary search for factors ending
    // with the prefix then look for the complete occurrences
    for (; i<ptnlen; i++)
    {
        pfxlen = i;
        suflen = ptnlen-i;

        // Copy the 3bpb version of the prefix
        Array intpfxptn(i, NUCLALPHASIZE);
        for (j=0; j<i; j++)
            intpfxptn.setField(j, nucl_to_int[(int)pattern[j]]);

        // Copy the 3bpb version of the suffix
        Array intsufptn(ptnlen-i, NUCLALPHASIZE);
        for (; j<ptnlen; j++)
            intsufptn.setField(j-i, nucl_to_int[(int)pattern[j]]);

        // Search for the positions in the reference sequence at which
        // the current prefix occurs
        sa_binary_search(intpfxptn, &lb, &rb);
        /*
        // The previous suffix tree branch taken covers more than
        // one symbol
        if (i-1 < depth)
        {
            // Couldn't extend current match
            if (pattern[i-1] != int_to_nucl[refseq->getField(saval+i-1)])
                break;
        }
        // Need to traverse to a new child
        else
        {
            st->Child(lb, rb, (unsigned char)pattern[i-1], &lb, &rb);
            // Couldn't extend current match
            if (lb == (uint64_t)(-1))
                break;
            // Set the suffix array boundaries to narrow the search
            // for the next symbol
            else
            {
                depth = st->SDepth(lb, rb);
                if (i+1 < depth) saval = st->Locate(lb, lb);
            }
        }
        */

        // The prefix doesn't occur in the reference sequence
        if (lb == (uint64_t)-1 || rb == (uint64_t)-1)
            break;

        // Search for pairs of factors where the first factor ends with
        // the prefix and the second factor starts with the suffix
        for (l=lb; l<=rb; l++)
        {
            pos = sa->getField(l);
            //pos = st->Locate(l, l);
            // Ignore end positions at which factors don't end
            if (!isend->access(pos+pfxlen))
                continue;
            for (j=0; j<numlevels; j++)
            {
                poslb = levelidx[j]; posrb = levelidx[j+1] - 1;
                factor_end_binary_search(pos+pfxlen, &poslb, &posrb);
                if (poslb == (uint32_t)-1 || posrb == (uint32_t)-1)
                    continue;

                for (k=poslb; k<=posrb; k++)
                {
                    // Get the position at which the next factor
                    // occurs in the positions array
                    facidx = nll->getField(k)+1;
                    // Get the position component of the next factor
                    nextpos = positions->getField(facidx);
                    // Ignore factors that are all Ns
                    if (nextpos == refseqlen) continue;

                    nextlen = factor_length(facidx);
                    // Ignore next factor since it's too short to
                    // contain the suffix we're looking for
                    // TODO: Modify algorithm to look for patterns
                    // occurring across multiple factors
                    if (nextlen < suflen) continue;

                    // Compare the suffix with the equivalent length
                    // prefix of the next factor and if they are
                    // equal then we have a match
                    if (compare_substr_to_refseq(intsufptn, nextpos,
                        suflen))
                    {
                        occurrences ++; 

                        // Store the seq,pos pairs at which the pattern
                        // occurs
                        if (!iscount)
                        {
                            seq = seqfacstart->rank1(facidx-1);
                            occpos = facstarts->select1(facidx+1)
                                     - pfxlen - cumseqlens[seq];
                            occ.seq = seq; occ.pos = occpos;
                            occs.push_back(occ);
                        }
                    }
                }
            }
        }
    }

    return occurrences;
}

void RLZ_index::sa_binary_search(Array &pattern, uint64_t *lb,
                                 uint64_t *cr)
{
    uint64_t low, high, mid, pl, pr; 
    int midval, midvalleft, midvalright, c;

    pl = 0; pr = refseqlen;

    uint64_t length = pattern.getLength();
    for (uint64_t i=0; i<length; i++)
    {
        c = pattern.getField(i);
        low = pl; high = pr;

        // Binary search left
        while (low <= high)
        {
            mid = (low + high) >> 1;

            midval = refseq->getField(sa->getField(mid)+i);
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
                    *lb = mid;
                    break;
                }
                midvalleft = refseq->getField(sa->getField(mid-1)+i);
                // Discard mid and values to the right of mid
                if(midvalleft == midval)
                    high = mid - 1;
                // Left-most occurrence found
                else
                {
                    *lb = mid;
                    break;
                }
            }
        }

        // Key not found so return not found symbols
        if (low > high)
        {
            *lb = (uint64_t)(-1);
            *cr = (uint64_t)(-1);
            return;
        }

        // Binary search right
        low = *lb; high = pr;
        while (low <= high)
        {
            mid = (low + high) >> 1;


            midval = refseq->getField(sa->getField(mid)+i);
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
                midvalright = refseq->getField(sa->getField(mid+1)+i);
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
            *lb = (uint64_t)(-1);
            *cr = (uint64_t)(-1);
            return;
        }

        pl = *lb;
        pr = *cr;
    }

    return;
}

void RLZ_index::facs_binary_search(uint64_t start, uint64_t end, 
                                   uint32_t *lb, uint32_t *rb)
{
    uint32_t low, high, facidx, mid;
    uint64_t pos, len;

    // Binary search left
    low = *lb; high = *rb; 
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);
        len = factor_length(facidx);

        // The factor is to the right of the current middle 
        if (end > pos+len)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (start < pos)
            high = mid-1;
        // The middle factor contains start 
        else
        { 
            // It's the left most so return mid 
            if(mid == *lb)
            {
                *lb = mid;
                break;
            }

            // Get the factor at the left of the middle
            facidx = nll->getField(mid-1);
            pos = positions->getField(facidx);
            len = factor_length(facidx);

            // mid - 1 factor is less than the current position so we've
            // reached the left most boundary 
            if (end > pos+len)
            {
                *lb = mid;
                break;
            }
            // Move the right boundary to the left of the middle factor 
            // since it matches 
            else
                high = mid-1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }

    // Binary search right
    low = *lb; high = *rb;
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);
        len = factor_length(facidx);

        // The factor is to the right of the current middle 
        if (end > pos+len)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (start < pos)
            high = mid-1;
        // The middle factor contains curridx 
        else
        { 
            // It's the right most so return mid 
            if(mid == *rb)
            {
                *rb = mid;
                break;
            }

            // Get the factor at the right of the middle
            facidx = nll->getField(mid+1);
            pos = positions->getField(facidx);
            len = factor_length(facidx);

            // mid + 1 factor is greater than the current position so we've
            // reached the right most boundary 
            if (start < pos)
            {
                *rb = mid;
                break;
            }
            // Move the left boundary to the right of the middle factor 
            // since it matches 
            else
                low = mid+1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }
}

void RLZ_index::factor_start_binary_search(uint64_t start, uint32_t *lb,
                                           uint32_t *rb)
{
    uint32_t low, high, facidx, mid;
    uint64_t pos;

    // Binary search left
    low = *lb; high = *rb;
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);

        // The factor is to the right of the current middle 
        if (start > pos)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (start < pos)
            high = mid-1;
        // The middle factor contains start 
        else
        { 
            // It's the left most so return mid 
            if(mid == *lb)
            {
                *lb = mid;
                break;
            }

            // Get the nucleotide at the left of the middle suffix +
            // offset 
            facidx = nll->getField(mid-1);
            pos = positions->getField(facidx);

            // mid - 1 factor is less than the current position so we've
            // reached the left most boundary 
            if (start > pos)
            {
                *lb = mid;
                break;
            }
            // Move the right boundary to the left of the middle factor 
            // since it matches 
            else
                high = mid-1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }

    // Binary search right
    low = *lb; high = *rb;
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);

        // The factor is to the right of the current middle 
        if (start > pos)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (start < pos)
            high = mid-1;
        // The middle factor contains start 
        else
        { 
            // It's the right most so return mid 
            if(mid == *rb)
            {
                *rb = mid;
                break;
            }

            // Get the factor at the right of the middle
            facidx = nll->getField(mid+1);
            pos = positions->getField(facidx);

            // mid + 1 factor is greater than the current position so we've
            // reached the right most boundary 
            if (start < pos)
            {
                *rb = mid;
                break;
            }
            // Move the left boundary to the right of the middle factor 
            // since it matches 
            else
                low = mid+1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }
}

void RLZ_index::factor_end_binary_search(uint64_t end, uint32_t *lb,
                                           uint32_t *rb)
{
    uint32_t low, high, facidx, mid;
    uint64_t pos, len;

    // Binary search left
    low = *lb; high = *rb;
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);
        len = factor_length(facidx);

        // The factor is to the right of the current middle 
        if (end > pos+len)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (end < pos+len)
            high = mid-1;
        // The middle factor contains end 
        else
        { 
            // It's the left most so return mid 
            if(mid == *lb)
            {
                *lb = mid;
                break;
            }

            // Get the nucleotide at the left of the middle suffix +
            // offset 
            facidx = nll->getField(mid-1);
            pos = positions->getField(facidx);
            len = factor_length(facidx);

            // mid - 1 factor is less than the current position so we've
            // reached the left most boundary 
            if (end > pos+len)
            {
                *lb = mid;
                break;
            }
            // Move the right boundary to the left of the middle factor 
            // since it matches 
            else
                high = mid-1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }

    // Binary search right
    low = *lb; high = *rb;
    while (low <= high)
    {
        // Get the middle index 
        mid = (low + high) >> 1;

        // Get the factor at the middle index 
        facidx = nll->getField(mid);
        pos = positions->getField(facidx);
        len = factor_length(facidx);

        // The factor is to the right of the current middle 
        if (end > pos+len)
            low = mid+1;
        // The factor is to the left of the current middle 
        else if (end < pos+len)
            high = mid-1;
        // The middle factor contains end 
        else
        { 
            // It's the right most so return mid 
            if(mid == *rb)
            {
                *rb = mid;
                break;
            }

            // Get the factor at the right of the middle
            facidx = nll->getField(mid+1);
            pos = positions->getField(facidx);
            len = factor_length(facidx);

            // mid + 1 factor is greater than the current position so we've
            // reached the right most boundary 
            if (end < pos+len)
            {
                *rb = mid;
                break;
            }
            // Move the left boundary to the right of the middle factor 
            // since it matches 
            else
                low = mid+1;
        } 
    }

    // Key not found 
    if (low > high)
    {
        *lb = (uint32_t)-1;
        *rb = (uint32_t)-1;
        return;
    }
}

void RLZ_index::st_search(unsigned char *pattern, unsigned int ptnlen, 
                          uint64_t& lb, uint64_t& rb)
{
    unsigned char c;
    uint64_t i, depth=0, saval=0;
    size_t pl, pr;

    st->Root(&pl, &pr);
    i = 0;
    for (i=0; i<ptnlen; i++)
    {
        c = pattern[i];
        // The previous suffix tree branch taken covers more than
        // one symbol
        if (i < depth)
        {
            // Couldn't extend current match so print factor
            if (c != int_to_nucl[refseq->getField(saval+i)])
            {
                lb = rb = (uint64_t)-1;
                return;
            }
        }
        // Need to traverse to a new child
        else
        {
            st->Child(pl, pr, c, &lb, &rb);
        
            // Couldn't extend current match so print factor
            if (lb == (uint64_t)(-1))
            {
                lb = rb = (uint64_t)-1;
                return;
            }
            // Set the suffix array boundaries to narrow the search
            // for the next symbol
            else
            {
                depth = st->SDepth(lb, rb);
                if (i+1 < depth)
                    saval = st->Locate(lb, lb);
                pl = lb; pr = rb;
            }
        }
    }
}

inline uint64_t RLZ_index::factor_length(uint32_t facidx)
{
    if (facidx+1 == numfacs)
        return cumseqlens[numseqs] - facstarts->select1(facidx+1);
        
    return facstarts->select1(facidx+2) - facstarts->select1(facidx+1);
}

inline bool RLZ_index::compare_substr_to_refseq(Array& substr,
                       uint64_t start, uint64_t len)
{
    uint64_t i, j;

    for (i=start, j=0; i<start+len; i++,j++)
    {
        if (refseq->getField(i) != substr.getField(j))
            return false;
    }

    return true;
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
    // Contents of cumseqlens
    size += numseqs*sizeof(uint64_t);
    cerr << "cumseqlens: " << numseqs*sizeof(uint64_t) << " bytes\n";
    // Size of the pointers
    size += (sizeof(cumseqlens)+sizeof(facstarts));
    // Size of the numseqs and numfacs variables
    size += (sizeof(numseqs)+sizeof(numfacs));
    totalsize += size;

    cerr << "-----\n";
    cerr << "Factors: " << size << " bytes\n";
    cerr << endl;

	if (!displayonly)
	{
		size = 0;
		// Contents of suffix array
		size += (unsigned int)sa->getSize();
		cerr << "suffix array: " << (unsigned int)sa->getSize() << " bytes\n";
		// Contents of suffix tree
		//size += (unsigned int)st->getSize();
		//cerr << "suffix tree: " << (unsigned int)st->getSize() << " bytes\n";
		// Contents of nested level lists
		size += ((unsigned int)nll->getSize() +
				 (numlevels+1)*sizeof(uint32_t));
		cerr << "nested level lists: "
			 << (unsigned int)nll->getSize() + 
				(numlevels+1)*sizeof(uint32_t)
			 << " bytes\n";
		// Contents of isstart and isend
		size += (unsigned int)isstart->getSize();
		cerr << "isstart: " << (unsigned int)isstart->getSize() << " bytes\n";
		size += (unsigned int)isend->getSize();
		cerr << "isstart: " << (unsigned int)isend->getSize() << " bytes\n";
		// Contents of seqfacstart
		size += (unsigned int)seqfacstart->getSize();
		cerr << "seqfacstart: " << (unsigned int)seqfacstart->getSize() << " bytes\n";
		// Size of the sa and nll variables
		size += (sizeof(sa)+sizeof(nll)+sizeof(levelidx)+sizeof(numlevels));
		// Size of isstart, isend and seqfacstart bit vectors
		size += (sizeof(isstart)+sizeof(isend)+sizeof(seqfacstart));
		totalsize += size;

		cerr << "-----\n";
		cerr << "Search: " << size << " bytes\n";
		cerr << endl;
	}

    cerr << "Total: " << totalsize << " bytes\n";
    cerr << endl;

    return totalsize;
}
