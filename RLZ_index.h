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

#include <stdio.h>
#include <vector>

#include <Array.h>
#include <BitSequenceSDArray.h>
#include <BitSequence.h>
#include <TextIndex.h>

typedef struct occurrence
{
    uint64_t seq;
    uint64_t pos;
} occ_t;

class RLZ_index
{
    public:
        /* Constructor */
        RLZ_index(char *filenames);

        /* Destructor */
        ~RLZ_index(); 

        /** Display function */
        void display();

        /** Count function */
        void count();

        /** Search function */
        void locate();

        /** Factor file decode function */
        void decode();

        /** Prints the space usage of the RLZ data structure */
        int size();

    private:
        /* Length of the reference sequence */
        uint64_t refseqlen;
        uint64_t logrefseqlen;

        /* Number of sequences */
        uint64_t numseqs;

        /* Number of factors */
        uint64_t numfacs;

        /* The names of the input files */
        char **filenames;

        // If this is true, can only read the data structures required
        // to implement display() query
        bool displayonly;

        /*-------------------------------------------------------------------*/
        /* Default data structures and methods                               */
        /*-------------------------------------------------------------------*/

        /* Reference sequence as a bit vector with 3bpb encoding
         * {a,c,g,t,n} */
        cds_utils::Array *refseq;

        /* Factor positions as a bit vector using logn bits per position */
        cds_utils::Array *positions;

        /* Factor facstarts in a rank-select data structure */
        cds_static::BitSequenceSDArray *facstarts;

        /* Sequence cumseqlens for numseqs sequences */
        cds_utils::Array *cumseqlens;

        /** Reads the reference sequence into memory and fills in the
         * refseqlen and refseq variables. Also creates the suffix array
         * of the reference sequence and stores it in SA.
         * @param filename Reference sequence file name
         */
        void read_reference_sequence(char *filename);

        /* Reads the factors into memory and fill the positions, facstarts
         * and cumseqlens arrays */
        void read_compressed_factors();

        /** Retrieves a substring starting at start and ending at end for
         * seq.
         * @param seq Sequence to retrieve from
         * @param start Start position to retrieve from
         * @param end End position to stop retrieving
         * @return Time taken to evaluate the query in microseconds
         */
        long display(uint64_t seq, uint64_t start, uint64_t end,
                     std::vector<uint> &substring);

        /** Retrieves a substring starting and ending at absolute
         * positions in the index.
         * @param start Start position in the index to retrieve from
         * @param end End position in the index to retrieve from
         * @param substring Vector to store the retrieved substring
         * @return Time taken to evaluate the query in microseconds
         */
        long display(uint64_t start, uint64_t end, vector <uint> &substring);

        /** Given a factor id, returns the length of that factor.
         * @param facidx Factor id
         * @return Length of the factor
         */
        inline uint64_t factor_length(uint32_t facidx);

        /*-------------------------------------------------------------------*/
        /* Locate related data structures and methods                        */
        /*-------------------------------------------------------------------*/

        // Suffix array of the reference sequence
        cds_utils::Array *sa;

        // Nested level lists of the sorted factors
        cds_utils::Array *nll;
        cds_utils::Array *levelidx;
        uint32_t numlevels;

        // Compressed bit vectors to indicate at which positions in the
        // reference sequence factors start and end at
        cds_static::BitSequenceRRR *isstart;
        cds_static::BitSequenceRRR *isend;

        cds_static::BitSequenceSDArray *seqfacstart;

        /** Implements the search functionality.
         * @param pattern Pattern to search for
         * @param ptnlen Length of the pattern being searched for
         * @param occs Vector to store occurrences if iscount=false
         * @param iscount True if just counting the occurrence, false otherwise
         * @return Number of occurrences of the pattern
         */
        uint64_t search(const char *pattern, unsigned int ptnlen, 
                        vector<occ_t>& occs, bool iscount=false);

        /** Returns the boundaries of the suffix array that contains the
         * given pattern. lb and rb are set to (uint64_t)-1 if pattern
         * does not occur in the suffix array.
         * @param pattern Pattern to search for
         * @param lb Variable to store the return left boundary
         * @param rb Variable to store the return right boundary
         */
        void sa_binary_search(cds_utils::Array &pattern, uint64_t *lb,
                              uint64_t *rb);


        /** Returns the boundaries of a level in the nested level list
         * that contains the given interval represented by the variables
         * start and end. lb and lr are set to (uint32_t)-1 if the
         * interval does not occur in the suffix array.
         * @param start Start position of the interval
         * @param end End position of the interval
         * @param lb Initial left boundary and the place to store the return left boundary
         * @param rb INitial right boundary and the place to store the return right boundary
         */
        void facs_binary_search(uint64_t start, uint64_t end, 
                                uint32_t *lb, uint32_t *rb);

        /** Returns the boundaries of a level in the nested level list 
         * that contains the factors that start with the given start
         * position. lb and lr are set to (uint32_t)-1 if the interval
         * does not occur in the suffix array.
         * @param start Start position of an interval to search for
         * @param lb Initial left boundary and the place to store the return left boundary
         * @param rb Initial right boundary and the place to store the return right boundary
         */
        void factor_start_binary_search(uint64_t start, uint32_t *lb,
                                        uint32_t *rb);

        /** Returns the boundaries of a level in the nested level list
         * that contains the factors that end with the given end
         * position. lb and lr are set to (uint32_t)-1 if the interval
         * does not occur in the suffix array.
         * @param end End position of an interval to search for
         * @param lb Initial left boundary and the place to store the return left boundary
         * @param rb Initial right boundary and the place to store the return right boundary
         */
        void factor_end_binary_search(uint64_t end, uint32_t *lb,
                                      uint32_t *rb);

        /** Compare the given substring to the reference sequence
         * starting from the start position.
         * @param substr Substring to search for
         * @param start Position in the reference sequence to search from
         * @param len Length of the substring to search for
         * @return Returns true if the substrings are equal and false otherwise
         */
        inline bool compare_substr_to_refseq(cds_utils::Array& substr, 
                                             uint64_t start, 
                                             uint64_t len);
};
