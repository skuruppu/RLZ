/* RLZ_index display()
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
 * RLZ_index - Relative Lempel Ziv
 * Implements search for a compressed RLZ index. So far only supports
 * random access (display()).
 * Authors: Simon Puglisi (simon.puglisi@rmit.edu.au)
 *          Shanika Kuruppu (kuruppu@csse.unimelb.edu.au)
 */

#include <stdio.h>
#include <vector>

#include <Array.h>
#include <BitSequenceSDArray.h>
#include <BitSequence.h>
#include <Mapper.h>
#include <Sequence.h>

class RLZ_index
{
    public:
        /* Constructor */
        RLZ_index(char *filenames);

        /* Destructor */
        ~RLZ_index(); 

        /* Display function */
        void display();

        /* Search function */
        void search();

        /* Factor file decode function */
        void decode();

        /* Prints the space usage of the RLZ data structure */
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

        /*-------------------------------------------------------------------*/
        /* Default data structures                                           */
        /*-------------------------------------------------------------------*/

        /* Reference sequence as a bit vector with 3bpb encoding
         * {a,c,g,t,n} */
        cds_utils::Array *refseq;

        /* Factor positions as a bit vector using logn bits per position */
        cds_utils::Array *positions;
        //cds_static::WaveletTreeNoptrs *positions;

        /* Factor facstarts in a rank-select data structure */
        cds_static::BitSequenceSDArray *facstarts;

        /* Sequence cumseqlens for numseqs sequences */
        uint64_t *cumseqlens;

        /* Reads the reference sequence into memory and fills in the
         * refseqlen and refseq variables. Also creates the suffix array
         * of the reference sequence and stores it in SA. */
        void read_reference_sequence(char *filename);

        /* Reads the suffix array of the refrence sequence into memory */
        //void read_suffix_array(char *filename);

        /* Reads the factors into memory and fill the positions, facstarts
         * and cumseqlens arrays */
        void read_compressed_factors();

        void decompress_factors(std::ifstream &facfile,
                                std::vector<uint64_t> &poss,
                                std::vector<uint64_t> &lens,
                                uint64_t &totseqlen);

        /* Retrieves a substring starting at start and ending at end for
         * seq. Returns the total number of microsecs taken to execute
         * the query */
        long display(uint64_t seq, uint64_t start, uint64_t end,
                     std::vector<uint> &substring);

        /*-------------------------------------------------------------------*/
        /* Locate data structures                                            */
        /*-------------------------------------------------------------------*/

        // Suffix array of the reference sequence
        cds_utils::Array *sa;

        void sa_binary_search(cds_utils::Array &pattern, uint64_t *cl,
                              uint64_t *cr);
};
