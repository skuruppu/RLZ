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
 * Wrapper for using various data structure libraries with RLZ.
 * Authors: Shanika Kuruppu
 */

#include <cstring>

#include "Array.h"
#include "BitString.h"
#include "BitSequenceSDArray.h"
#include "BitSequenceRRR.h"
#include "libcdsTrees.h"

namespace lib_wrapper
{

// Libraries supported by wrapper
enum Library
{
    LIBCDS,
    LIBSDSL
};

class CDSArrayReference
{
    public:

        CDSArrayReference(cds_utils::Array* cdsarray, const uint64_t idx) :
            cdsarray(cdsarray), idx(idx) {}

        ~CDSArrayReference() {}

        operator uint64_t()
        {
            return (*cdsarray)[idx];
        }

        CDSArrayReference& operator=(const uint64_t val)
        {
            cdsarray->setField(idx, val);
            return *this;
        }

    private:

        cds_utils::Array *cdsarray;

        uint64_t idx;
};

class CDSArray
{
    public:

        CDSArray() : cdsarray(NULL) {}

        CDSArray(const uint64_t size, const uint64_t maxval)
        {
            cdsarray = new cds_utils::Array(size, maxval);
        }

        CDSArray(ifstream& file)
        {
            cdsarray = new cds_utils::Array(file);
        }

        ~CDSArray()
        {
            delete cdsarray;
        }

        void save(ofstream& file)
        {
            cdsarray->save(file);
        }

        size_t getSize()
        {
            return cdsarray->getSize();
        }

        uint64_t getLength()
        {
            return cdsarray->getLength();
        }

        CDSArrayReference operator[](const uint64_t idx)
        {
            return CDSArrayReference(cdsarray, idx);
        }

    private:

        cds_utils::Array *cdsarray;

        // Disable copy constructor and assignment operator
        CDSArray(const CDSArray& other);

        CDSArray& operator=(const CDSArray& other);
};

class CDSBitString
{
    public:

        CDSBitString() : cdsbitstring(NULL) {}

        CDSBitString(const uint64_t length)
        {
            cdsbitstring = new cds_utils::BitString(length);
        }

        ~CDSBitString()
        {
            delete cdsbitstring;
        }

        void setBit(const uint64_t idx)
        {
            cdsbitstring->setBit(idx);
        }

    private:

        cds_utils::BitString *cdsbitstring;

        // Disable copy constructor and assignment operator
        CDSBitString(const CDSBitString& other);

        CDSBitString& operator=(const CDSBitString& other);

        friend class CDSBitSequenceSDArray;
        friend class CDSBitSequenceRRR;
};

class CDSBitSequenceSDArray
{
    public:

        CDSBitSequenceSDArray() : cdsbitsequencesdarray(NULL) {}

        CDSBitSequenceSDArray(CDSBitString& bitstring)
        {
            cdsbitsequencesdarray =
                new cds_static::BitSequenceSDArray(*(bitstring.cdsbitstring));
        }

        CDSBitSequenceSDArray(ifstream& file)
        {
            cdsbitsequencesdarray = cds_static::BitSequenceSDArray::load(file);
        }

        ~CDSBitSequenceSDArray()
        {
            delete cdsbitsequencesdarray;
        }

        uint64_t select1(const uint64_t idx)
        {
            return cdsbitsequencesdarray->select1(idx);
        }

        uint64_t rank1(const uint64_t idx)
        {
            return cdsbitsequencesdarray->rank1(idx);
        }

        bool operator[](const uint64_t idx)
        {
            return cdsbitsequencesdarray->access(idx);
        }

        void save(ofstream& file)
        {
            cdsbitsequencesdarray->save(file);
        }

        uint64_t getSize()
        {
            return cdsbitsequencesdarray->getSize();
        }

    private:

        cds_static::BitSequenceSDArray *cdsbitsequencesdarray;

        // Disable copy constructor and assignment operator
        CDSBitSequenceSDArray(const CDSBitSequenceSDArray& other);

        CDSBitSequenceSDArray& operator=(const CDSBitSequenceSDArray& other);
};

class CDSBitSequenceRRR
{
    public:

        CDSBitSequenceRRR() : cdsbitsequencerrr(NULL) {}

        CDSBitSequenceRRR(CDSBitString& bitstring)
        {
            cdsbitsequencerrr =
                new cds_static::BitSequenceRRR(*(bitstring.cdsbitstring));
        }

        CDSBitSequenceRRR(ifstream& file)
        {
            cdsbitsequencerrr = cds_static::BitSequenceRRR::load(file);
        }

        ~CDSBitSequenceRRR()
        {
            delete cdsbitsequencerrr;
        }

        uint64_t select1(const uint64_t idx)
        {
            return cdsbitsequencerrr->select1(idx);
        }

        uint64_t rank1(const uint64_t idx)
        {
            return cdsbitsequencerrr->rank1(idx);
        }

        bool operator[](const uint64_t idx)
        {
            return cdsbitsequencerrr->access(idx);
        }

        void save(ofstream& file)
        {
            cdsbitsequencerrr->save(file);
        }

        uint64_t getSize()
        {
            return cdsbitsequencerrr->getSize();
        }

    private:

        cds_static::BitSequenceRRR *cdsbitsequencerrr;

        // Disable copy constructor and assignment operator
        CDSBitSequenceRRR(const CDSBitSequenceRRR& other);

        CDSBitSequenceRRR& operator=(const CDSBitSequenceRRR& other);
};

class CDSIntVector
{
    public:

        CDSIntVector() : array(NULL), width(0), length(1024) {}

        CDSIntVector(uint64_t width)
        {
            this->width = width;
            this->length = 10;
            this->array = new unsigned int[length];
        }

        ~CDSIntVector()
        {
            delete [] array;
        }

        void setField(uint64_t idx, uint64_t val)
        {
            // Double memory if we run out of space
            if (idx*width/(sizeof(unsigned int)*8)+1 >= length)
            {
                unsigned int *newarray = new unsigned int[length*2];
                memcpy(newarray, array, length*sizeof(unsigned int));
                delete [] array;
                array = newarray;
                length *= 2;
            }

            set_field_64(array, width, idx, val);
        }

        uint64_t operator[](const uint64_t idx)
        {
            return get_field_64(array, width, idx);
        }

    private:

        unsigned int* array;

        uint64_t width;

        uint64_t length;

        // Disable copy constructor and assignment operator
        CDSIntVector(const CDSIntVector& other);

        CDSIntVector& operator=(const CDSIntVector& other);
};

}
