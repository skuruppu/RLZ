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

#include "Array.h"

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
            return this->cdsarray->getField(this->idx);
        }

        CDSArrayReference& operator=(const uint64_t val)
        {
            this->cdsarray->setField(this->idx, val);
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
            this->cdsarray = new cds_utils::Array(size, maxval);
        }

        CDSArray(ifstream& file)
        {
            this->cdsarray = new cds_utils::Array(file);
        }

        ~CDSArray()
        {
            delete this->cdsarray;
        }

        uint64_t getField(const uint64_t idx)
        {
            return this->cdsarray->getField(idx);
        }

        uint64_t setField(const uint64_t idx, const uint64_t val)
        {
            return this->cdsarray->setField(idx, val);
        }

        void save(ofstream& file)
        {
            this->cdsarray->save(file);
        }

        size_t getSize()
        {
            return this->cdsarray->getSize();
        }

        uint64_t getLength()
        {
            return this->cdsarray->getLength();
        }

        CDSArrayReference operator[](const uint64_t idx)
        {
            return CDSArrayReference(this->cdsarray, idx);
        }

    private:

        cds_utils::Array *cdsarray;

        // Disable copy constructor and assignment operator
        CDSArray(const CDSArray& other);

        CDSArray& operator=(const CDSArray& other);
};

}
