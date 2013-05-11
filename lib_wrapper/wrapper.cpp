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

#include "wrapper.h"

namespace lib_wrapper
{

Library GLOBAL_LIB_TYPE;

Array* Array::create(const uint64_t size, const uint64_t maxval)
{
    switch (GLOBAL_LIB_TYPE)
    {
        case LIBCDS:
            return new CDSArray(size, maxval);
            break;
        case LIBSDSL:
            // TODO: Import SDSL
            break;
    }

    return NULL;
}

Array* Array::create(ifstream& file)
{
    switch (GLOBAL_LIB_TYPE)
    {
        case LIBCDS:
            return new CDSArray(file);
            break;
        case LIBSDSL:
            // TODO: Import SDSL
            break;
    }

    return NULL;
}

Array::~Array() {}

CDSArray::CDSArray() : cdsarray(NULL) {}

CDSArray::CDSArray(const uint64_t size, const uint64_t maxval)
{
    this->cdsarray = new cds_utils::Array(size, maxval);
}

CDSArray::CDSArray(ifstream& file)
{
    this->cdsarray = new cds_utils::Array(file);
}

CDSArray::~CDSArray()
{
    delete this->cdsarray;
}

uint64_t CDSArray::getField(const uint64_t idx)
{
    return this->cdsarray->getField(idx);
}

uint64_t CDSArray::setField(const uint64_t idx, const uint64_t val)
{
    return this->cdsarray->setField(idx, val);
}

void CDSArray::save(ofstream& file)
{
    this->cdsarray->save(file);
}

size_t CDSArray::getSize()
{
    return this->cdsarray->getSize();
}

uint64_t CDSArray::getLength()
{
    return this->cdsarray->getLength();
}

ArrayReference::~ArrayReference() {}

CDSArrayReference CDSArray::operator[](const uint64_t idx)
{
    return CDSArrayReference(this->cdsarray, idx);
}

CDSArrayReference::CDSArrayReference(cds_utils::Array* cdsarray,
                                     const uint64_t idx) :
                                     cdsarray(cdsarray), idx(idx) {}

CDSArrayReference::~CDSArrayReference() {}


CDSArrayReference::operator uint64_t()
{
    return this->cdsarray->getField(this->idx);
}

CDSArrayReference& CDSArrayReference::operator=(const uint64_t val)
{
    this->cdsarray->setField(this->idx, val);
    return *this;
}

}
