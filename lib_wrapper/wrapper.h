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

extern Library GLOBAL_LIB_TYPE;

class ArrayReference
{
    public:

        virtual ~ArrayReference();

        virtual operator uint64_t() = 0;

        virtual ArrayReference& operator=(const uint64_t) = 0;
};

class CDSArrayReference : public ArrayReference
{
    public:

        CDSArrayReference(cds_utils::Array*, uint64_t);

        ~CDSArrayReference();

        operator uint64_t();

        CDSArrayReference& operator=(uint64_t);

    private:

        cds_utils::Array *cdsarray;

        uint64_t idx;
};

class Array
{
    public:

        static Array* create(const uint64_t, const uint64_t);

        static Array* create(ifstream&);

        virtual ~Array();

        virtual uint64_t getField(const uint64_t) = 0;

        virtual uint64_t setField(const uint64_t, const uint64_t) = 0;

        virtual void save(ofstream&) = 0;

        virtual size_t getSize() = 0;

        virtual uint64_t getLength() = 0;

        virtual CDSArrayReference operator[](const uint64_t) = 0;
};

class CDSArray : public Array
{
    public:

        CDSArray();

        CDSArray(const uint64_t, const uint64_t);

        CDSArray(ifstream&);

        ~CDSArray();

        uint64_t getField(const uint64_t);

        uint64_t setField(const uint64_t, const uint64_t);

        void save(ofstream&);

        size_t getSize();

        uint64_t getLength();

        CDSArrayReference operator[](const uint64_t);

    private:

        cds_utils::Array *cdsarray;
};

}
