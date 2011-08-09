/* RLZ
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
 * Bits.cpp: Implements a bit writer that writes to an output file and a
 * bit reader that reads from an input stream of bits. Implements Golomb
 * coding and binary codig on top of the bit reader and bit writer.
 * Author: Shanika Kuruppu (kuruppu@csse.unimelb.edu.au)
 */

#include <iostream>
#include <cstdlib>

#include "Bits.h"

#define BITSPERBYTE 8
#define BITSPERWORD 32

using namespace std;

BitWriter::BitWriter() :
    outfile(*(ofstream*)NULL) { }

BitWriter::BitWriter(ofstream &outfile) :
    outfile(outfile)
{
    // Set all the bits to zeros
    this->currbyte = 0;
    // All bits in currbyte are available
    this->remaining = BITSPERBYTE;
}

void BitWriter::write_bit(unsigned int bit)
{
    if (bit > 2)
    {
        cerr << "write_bit(): Bit can only have either a 0 or 1. Input " << bit << endl;
        exit(1);
    }
   
    // Put the current bit at the back of the buffer
    currbyte <<= 1;
    currbyte |= (unsigned char)bit;
    remaining--;

    // If no space left on the buffer, flush it to the output stream
    if (remaining == 0)
    {
        flush();
    }
}

void BitWriter::flush()
{
    if (remaining < BITSPERBYTE)
    {
        // Shifts the bits to the left so that the first bit that arrived is
        // at the left-most end of the byte
        currbyte <<= remaining;

        // Write the byte
        outfile.write((char*)&currbyte, sizeof(unsigned char));

        // Re-initialize the buffer and the remaining bit counter
        currbyte = 0;
        remaining = BITSPERBYTE;
    }
}

void BitWriter::int_to_binary(uint64_t n, uint64_t length)
{
    uint64_t i;
    unsigned int bit;

    // Check if the int is possible to be represented with length bits
    if (n >= ((unsigned)0x01 << length))
    {
        cerr << "Warning: " << length << " bits is not enough to encode " 
             << n << '.' << endl;
    }

    // Write out the bit representation of n
    for (i=length; i>0; i--)
    {
        bit = n & (0x01 << (i-1));
        bit >>= (i-1);
        write_bit(bit);
    }
}

BitReader::BitReader() :
    infile(*(ifstream*)NULL) {}

BitReader::BitReader(ifstream &infile) :
    infile(infile)
{
    // Set all the bits to zeros
    this->currbyte = 0;
    // All bits in currbyte are available
    this->remaining = 0;
}

unsigned int BitReader::read_bit()
{
    BitsUnexpectedException uexp;
    BitsEOFException eofexp;

    // Check if it's EOF
    if (infile.eof())
    {
        currbyte = 0;
        remaining = 0;
        // Throw an exception to indicate that it's EOF
        throw eofexp;
    }
    // Check if the file is readable
    if (!infile.good())
    {
        currbyte = 0;
        remaining = 0;
        // Throw an exception to indicate that there's an error
        throw uexp;
    }

    // Check if currbyte has any more bits left to be read. If not then
    // read a new byte from the input file.
    if (remaining == 0)
    {
        // The infile object may have exceptions set up to be raised
        // when an IO error occurs so try to catch them if possible
        try
        {
            infile.read((char*)&currbyte, sizeof(unsigned char));
            remaining = BITSPERBYTE;
        }
        catch (ifstream::failure e)
        {
            if (infile.eof())
                throw eofexp;
            else
                throw uexp;
        }
    }

    // Store the next bit in the stream into the output variable
    unsigned int bit;
    bit = currbyte & mask;
    bit >>= BITSPERBYTE-1;

    // Shift by one bit to the left so that the next bit in the stream
    // is in front
    currbyte <<= 1;
    
    // Update remaining to the number of bits yet to be read from
    // currbyte
    remaining --;

    return bit;
}

uint64_t BitReader::binary_to_int(uint64_t length)
{
    uint64_t i, n = 0;
    unsigned int bit;

    // Read the integer
    for (i=length; i>0; i--)
    {
        bit = read_bit();
        n |= (bit << (i-1)); 
    }
    return n;
}

GolombCoder::GolombCoder(BitWriter &bwriter) :
    bwriter(bwriter), breader(*(BitReader*)NULL)
{
    // Set the default divisor
    this->b = 64;
    // Set the number of bits needed to represent the default divisor
    this->log2b = 6;
}

GolombCoder::GolombCoder(BitWriter &bwriter, unsigned int b) :
    bwriter(bwriter), breader(*(BitReader*)NULL)
{
    // Check to see if the divisor given is a power of two
    if ((this->log2b = is_power_of_two(b)) > 0)
    {
        this->b = b;
    }
    // Otherwise it's an error
    else
    {
        cerr << "The divisor b must be a power of 2." << endl;
        exit(1);
    }
}

GolombCoder::GolombCoder(BitReader &breader) :
   bwriter(*(BitWriter*)NULL), breader(breader) 
{
    // Set the default divisor
    this->b = 64;
    // Set the number of bits needed to represent the default divisor
    this->log2b = 6;
}

GolombCoder::GolombCoder(BitReader &breader, unsigned int b) :
    bwriter(*(BitWriter*)NULL), breader(breader)
{
    // Check to see if the divisor given is a power of two
    if ((this->log2b = is_power_of_two(b)) > 0)
    {
        this->b = b;
    }
    // Otherwise it's an error
    else
    {
        cerr << "The divisor b must be a power of 2." << endl;
        exit(1);
    }
}

void GolombCoder::golomb_encode(uint64_t n, unsigned int b, unsigned int log2b)
{
    uint64_t quotient, remainder, i;

    // Get the quotient and the remainder from dividing n by b
    quotient = n >> log2b;
	remainder = (n << ((sizeof(uint64_t)<<3)-log2b)) >>
			    ((sizeof(uint64_t)<<3)-log2b);
    
    // Encode the quotient
    for (i = 0; i<quotient; i++)
    {
        bwriter.write_bit(1);
    }
    bwriter.write_bit(0);

    // Encode the remainder
    bwriter.int_to_binary(remainder, log2b);

    return;
}

void GolombCoder::golomb_encode(uint64_t n, unsigned int b)
{
    unsigned int log2b;

    // Check if the divisor is large enough and its a power of two
    if (b < 2 || (log2b = is_power_of_two(b)) == 0)
    {
        cerr << "The divisor b must be a power of 2." << endl;
        exit(1);
    }

    golomb_encode(n, b, log2b);

    return;
}

void GolombCoder::golomb_encode(uint64_t n)
{
    golomb_encode(n, this->b, this->log2b);
    return;
}

int GolombCoder::is_power_of_two(unsigned int b)
{
    unsigned int i, bit;
    
    // Check that b is a power of 2
    for (i=0; i<BITSPERWORD; i++)
    {
        bit = b & 0x01;
        b >>= 1;
        if (bit)
        {
            // b is a power of 2
            if (i>0 && !b)
                return i;
            else
                return 0;
        }
    }
    return false;
}

uint64_t GolombCoder::golomb_decode(unsigned int b, unsigned int log2b)
{
    uint64_t quotient, remainder;
	unsigned int bit;

    // Decode the quotient
    quotient = 0;
    bit = breader.read_bit();
    while (bit == 1)
    {
        quotient ++;
        bit = breader.read_bit();
    }

    // Decode the remainder
    remainder = breader.binary_to_int(log2b);

    // Calculate the original value
    return quotient*b+remainder;
}

uint64_t GolombCoder::golomb_decode(unsigned int b)
{
    unsigned int log2b;

    // Check if the divisor is large enough and its a power of two
    if (b < 2 || (log2b = is_power_of_two(b)) == 0)
    {
        cerr << "The divisor b must be a power of 2." << endl;
        exit(1);
    }

    return golomb_decode(b, log2b);
}

uint64_t GolombCoder::golomb_decode()
{
    return golomb_decode(this->b, this->log2b);
}


const char* BitsUnexpectedException::what() const throw()
{
    return "Unexpected error in file.\n";
}

const char* BitsEOFException::what() const throw()
{
    return "End of file.\n";
}
