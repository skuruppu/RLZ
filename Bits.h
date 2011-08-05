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

#include <fstream>
#include <exception>

class BitReader
{
    public:
        // Empty constructor
        BitReader();

        // Constructor takes in an input stream to read from
        BitReader(std::ifstream &infile);

        // Returns a bit from the input stream
        unsigned int read_bit();

        // Read the binary representation of an integer that takes
        // length bits
        unsigned int binary_to_int(unsigned int length);

    private:
        std::ifstream &infile;

        // Buffer that stores the bits to be read 
        unsigned char currbyte;

        // The number of bits not occupied in currbyte. Starts at 8 and
        // when remaining gets to zero, need to read another byte from
        // the input stream.
        unsigned short remaining;

        // Initialising the mask for getting the most significant bit from
        // currbyte
        const static unsigned char mask = 0x80;
};

class BitWriter
{
    public:
        // Empty constructor
        BitWriter();

        // Constructor takes in an output stream to write to
        BitWriter(std::ofstream &outfile);

        // Writes the given bit to the output stream
        void write_bit(unsigned int bit);

        // Writes the current buffered bits to output stream
        void flush();

        // Writes the given integer in binary representation using
        // length bits
        void int_to_binary(unsigned int n, unsigned int length);

    private:
        // Output stream
        std::ofstream &outfile;

        // Buffer that stores the bits to be written
        unsigned char currbyte;

        // The number of bits not occupied in currbyte. Starts at 8 and
        // when remaining gets to zero, need to write currbyte.
        unsigned short remaining;
};

class GolombCoder
{
    public:
        // The constructor for Golomb encoding
        GolombCoder(BitWriter &bwriter);

        // The constructor for Golomb encoding that allows you to set
        // the divisor
        GolombCoder(BitWriter &bwriter, unsigned int b);

        // The constructor for Golomb decoding
        GolombCoder(BitReader &breader);

        // The constructor for Golomb decoding that allows you to set
        // the divisor
        GolombCoder(BitReader &breader, unsigned int b);

        // Golomb encoding either using the default divisor or when the
        // divisor was set by the constructor
        void golomb_encode(unsigned int n);

        // Golomb encoding when a custom divisor is used
        void golomb_encode(unsigned int n, unsigned int b);

        // Golomb decoding either using the default divisor or when the
        // divisor was set by the constructor
        unsigned int golomb_decode();

        // Golomb decoding when a custom divisor is used
        unsigned int golomb_decode(unsigned int b);

    private:
        // Output stream for Golomb encoding
        BitWriter &bwriter;

        // Output stream for Golomb decoding
        BitReader &breader;

        // The divisor
        unsigned int b;

        // The number of bits for the remainder
        unsigned int log2b;

        // Checks if b (the divisor) is a power of two and returns
        // log_2(b) if b is a power of two. Otherwise returns zero.
        int is_power_of_two(unsigned int b);

        // The method that does the Golomb encoding
        void golomb_encode(unsigned int n, unsigned int b, unsigned int log2b);

        // The method that does the Golomb decoding
        unsigned int golomb_decode(unsigned int b, unsigned int log2b);
};

class BitsUnexpectedException: public std::exception
{
    public:
        virtual const char* what() const throw();
};

class BitsEOFException: public std::exception
{
    public:
        virtual const char* what() const throw();
};
