#include <iostream>
#include <fstream>
#include <stdint.h>

#include "util.h"

bool read_sequence(const char* filename, char** sequence, uint64_t* length)
{
    std::ifstream file;
    file.open(filename, std::ifstream::in);

    // File could not be opened
    if (!file)
    {
        return false;
    }

    // Get length of file
    file.seekg(0, file.end);
    *length = file.tellg();
    file.seekg (0, file.beg);

    *sequence = new char[(*length) + 1];

    // Read whole file into sequence
    file.read(*sequence, *length);
    (*sequence)[*length] = '\0';

    file.close();

    return true;
}
