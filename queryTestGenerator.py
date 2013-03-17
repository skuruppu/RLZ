#! /usr/bin/env python

# queryTestGenerator.py generates tests for the display, search and locate
# queries.
#
# Author: Shanika Kuruppu
# Email: shanika.kuruppu@gmail.com
# Date: 17/03/2013

import sys, optparse, random

# Read in the input sequences and their lengths
def readSequences(files):
    sequences = []
    sequencelens = []
    for filename in files:
        file = open(filename)
        sequence = file.read()
        file.close()

        sequences.append(sequence)
        sequencelens.append(len(sequence))

    return (sequences, sequencelens)

def generateDisplayTests(sequences, sequencelens, numqueries, querylen,
                         testfile, expoutfile):
    numseqs = len(sequences)
    i = 0
    while i < numqueries:
        # Get the sequence number
        s = random.randint(0, numseqs - 1)

        if sequencelens[s] < querylen:
            continue

        # Set the boundaries
        min = 0
        max = sequencelens[s]

        # Get the left boundary for the query
        x = random.randint(min, max - querylen)
        y = x + querylen

        # Write to file
        testfile.write("%d %d %d\n" % (s, x, y))
        expoutfile.write(sequences[s][x:x + querylen] + '\n')

        i += 1

def generatePattern(querylen):
    alpha = ['a','c','g','t']

    pattern = ''
    for j in range(0, querylen):
        pattern += alpha[random.randint(0, len(alpha) - 1)]
    return pattern

def generateCountTests(sequences, sequencelens, numqueries, querylen,
                       testfile, expoutfile):

    for i in range(0, numqueries):
        pattern = generatePattern(querylen)
        testfile.write("%s\n" % (pattern))

        expocc = 0
        for sequence in sequences:
            start = 0
            while start < len(sequence):
                pos = sequence.find(pattern, start)
                if pos != -1:
                    expocc += 1
                    start = pos+1
                else:
                    break

        expoutfile.write("%s : %d\n" % (pattern, expocc))

# Main function
def main(args=None):
    usage = """%prog query queries querylen prefix ref file ...
    query: Query to execute ('d' for display, 'c' for count, 'l' for locate)
    queries: Number of queries to output
    querylen: Length of the query
    out_file_prefix: Prefix of output files
    ref: Reference file name
    file ...: Input file names"""

    if args is None:
        args = sys.argv[1:]

    # Parse arguments and options
    parser = optparse.OptionParser(usage)
    (options, args) = parser.parse_args(args)

    # Check if there's at least four arguments
    if len(args) < 6:
        parser.print_help()
        parser.exit()

    # Get query type
    qrytype = args[0]
    if qrytype not in ['d', 'c', 'l']:
        parser.print_help()
        parser.exit()

    # Get number of queries to output
    numqueries = int(args[1])

    # Length of queries
    querylen = int(args[2])

    # Output file prefix
    prefix = args[3]

    # Get file names where the first file is the reference
    files = args[4:]

    # Read in the contents of the files
    sequences, sequencelens = readSequences(files)

    # Open output files
    testfile = open(prefix + '.test', 'w')
    expoutfile = open(prefix + '.exp', 'w')

    # Seed the random number generator with the output file prefix
    random.seed(hash(prefix))

    # Generate tests for display query
    if qrytype == 'd':
        generateDisplayTests(sequences, sequencelens, numqueries, querylen,
                             testfile, expoutfile)
    elif qrytype == 'c':
        generateCountTests(sequences, sequencelens, numqueries, querylen,
                           testfile, expoutfile)
        pass
    elif qrytype == 'l':
        pass

    testfile.close()
    expoutfile.close()

if __name__ == "__main__":
    sys.exit(main())
