#!/usr/bin/env sh

# runTests.sh runs tests to check if RLZ and RLZindex are operating as
# expected.
#
# Author: Shanika Kuruppu
# Email: shanika.kuruppu@gmail.com
# Date: 17/03/2013

TESTDIR="`pwd`/test"
FILES="`cat $TESTDIR/files`"

# Build the executables
make clobber
make rlz rlzindex
echo ""

# Test if decompressed output is the same as the input
testDecompress()
{
    I=0
    for FILE in $FILES
    do
        if [ $I -gt 0 ]
        then
            OUTPUT=`diff $FILE "$FILE.dec"`
            if [ ! -z "$OUTPUT" ]
            then
                echo "\t[FAILED] Decompressed $FILE does not match original."
            fi
            rm $FILE.fac $FILE.dec
        else
            rm $FILE.sa
        fi
        I=`expr $I + 1`
    done
}

runCompress()
{
    # Test standard compression and decompression
    echo "[TEST] Standard RLZ compression and decompression"
    ./rlz $FILES
    ./rlz -d $FILES
    echo ""

    testDecompress

    # Test with LISS option
    echo "[TEST] RLZ compression with LISS option and decompression"
    ./rlz -l $FILES
    ./rlz -d $FILES
    echo ""

    testDecompress

    # Test with short factor encoding option
    echo "[TEST] RLZ compression with short factor option and decompression"
    ./rlz -s $FILES
    ./rlz -d $FILES
    echo ""

    testDecompress

    # Test with short factor encoding and LISS option
    echo "[TEST] RLZ compression with short factor and LISS options, and decompression"
    ./rlz -l -s $FILES
    ./rlz -d $FILES
    echo ""

    testDecompress
}

DISPLAYTESTDIR="$TESTDIR/display"
DISPLAYFILES="`cat $DISPLAYTESTDIR/files`"
IDXFILE="$TESTDIR/test.idx"

testDisplay()
{
    for FILE in $DISPLAYFILES
    do
        ./rlzindex d $IDXFILE < $FILE.test > $FILE.out 2> /dev/null
        OUTPUT=`diff $FILE.exp "$FILE.out"`
        if [ ! -z "$OUTPUT" ]
        then
            echo "\t[FAILED] Output for $FILE does not match expected output."
        fi
        rm $FILE.out
    done
}

runDisplay()
{
    # Test display with display only
    echo "[TEST] RLZ display with display only mode"
    ./rlz -i $IDXFILE -r $FILES > /dev/null 2>&1
    testDisplay
    rm $IDXFILE
    echo ""

    # Test display with display only
    echo "[TEST] RLZ display with whole index"
    ./rlz -i $IDXFILE $FILES > /dev/null 2>&1
    testDisplay
    rm $IDXFILE
    echo ""
}

COUNTTESTDIR="$TESTDIR/count"
COUNTFILES="`cat $COUNTTESTDIR/files`"
IDXFILE="$TESTDIR/test.idx"

testCount()
{
    for FILE in $COUNTFILES
    do
        ./rlzindex c $IDXFILE < $FILE.test > $FILE.out 2> /dev/null
        OUTPUT=`diff $FILE.exp "$FILE.out"`
        if [ ! -z "$OUTPUT" ]
        then
            echo "\t[FAILED] Output for $FILE does not match expected output."
        fi
        rm $FILE.out
    done
}

runCount()
{
    # Test count
    echo "[TEST] RLZ count"
    ./rlz -i $IDXFILE $FILES > /dev/null 2>&1
    testCount
    rm $IDXFILE
    echo ""
}

LOCATETESTDIR="$TESTDIR/locate"
LOCATEFILES="`cat $LOCATETESTDIR/files`"
IDXFILE="$TESTDIR/test.idx"

testLocate()
{
    for FILE in $LOCATEFILES
    do
        ./rlzindex l $IDXFILE < $FILE.test > $FILE.out 2> /dev/null
        # This is a bit of dodgy trick to account for the fact that the output
        # is not in the same order as the expected output
        sort $FILE.exp > $LOCATETESTDIR/e
        sort $FILE.out > $LOCATETESTDIR/o
        OUTPUT=`diff "$LOCATETESTDIR/e" "$LOCATETESTDIR/o"`
        if [ ! -z "$OUTPUT" ]
        then
            echo "\t[FAILED] Output for $FILE does not match expected output."
        fi
        rm $FILE.out $LOCATETESTDIR/e $LOCATETESTDIR/o
    done
}

runLocate()
{
    # Test locate
    echo "[TEST] RLZ locate"
    ./rlz -i $IDXFILE $FILES > /dev/null 2>&1
    testLocate
    rm $IDXFILE
    echo ""
}

runCompress
runDisplay
runCount
runLocate

# Clean up
make clobber
