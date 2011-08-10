#include <getopt.h>
#include "RLZ.h"

using namespace std;

int main(int argc, char **argv)
{
	char usage[] = "Usage: %s [OPTIONS] REF FILE1 FILE2 ...\n\
    -d: Decompress (all other options ignored)\n\
	-e: Type of encoding (t: text, b: binary) (default: b)\n\
    -l: Enable LISS encoding\n\
    -s: Output short factors as substring and length pairs\n\
	REF: Name of reference sequence\n\
	FILE1 ...: Names of files to be compressed\n";
	char option, encoding='b';
    bool isdecomp = false, isshort = false, isliss = false;

	if (argc < 3)
	{
		cerr << usage;
		exit(1);
	}

    while ((option = getopt(argc, argv, "e:dls")) != EOF)
    {
        switch (option)
        {
            case 'e':
                sscanf(optarg, "%c", &encoding);
				if (encoding != 't' && encoding != 'b' &&
					encoding != 'l')
				{
					cerr << usage;
					exit(1);
				}
                break;
            case 's':
                isshort = true;
                break;
            case 'd':
                isdecomp = true;
                break;
            case 'l':
                isliss = true;
                break;
            default:
                cerr << usage;
                exit(1);
        }
    }

    if (!isdecomp)
    {
        RLZCompress rlz(argv+optind, argc-optind, encoding, isshort,
                        isliss);

        rlz.compress();
    }
    else
    {
        RLZDecompress rlzd(argv+optind, argc-optind);

        rlzd.decompress();
    }

    return 0;
}
