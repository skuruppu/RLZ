#include <getopt.h>
#include "RLZ.h"

using namespace std;

int main(int argc, char **argv)
{
	char usage[] = "Usage: %s [OPTIONS] REF FILE1 FILE2 ...\n\
	-e: Type of encoding (t: text, b: binary, l: liss) (default: b)\n\
	REF: Name of reference sequence\n\
	FILE1 ...: Names of files to be compressed\n";
	char option, encoding='b';

	if (argc < 3)
	{
		cerr << usage;
		exit(1);
	}

    while ((option = getopt(argc, argv, "e:")) != EOF)
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
            default:
                cerr << usage;
                exit(1);
        }
    }

	RLZ rlz(argv+optind, argc-optind, encoding);

    rlz.compress();

    return 0;
}
