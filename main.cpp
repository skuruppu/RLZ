#include "RLZ.h"

using namespace std;

int main(int argc, char **argv)
{
	ifstream infile;
	infile.open(argv[1], ifstream::in);

	RLZ rlz(infile);

    return 0;
}
