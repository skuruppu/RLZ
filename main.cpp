#include "RLZ.h"

using namespace std;

int main(int argc, char **argv)
{
	RLZ rlz(argv+1, argc-1);

    rlz.compress();

    return 0;
}
