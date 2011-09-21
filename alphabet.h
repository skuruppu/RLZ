const char *NUCLALPHA = "acgnt";
const uint64_t BITSPERBASE = 3;
const uint64_t NUCLALPHASIZE = 5;

unsigned char nucl_to_int[256] = {0};
unsigned char int_to_nucl[NUCLALPHASIZE+1];
unsigned char int_to_2bpb[NUCLALPHASIZE+1];
unsigned char bpb_to_char[NUCLALPHASIZE+1];

void initialise_nucl_converters()
{
    uint64_t i;

    for (i=0; i<NUCLALPHASIZE; i++)
    {
        nucl_to_int[tolower(NUCLALPHA[i])] = i+1;
        nucl_to_int[toupper(NUCLALPHA[i])] = i+1;
    }

    for (i=0; i<NUCLALPHASIZE; i++)
    {
        int_to_nucl[i+1] = NUCLALPHA[i];
    }

    int_to_2bpb[1] = 0; // a
    int_to_2bpb[2] = 1; // c
    int_to_2bpb[3] = 2; // g
    int_to_2bpb[4] = -1; // n (invalid)
    int_to_2bpb[5] = 3; // t

    bpb_to_char[0] = 'a';
    bpb_to_char[1] = 'c';
    bpb_to_char[2] = 'g';
    bpb_to_char[3] = 't';
}
