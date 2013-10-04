#include <cstdlib>
#include <iostream>
#include "arithxx.h"
#include "fmpzxx.h"

using namespace flint;
using namespace std;

int
main(int argc, char * argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: partitions n\n";
        return 1;
    }

    ulong n;
    flint_sscanf(argv[1], "%wu", &n);

    std::cout << "p(" << n << ") =\n" << number_of_partitions(n) << '\n';

    return 0;
}
