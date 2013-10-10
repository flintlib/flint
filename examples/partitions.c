#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"

int
main(int argc, char * argv[])
{
    fmpz_t x;
    ulong n;

    if (argc != 2)
    {
        flint_printf("usage: partitions n\n");
        return 1;
    }

    flint_sscanf(argv[1], "%wu", &n);

    flint_printf("p(%wu) = \n", n);

    fmpz_init(x);
    arith_number_of_partitions(x, n);
    fmpz_print(x);
    flint_printf("\n");
    fmpz_clear(x);

    return 0;
}
