/* 
   Check that n_is_prime returns 'composite' for all base-2 Fermat pseudoprimes
   341, 561, 645, 1105, 1387, 1729, 1905, 2047, ... (OEIS A001567)
   less than 2^FLINT_BITS. In effect, this verifies correctness of the BPSW
   test in n_is_prime assuming that n_is_prime always does at least a
   base-2 Fermat test (or a base-2 Miller-Rabin test, which is strictly
   stronger).

   The user must provide the path to a text file containing all the
   pseudoprimes, obtainable from http://www.cecm.sfu.ca/Pseudoprimes/
   (885 MB bz2-compressed; 2.35 GB uncompressed).

   This will typically take about a minute to run.

   This file is public domain. Author: Fredrik Johansson.
*/

#include <stdio.h>
#include "flint/flint.h"
#include "flint/ulong_extras.h"
#include "flint/profiler.h"

int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        flint_printf("usage: check_n_is_prime <psps-below-2-to-64.txt>\n");
        return 1;
    }

    FILE * fp = fopen(argv[1], "r");
    ulong n, sum = 0;
    slong count = 0;

    if (fp == NULL)
    {
        flint_printf("unable to open file\n");
        return 1;
    }
    TIMEIT_ONCE_START

    while (!feof(fp) && fscanf(fp, "%lu", &n) == 1)
    {
        count++;

        if (count % 1000000 == 0)
            flint_printf("%wd: %wu  (%d bits)\n", count, n, FLINT_BIT_COUNT(n));

        sum += n;

        if (n_is_prime(n))
        {
            flint_printf("FAIL: %wu claimed prime\n", n);
            flint_abort();
        }

        if (FLINT_BITS == 32 && n == UWORD(4294901761))
            break;
    }

    TIMEIT_ONCE_STOP

    slong expect_count = (FLINT_BITS == 64) ? 118968378 : 10403;
    ulong expect_sum = (FLINT_BITS == 64) ? UWORD(6235045495123121954) : UWORD(882877973);

    if (count != expect_count || sum != expect_sum)
    {
        flint_printf("FAIL: read %wd pseudoprimes with checksum %wu, expected %wd with checksum %wu\n",
            count, sum, expect_count, expect_sum);
        flint_abort();
    }

    fclose(fp);        
    return 0;
}

