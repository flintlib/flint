#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "profiler.h"

#define cpumin 2

FLINT_DLL extern slong n_factor_pp1_table[][2];

void n_factor_pp1_table_insert(slong bits, slong B1, slong count)
{
    n_factor_pp1_table[bits][0] = B1;
    n_factor_pp1_table[bits][1] = count;
}

int
main(int argc, char** argv)
{
    double tbest = 1.0e300;
    mp_limb_t nums[1000];

    slong i;
    slong bits, B1, count;
    mp_limb_t n, cofactor;
    n_factor_t fac;

    FLINT_TEST_INIT(state);

    bits = atol(argv[1]);

    flint_printf("Looking for 1000 numbers with %ld bits cofactors\n", bits);

    for (i = 0; i < 1000; )
    {
       n_factor_init(&fac);
       n = n_randbits(state, bits + n_randint(state, FLINT_BITS - bits + 1));
       cofactor = n_factor_trial(&fac, n, FLINT_FACTOR_TRIAL_PRIMES);
       if (FLINT_BIT_COUNT(cofactor) == bits && !n_is_prime(cofactor))
       {
          nums[i++] = n;
          if (i % 100 == 0)
 	     printf("i = %ld\n", i);
       }
    }

    printf("Done computing table\n");

    for (count = 0; count < 6; count++)
    {
	for (B1 = 1; ((count == 0 && B1 == 1) || count > 0) && B1 < 10000; B1 = (slong) (1.05*(double)B1) + 1)
	{
            double t;
            int l, loops = 1;
	    n_factor_pp1_table_insert(bits, B1, count);
        loop:
            t = 0.0;
            init_clock(0);
            prof_start();
            for (l = 0; l < loops; l++)
            {
               for (i = 0; i < 1000; i++)
	       {
	           n_factor_init(&fac);
	           n_factor(&fac, nums[i], 0);
	       }
	    }
            prof_stop();
            t += get_clock(0);

            if (t * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }

	    if (t/loops < tbest)
	    {
	        flint_printf("%wd, %wd, %wd, %20f\n", bits, B1, count, t/loops);
		tbest = t/loops;
	    }
	}
    }

    FLINT_TEST_CLEANUP(state);

    return 0;
}
