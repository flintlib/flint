
#include "acb_theta.h"

void
arb_randtest_pos(arb_t x, flint_rand_t state, slong prec, slong mag_bits)
{
    int pos = 0;
    while (!pos)
    {
	arb_randtest_precise(x, state, prec, mag_bits);
	pos = arb_is_positive(x);
    }
}
