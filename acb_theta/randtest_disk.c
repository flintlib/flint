
#include "acb_theta.h"

void acb_randtest_disk(acb_t x, const acb_t ctr, const arf_t rad, flint_rand_t state, slong prec)
{
  arb_t half;
  arb_t err;
  acb_t diff;

  arb_init(half);
  arb_init(err);
  acb_init(diff);

  arb_one(half);
  arb_mul_2exp_si(half, half, -1);

  acb_get_mid(x, ctr);

  arb_urandom(err, state, prec);
  arb_sub(err, err, half, prec);
  arb_mul_arf(err, err, rad, prec);
  arb_add(acb_realref(x), acb_realref(x), err, prec);

  arb_urandom(err, state, prec);
  arb_sub(err, err, half, prec);
  arb_mul_arf(err, err, rad, prec);
  arb_add(acb_imagref(x), acb_imagref(x), err, prec);

  /* Check that x is indeed included in disk */
  acb_sub(diff, ctr, x, prec);
  acb_abs(err, diff, prec);
  arb_set_arf(half, rad);
  if (!arb_lt(err, half))
    {
      /* Simply set x to ctr; abort if radius is larger than rad */
      acb_get_mid(x, ctr);
      acb_sub(diff, x, ctr, prec);
      acb_abs(err, diff, prec);
      arb_set_arf(half, rad);
      if (!arb_lt(err, half))
	{
	  flint_printf("acb_randtest_disk: Error (center is too imprecise)\n");
	  fflush(stdout);
	  flint_abort();
	}
      acb_set(x, ctr);
    }
  
  arb_clear(half);
  arb_clear(err);
  acb_clear(diff);
}
