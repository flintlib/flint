
#include "acb_theta.h"

void acb_siegel_randtest_fund(acb_mat_t tau, flint_rand_t state, slong prec)
{
  slong g = arb_mat_nrows(tau);
  arf_t rad;
  acb_t err;
  acb_t c;
  slong k, j;
  
  arf_init(rad);
  acb_init(err);
  acb_init(c);

  arf_one(rad);
  arf_div_si(rad, rad, 2*g, prec, ARF_RND_FLOOR);

  acb_mat_zero(tau);
  for (k = 0; k < g; k++)
    {
      acb_onei(c);
      acb_mul_si(c, c, k+3, prec);
      acb_mul_2exp_si(c, c, -1);
      acb_set(acb_mat_entry(tau, k, k), c);
    }

  acb_zero(c);
  for (k = 0; k < g; k++)
    {
      for (j = 0; j < g; j++)
	{
	  acb_randtest_disk(err, c, rad, state, prec);
	  acb_add(acb_mat_entry(tau, k, j), acb_mat_entry(tau, k, j), err, prec);
	}
    }
  
  arf_clear(rad);
  acb_clear(err);
  acb_clear(c);
}
