
#include "acb_theta.h"

void acb_mat_ninf(arb_t norm, const acb_mat_t m, slong prec)
{
  slong n = acb_mat_nrows(m);
  slong k = acb_mat_ncols(m);
  arb_t abs, sum;
  slong i, j;

  arb_init(abs);
  arb_init(sum);

  arb_zero(norm);
  for (i = 0; i < n; i++)
    {
      arb_zero(sum);
      for (j = 0; j < k; j++)
	{
	  acb_abs(abs, acb_mat_entry(m,i,j), prec);
	  arb_add(sum, sum, abs, prec);
	}
      arb_max(norm, norm, sum, prec);
    }
  
  arb_clear(abs);
  arb_clear(sum);
}
