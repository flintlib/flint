
#include "acb_theta.h"

void acb_theta_duplication(acb_ptr th2, acb_srcptr th, slong g, slong prec)
{
  /* To be replaced by a more efficient algorithm using dyadic convolutions */
  acb_mat_t row, col, prod;
  slong nb = n_pow(2,g);
  ulong b, bp, bpp;
  
  acb_mat_init(row, 1, nb);
  acb_mat_init(col, nb, 1);
  acb_mat_init(prod, nb, nb);

  for (b = 0; b < nb; b++)
    {
      acb_set(acb_mat_entry(row, 0, b), &th[b]);
      acb_set(acb_mat_entry(col, b, 0), &th[b]);
    }
  acb_mat_mul(prod, col, row, prec);

  for (b = 0; b < nb; b++)
    {
      for (bp = 0; bp < nb; bp++)
	{
	  bpp = b ^ bp; /* bitwise xor */
	  acb_add(&th2[b], &th2[b], acb_mat_entry(prod, bp, bpp), prec);
	}
      acb_mul_2exp_si(&th2[b], &th2[b], -g);
    }

  acb_mat_clear(row);
  acb_mat_clear(col);
  acb_mat_clear(prod);  
}
