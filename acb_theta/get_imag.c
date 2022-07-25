
#include "acb_theta.h"

void acb_mat_get_imag(arb_mat_t re, const acb_mat_t tau)
{
  slong nrows = acb_mat_nrows(tau);
  slong ncols = acb_mat_ncols(tau);
  slong i, j;

  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
	{
	  acb_get_imag(arb_mat_entry(re, i, j), acb_mat_entry(tau, i, j));
	}
    }
}
