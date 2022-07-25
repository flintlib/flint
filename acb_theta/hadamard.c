
#include "acb_theta.h"

void acb_theta_hadamard(fmpz_mat_t h)
{
  slong k, n;
  fmpz_mat_t src, dest;
  slong g = n_clog(fmpz_mat_nrows(h), 2);

  fmpz_mat_zero(h);
  fmpz_one(fmpz_mat_entry(h, 0, 0));
  fmpz_one(fmpz_mat_entry(h, 0, 1));
  fmpz_one(fmpz_mat_entry(h, 1, 0));
  fmpz_one(fmpz_mat_entry(h, 1, 1));
  fmpz_neg(fmpz_mat_entry(h, 1, 1));

  for (k = 1; k < g; k++)
    {
      n = n_pow(2,k);
      fmpz_mat_window_init(src, h, 0, 0, n, n);
      
      fmpz_mat_window_init(dest, h, 0, n, n, 2*n);
      fmpz_mat_set(dest, src);
      fmpz_mat_window_clear(dest);

      fmpz_mat_window_init(dest, h, n, 0, 2*n, n);
      fmpz_mat_set(dest, src);
      fmpz_mat_window_clear(dest);
      
      fmpz_mat_window_init(dest, h, n, n, 2*n, 2*n);
      fmpz_mat_neg(dest, src);
      fmpz_mat_window_clear(dest);
    }
}
