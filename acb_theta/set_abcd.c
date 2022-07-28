
#include "acb_theta.h"

void fmpz_mat_set_abcd(fmpz_mat_t m,
		       const fmpz_mat_t a, const fmpz_mat_t b,
		       const fmpz_mat_t c, const fmpz_mat_t d)
{
  slong g = fmpz_mat_nrows(m)/2;
  slong j, k;
  for (j = 0; j < g; j++)
    {
      for (k = 0; k < g; k++)
	{
	  fmpz_set(fmpz_mat_entry(m, j, k),
		   fmpz_mat_entry(a, j, k));
	  fmpz_set(fmpz_mat_entry(m, j, k+g),
		   fmpz_mat_entry(b, j, k));
	  fmpz_set(fmpz_mat_entry(m, j+g, k),
		   fmpz_mat_entry(c, j, k));
	  fmpz_set(fmpz_mat_entry(m, j+g, k+g),
		   fmpz_mat_entry(d, j, k));
	}
    }
}
