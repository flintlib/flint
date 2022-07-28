
#include "acb_theta.h"

void fmpz_mat_get_c(fmpz_mat_t c, const fmpz_mat_t m)
{
  slong g = fmpz_mat_nrows(m)/2;
  slong j, k;
  for (j = 0; j < g; j++)
    {
      for (k = 0; k < g; k++)
	{
	  fmpz_set(fmpz_mat_entry(c, j, k),
		   fmpz_mat_entry(m, j+g, k));
	}
    }
}
