
#include "acb_theta.h"

void acb_siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec)
{
  slong g = fmpz_mat_nrows(m)/2;
  fmpz_mat_t a;
  acb_mat_t x, num, den, invden;
  arb_mat_t im;
  int res;
  slong j, k;

  fmpz_mat_init(a, g, g);  
  acb_mat_init(x, g, g);
  acb_mat_init(num, g, g);
  acb_mat_init(den, g, g);
  acb_mat_init(invden, g, g);
  arb_mat_init(im, g, g);

  fmpz_mat_get_a(a, m);
  acb_mat_set_fmpz_mat(x, a);
  acb_mat_mul(num, x, z, prec);
  fmpz_mat_get_b(a, m);
  acb_mat_set_fmpz_mat(x, a);
  acb_mat_add(num, num, x, prec);

  acb_siegel_cocycle(den, m, z, prec);
  res = acb_mat_inv(invden, den, prec);
  if (!res)
    {
      for (j = 0; j < g; j++)
	{
	  for (k = 0; k < g; j++) acb_indeterminate(acb_mat_entry(invden,j,k));
	}
    }

  acb_mat_mul(w, num, invden, prec);
  acb_mat_get_imag(im, w);    

  fmpz_mat_clear(a);
  acb_mat_clear(x);
  acb_mat_clear(num);
  acb_mat_clear(den);
  acb_mat_clear(invden);
  arb_mat_clear(im);
}
