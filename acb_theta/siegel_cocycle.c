
#include "acb_theta.h"

void acb_siegel_cocycle(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec)
{
  slong g = fmpz_mat_nrows(m)/2;
  fmpz_mat_t cd;
  acb_mat_t r, s;

  fmpz_mat_init(cd, g, g);
  acb_mat_init(r, g, g);
  acb_mat_init(s, g, g);

  fmpz_mat_get_c(cd, m);
  acb_mat_set_fmpz_mat(r, cd);
  acb_mat_mul(r, r, z, prec);
  fmpz_mat_get_d(cd, m);
  acb_mat_set_fmpz_mat(s, cd);
  acb_mat_add(r, r, s, prec);

  acb_mat_set(w, r);

  fmpz_mat_clear(cd);
  acb_mat_clear(r);
  acb_mat_clear(s);
}
