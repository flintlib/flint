
#include "acb_theta.h"

void fmpz_mat_trig_sp(fmpz_mat_t m, const fmpz_mat_t s)
{  
  slong g = fmpz_mat_nrows(m)/2;
  fmpz_mat_t zero, one;
  
  fmpz_mat_init(zero, g, g);
  fmpz_mat_init(one, g, g);
  
  fmpz_mat_one(one);
  fmpz_mat_set_abcd(m, one, s, zero, one);

  fmpz_mat_clear(zero);
  fmpz_mat_clear(one);  
}
