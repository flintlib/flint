
#include "acb_theta.h"

void acb_theta_transform_sqr_proj(acb_ptr r, acb_srcptr th, const fmpz_mat_t N, slong prec)
{
  acb_ptr res;
  slong g = fmpz_mat_nrows(N)/2;
  ulong n = 1<<g;
  ulong ab;
  ulong image_ab;
  fmpz_t epsilon;
  acb_t c;

  res = _acb_vec_init(n);
  fmpz_init(epsilon);
  acb_init(c);

  for (ab = 0; ab < n; ab++)
    {
      image_ab = acb_theta_transform_image_char(epsilon, ab, N);
      acb_unit_root(c, 4, prec);
      acb_pow_fmpz(c, c, epsilon, prec);
      acb_mul(c, c, &th[image_ab], prec);
      acb_set(&res[ab], c);
    }

  _acb_vec_set(r, res, n);

  _acb_vec_clear(res, n);
  fmpz_clear(epsilon);
  acb_clear(c);
}
