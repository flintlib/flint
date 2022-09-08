
#include "acb_theta.h"

/* Work in dimension 1: compute exponentiel terms with two
   multiplications per term only, at just the necessary precision.
   Each term is: cofactor * lin^k * x^(k^2), and square
   powers of x are precomputed.
*/
void acb_theta_naive_worker_dim1(acb_ptr th,
				 const acb_theta_eld_t E, const acb_theta_precomp_t D,
				 const acb_t lin, const acb_t cofactor,
				 ulong ab, slong ord, slong prec, slong fullprec,
				 acb_theta_naive_worker_t worker_dim0)
{
  acb_t start, diff, aff, term;
  slong* coords;
  slong g = acb_theta_eld_ambient_dim(E);
  slong min = acb_theta_eld_min(E);
  slong mid = acb_theta_eld_mid(E);
  slong max = acb_theta_eld_max(E);
  slong newprec;
  slong k;

  if (acb_theta_eld_nb_pts(E) == 0) {return;}

  acb_init(start);
  acb_init(diff);
  acb_init(aff);
  acb_init(term);
  coords = flint_malloc(g * sizeof(slong));

  for (k = 1; k < g; k++)
    {
      coords[k] = acb_theta_eld_coord(E,k);
    }

  acb_pow_si(start, lin, mid, prec);
  acb_mul(start, start, cofactor, prec);
  acb_pow_si(diff, lin, 2, prec);
  
  acb_set(aff, start);
  for (k = mid; k <= max; k += 2)
    {
      coords[0] = k;
      newprec = acb_theta_naive_newprec(prec, k, k-mid, max-mid, ord);
      if (k > mid) acb_mul(aff, aff, diff, newprec);
      
      acb_mul(term, aff, acb_theta_precomp_sqr_pow(D, 0, FLINT_ABS(k)/2), newprec);
      worker_dim0(th, term, coords, g, ab, ord, newprec, fullprec);
    }

  acb_set(aff, start);
  acb_inv(diff, diff, prec);
  for (k = mid - 2; k >= min; k -= 2)
    {
      coords[0] = k;
      newprec = acb_theta_naive_newprec(prec, k, mid-k, mid-min, ord);      
      acb_mul(aff, aff, diff, newprec);
      
      acb_mul(term, aff, acb_theta_precomp_sqr_pow(D, 0, FLINT_ABS(k)/2), newprec);
      worker_dim0(th, term, coords, g, ab, ord, newprec, fullprec);      
    }

  acb_clear(start);
  acb_clear(diff);
  acb_clear(aff);
  acb_clear(term);
  flint_free(coords);
}
