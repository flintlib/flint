
#include "acb_theta.h"

/* Assuming a >= 0, return R such that x - (a/2)*log(x)\geq b for all
   x\geq R, and R is close to the smallest possible */

static void invert_lin_plus_log(arf_t R, slong a, const arb_t b, slong prec)
{
  arb_t x, y;
  arf_t z;
  slong k;
    
  arb_init(x);
  arb_init(y);
  arf_init(z);
  
  if (a == 0)
    {
      arb_get_ubound_arf(R, b, prec);
      goto exit;
    }
  if (!arb_is_finite(b))
    {
      arf_nan(R);
      goto exit;
    }
  
  /* Now a>0 and b finite; minimum is at x=a/2 */  
  arb_set_si(x, a);
  arb_div_si(x, x, 2, prec);
  arb_log(y, x, prec);
  arb_mul(y, y, x, prec);
  arb_sub(y, x, y, prec);
  
  if (arb_lt(b, y))
    {
      arf_zero(R);
      goto exit;
    }
  
  /* Otherwise, x = max(a, 2*(b - min)) is always large enough; then
     iterate function a few times */
  arb_sub(y, b, y, prec);
  arb_max(y, y, x, prec);
  arb_mul_si(x, y, 2, prec);
  arb_get_ubound_arf(z, x, prec);
  arb_set_arf(x, z);

  for (k = 0; k < 4; k++)
    {
      arb_log(y, x, prec);
      arb_mul_si(y, y, a, prec);
      arb_div_si(y, y, 2, prec);
      arb_add(x, b, y, prec);
      arb_get_ubound_arf(z, x, prec);
      arb_set_arf(x, z);
    }

  arb_get_ubound_arf(R, x, prec);
  goto exit;
  
 exit:
  {
    arb_clear(x);
    arb_clear(y);
    arf_clear(z);
  }
}

void acb_theta_naive_radius(arf_t R, const arb_mat_t Y, slong p, const arf_t epsilon, slong prec)
{
  arb_t b, temp;
  arf_t cmp;
  slong g = arb_mat_nrows(Y);
  slong k;

  arb_init(b);
  arb_init(temp);
  arf_init(cmp);

  /* Divide epsilon by right factors to reduce to invert_lin_plus_log */
  arb_set_arf(b, epsilon);
  arb_mul_2exp_si(b, b, -2*g-2);

  for (k = 0; k < g; k++)
    {
      arb_inv(temp, arb_mat_entry(Y, k, k), prec);
      arb_add_si(temp, temp, 1, prec);
      arb_div(b, b, temp, prec);
    }

  /* Solve R^((g-1)/2+p) exp(-R) \leq b */
  arb_log(b, b, prec);
  arb_neg(b, b);
  invert_lin_plus_log(R, g-1+2*p, b, prec);
  
  /* Max with 4, 2*p for formula to be valid */
  arf_set_si(cmp, FLINT_MAX(4, 2*p));
  arf_max(R, R, cmp);

  arb_clear(b);
  arb_clear(temp);
  arf_clear(cmp);
}
