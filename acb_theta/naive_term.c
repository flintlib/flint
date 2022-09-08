
#include "acb_theta.h"

void acb_theta_naive_term(acb_t exp, const acb_mat_t tau, acb_srcptr z,
        ulong ab, slong* coords, slong prec)
{
    slong j, k;
    slong g = acb_mat_nrows(tau);
    slong* a;
    slong* b;
    acb_t x, t;
    arb_t pi;

    acb_init(x);
    acb_init(t);
    a = flint_malloc(g * sizeof(slong));
    b = flint_malloc(g * sizeof(slong));
    arb_init(pi);

    for (k = 0; k < g; k++)
    {
        b[g-1-k] = ab % 2;
        ab = ab >> 1;
    }
    for (k = 0; k < g; k++)
    {
        a[g-1-k] = ab % 2;
        ab = ab >> 1;
    }

    acb_zero(x);
    for (k = 0; k < g; k++)
    {
        acb_mul_si(t, acb_mat_entry(tau,k,k), n_pow(2*coords[k] + a[k], 2), prec);
        acb_add(x, x, t, prec);
      
        acb_mul_2exp_si(t, &z[k], 1);
        acb_add_si(t, t, b[k], prec);
        acb_mul_si(t, t, 2*coords[k] + a[k], prec);
        acb_add(x, x, t, prec);

        for (j = k+1; j < g; j++)
	{
            acb_mul_si(t, acb_mat_entry(tau,k,j), 2*(2*coords[k] + a[k])*(2*coords[j] + a[j]),
                    prec);
            acb_add(x, x, t, prec);
	}
    }
    acb_mul_2exp_si(x, x, -2);

    acb_mul_onei(x, x);
    arb_const_pi(pi, prec);
    acb_mul_arb(x, x, pi, prec);
    acb_exp(exp, x, prec);

    acb_clear(x);
    acb_clear(t);
    flint_free(a);
    flint_free(b);
    arb_clear(pi);
}
  

  
    

  
  
