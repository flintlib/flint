
#include "acb_theta.h"

static void
worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g,
        ulong ab, slong ord, slong prec, slong fullprec)
{
    acb_t x;
    slong sgn;
    slong k;
    ulong b, b_shift;

    acb_init(x);
  
    for (b = 0; b < n_pow(2,g); b++)
    {
        b_shift = b;
        sgn = 0;
        for (k = 0; k < g; k++)
	{
            if (b_shift & 1)
	    {
                sgn += 4 + coords[g-1-k] % 4;
	    }
            b_shift = b_shift >> 1;
	}
        sgn = sgn % 4;
      
        acb_set(x, term);
        if (sgn == 1) acb_mul_onei(x, x);
        else if (sgn == 2) acb_neg(x, x);
        else if (sgn == 3) acb_div_onei(x, x);
      
        acb_add(&th[b], &th[b], x, fullprec);
    }
  
    acb_clear(x);
}

void
acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_t epsilon;
    int all = 0;
    int unif = 0;
    slong ord = 0;
    ulong ab = 0;
    acb_ptr exp_z;
    acb_mat_t lin_powers;
    acb_t cofactor;
    slong fullprec;
    slong g = acb_mat_nrows(tau);
    slong k;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, g);
    arf_init(epsilon);
    exp_z = _acb_vec_init(g);
    acb_mat_init(lin_powers, g, g);
    acb_init(cofactor);

    acb_theta_naive_ellipsoid(E, epsilon, ab, all, unif, ord, z, tau, prec);  
    fullprec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, tau, E, fullprec);

    acb_mat_set(lin_powers, acb_theta_precomp_exp_mat(D));
    for (k = 0; k < g; k++)
    {
        acb_exp_pi_i(&exp_z[k], &z[k], prec);
    }
    acb_one(cofactor);

    for (k = 0; k < n_pow(2,g); k++) acb_zero(&th[k]);
  
    acb_theta_naive_worker_rec(th, lin_powers, E, D, exp_z, cofactor, ab, ord,
            fullprec, fullprec, worker_dim0);

    for (k = 0; k < n_pow(2,g); k++) acb_add_error_arf(&th[k], epsilon);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    arf_clear(epsilon);
    _acb_vec_clear(exp_z, g);
    acb_mat_clear(lin_powers);
    acb_clear(cofactor);  
}
