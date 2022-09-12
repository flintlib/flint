
#include "acb_theta.h"

static void
worker_dim0(acb_ptr th, const acb_t term, slong* coords, slong g,
        ulong ab, slong ord, slong prec, slong fullprec)
{
    acb_t x;
    slong sgn;
    ulong b;

    acb_init(x);
  
    for (b = 0; b < n_pow(2,g); b++)
    {
        sgn = acb_theta_dot(b, coords, g) % 4;
      
        acb_set(x, term);
        if (sgn == 1) acb_mul_onei(x, x);
        else if (sgn == 2) acb_neg(x, x);
        else if (sgn == 3) acb_div_onei(x, x);
      
        acb_add(&th[b], &th[b], x, fullprec);
    }

    flint_printf("Coordinates");
    for (b = 0; b < g; b++) flint_printf(" %wd", coords[b]);
    flint_printf("\n");
    acb_printd(term, 10);
    flint_printf("\n");
  
    acb_clear(x);
}

void
acb_theta_naive(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_theta_eld_t E;
    acb_theta_precomp_t D;
    arf_t epsilon;
    int all = 0;
    int unif = 0;
    slong ord = 0;
    slong k = 0;
    ulong ab = 0;
    slong nb = 1<<g;

    acb_theta_eld_init(E, g, g);
    acb_theta_precomp_init(D, 1, g);
    arf_init(epsilon);

    acb_theta_naive_ellipsoid(E, epsilon, ab, all, unif, ord, z, tau, prec);  
    prec = acb_theta_naive_fullprec(E, prec);
    acb_theta_precomp_set(D, z, tau, E, prec);
  
    acb_theta_naive_worker(th, nb, epsilon, E, D, k, ab, ord, prec,
            worker_dim0);

    acb_theta_eld_clear(E);
    acb_theta_precomp_clear(D);
    arf_clear(epsilon);
}
