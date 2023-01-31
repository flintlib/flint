
#include "acb_theta.h"

static void
acb_theta_naive_red_z(arb_ptr offset, arf_struct * eps, acb_ptr new_z,
                      acb_ptr c, acb_srcptr z, slong nb_z, const acb_mat_t tau,
                      const arb_mat_t cho, slong g, slong prec)
{
    arb_mat_t x, y, vec, r;
    arb_mat_t X, Y, Yinv;
    arb_mat_t tp, prod;
    arb_t bound;
    slong *v;
    slong k, j;

    arb_mat_init(vec, g, 1);
    arb_mat_init(x, g, 1);
    arb_mat_init(y, g, 1);
    arb_mat_init(r, g, 1);
    arb_mat_init(X, g, g);
    arb_mat_init(Y, g, g);
    arb_mat_init(Yinv, g, g);
    arb_mat_init(tp, 1, g);
    arb_mat_init(prod, 1, 1);
    arb_init(bound);
    v = flint_malloc(g * sizeof(slong));

    acb_mat_get_real(X, tau);
    acb_mat_get_imag(Y, tau);
    arb_mat_inv(Yinv, Y, prec);

    for (k = 0; k < nb_z; k++)
    {
        for (j = 0; j < g; j++)
        {
            arb_set(arb_mat_entry(x, j, 0), acb_realref(&z[k * g + j]));
            arb_set(arb_mat_entry(y, j, 0), acb_imagref(&z[k * g + j]));
        }
        /* Get center of ellipsoid, update multiplicative factor */
        arb_mat_mul(vec, Yinv, y, prec);
        acb_zero(&c[k]);
        arb_mat_transpose(tp, y);
        arb_mat_mul(prod, tp, vec, prec);
        arb_sub(acb_imagref(&c[k]), acb_imagref(&c[k]),
                arb_mat_entry(prod, 0, 0), prec);

        /* Get multiplier on error bound */
        arb_const_pi(bound, prec);
        arb_mul(bound, bound, arb_mat_entry(prod, 0, 0), prec);
        arb_exp(bound, bound, prec);
        arb_get_ubound_arf(&eps[k], bound, prec);

        /* Round to nearest even integer vector v */
        arb_mat_scalar_mul_2exp_si(vec, vec, -1);
        acb_theta_eld_round(v, vec);
        for (j = 0; j < g; j++)
            v[j] *= 2;
        arb_mat_scalar_mul_2exp_si(vec, vec, 1);

        /* Get r and uniform offset */
        for (j = 0; j < g; j++)
        {
            arb_sub_si(arb_mat_entry(r, j, 0), arb_mat_entry(vec, j, 0),
                       v[j], prec);
        }
        arb_mat_mul(vec, cho, r, prec);
        for (j = 0; j < g; j++)
        {
            if (k == 0)
            {
                arb_set(&offset[j], arb_mat_entry(vec, j, 0));
            }
            else
            {
                arb_union(&offset[j], &offset[j],
                          arb_mat_entry(vec, j, 0), prec);
            }
        }

        /* Complete multiplicative factor and new_z */
        for (j = 0; j < g; j++)
        {
            acb_set_arb(&new_z[k * g + j], arb_mat_entry(x, j, 0));
            arb_set_si(arb_mat_entry(tp, 0, j), v[j]);
        }
        arb_mat_mul(prod, tp, x, prec);
        arb_mul_2exp_si(arb_mat_entry(prod, 0, 0),
                        arb_mat_entry(prod, 0, 0), 1);
        acb_sub_arb(&c[k], &c[k], arb_mat_entry(prod, 0, 0), prec);

        for (j = 0; j < g; j++)
        {
            arb_set_si(arb_mat_entry(vec, j, 0), v[j]);
        }
        arb_mat_transpose(tp, vec);
        arb_mat_mul(vec, X, vec, prec);
        for (j = 0; j < g; j++)
        {
            acb_sub_arb(&new_z[k * g + j], &new_z[k * g + j],
                        arb_mat_entry(vec, j, 0), prec);
        }
        arb_mat_mul(prod, tp, vec, prec);
        acb_add_arb(&c[k], &c[k], arb_mat_entry(prod, 0, 0), prec);

        arb_mat_mul(vec, Y, r, prec);
        for (j = 0; j < g; j++)
        {
            arb_add(acb_imagref(&new_z[k * g + j]),
                    acb_imagref(&new_z[k * g + j]), arb_mat_entry(vec, j, 0),
                    prec);
        }
        arb_mat_transpose(tp, r);
        arb_mat_mul(prod, tp, vec, prec);
        arb_add(acb_imagref(&c[k]), acb_imagref(&c[k]),
                arb_mat_entry(prod, 0, 0), prec);
        acb_exp_pi_i(&c[k], &c[k], prec);
    }

    arb_mat_clear(vec);
    arb_mat_clear(x);
    arb_mat_clear(y);
    arb_mat_clear(r);
    arb_mat_clear(X);
    arb_mat_clear(Y);
    arb_mat_clear(Yinv);
    arb_mat_clear(tp);
    arb_mat_clear(prod);
    arb_clear(bound);
    flint_free(v);
}

void
acb_theta_naive_ellipsoid(acb_theta_eld_t E, arf_struct * eps, acb_ptr c,
                          acb_ptr new_z, ulong ab, int all, slong ord,
                          acb_srcptr z, slong nb_z, const acb_mat_t tau,
                          slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong eld_prec = ACB_THETA_ELD_DEFAULT_PREC;
    arb_t pi;
    arf_t R2, bound;
    slong scl = -1;
    arb_mat_t cho;
    arb_ptr offset;
    int res;
    slong k;

    arb_init(pi);
    arf_init(R2);
    arf_init(bound);
    arb_mat_init(cho, g, g);
    offset = _arb_vec_init(g);

    arf_one(bound);
    arf_mul_2exp_si(bound, bound, -prec + ACB_THETA_NAIVE_EPS_2EXP);

    if (all)
    {
        ab = 0;
        scl = -2;
    }

    acb_mat_get_imag(cho, tau);
    arb_const_pi(pi, prec);
    arb_mat_scalar_mul_arb(cho, cho, pi, prec);

    res = arb_mat_cho(cho, cho, eld_prec);
    if (!res)
    {
        eld_prec = prec;
        res = arb_mat_cho(cho, cho, eld_prec);
    }
    if (!res)
    {
        flint_printf("acb_theta_naive_ellipsoid: Error ");
        flint_printf("(imaginary part is not positive definite)\n");
        fflush(stdout);
        flint_abort();
    }

    arb_mat_transpose(cho, cho);
    acb_theta_naive_radius(R2, cho, ord, bound, eld_prec);

    /* Set offset in terms of z */
    acb_theta_naive_red_z(offset, eps, new_z, c, z, nb_z, tau, cho, g, prec);
    for (k = 0; k < nb_z; k++)
    {
        arf_mul(&eps[k], &eps[k], bound, prec, ARF_RND_CEIL);
    }

    /* Fill ellipsoid */
    arb_mat_scalar_mul_2exp_si(cho, cho, scl);
    acb_theta_eld_fill(E, cho, R2, offset, NULL, ab >> g, eld_prec);

    arb_clear(pi);
    arf_clear(R2);
    arf_clear(bound);
    arb_mat_clear(cho);
    _arb_vec_clear(offset, g);
}
