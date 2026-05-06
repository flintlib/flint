#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_poly.h"
#include "acb_ode.h"
#include "fmpz.h"
#include "fmpz_vec.h"

static void
next_binom_n(fmpz * binom_n, slong n, slong len)
{
    if (len == 0)
        return;

    if (n == 0)
    {
        _fmpz_vec_zero(binom_n, len);
        fmpz_one(binom_n);
    }
    else
    {
        for (slong i = len - 1; i > 0; i--)
            fmpz_add(binom_n + i, binom_n + i, binom_n + i - 1);
    }
}

static void
next_coeff(acb_ode_sol_t sol, slong base, slong n,
           acb_poly_struct * rhs, slong rhs_nlogs,
           slong mult, slong ini_i, const acb_poly_t ind_n,
           int approx, slong prec)
{
    slong k;
    TMP_INIT;
    TMP_START;

    while (rhs_nlogs > 0
           && acb_is_zero((rhs + rhs_nlogs - 1)->coeffs + n - base))
        rhs_nlogs--;

    /* Compute the high-log part of the new term from the coefficient of x^{λ+n}
     * on the rhs and the indicial polynomial */

    acb_ptr new_term = TMP_ALLOC(sizeof(acb_t) * rhs_nlogs);
    acb_ptr rhs_term = TMP_ALLOC(sizeof(acb_t) * rhs_nlogs);

    for (slong k = 0; k < rhs_nlogs; k++) {
        acb_init(new_term + k);
        rhs_term[k] = (rhs + k)->coeffs[n - base];
    }

    _acb_ode_poly_negdivrevhigh(new_term, ind_n->coeffs + mult, NULL, rhs_term,
                                rhs_nlogs, prec);

    /* todo: also do this in rigorous mode but keep the radii for linear error
       analysis */
    if (approx)
    {
        for (slong k = 0; k < rhs_nlogs; k++)
        {
            mag_zero(arb_radref(acb_realref(new_term + k)));
            mag_zero(arb_radref(acb_imagref(new_term + k)));
        }
    }

    /* Write the new term to sol, store new extini */

    for (k = 0; k < mult; k++)  /* initial conditions */
    {
        acb_set((sol->series + k)->coeffs + n - base,
                acb_mat_entry(sol->extini, ini_i, k));
    }
    for (; k < rhs_nlogs + mult; k++)
    {
        acb_swap((sol->series + k)->coeffs + n - base, new_term + k - mult);

        if (mult > 0)
            acb_set(acb_mat_entry(sol->extini, ini_i, k),
                    (sol->series + k)->coeffs + n - base);
    }
    for (; k < sol->nlogs; k++)
    {
        acb_zero((sol->series + k)->coeffs + n - base);
    }

    /* Update log-degree and used initial values */

    sol->nlogs = FLINT_MAX(sol->nlogs, rhs_nlogs + mult);

    for (slong k = 0; k < rhs_nlogs; k++)
        acb_clear(new_term + k);
    TMP_END;
}

static void
next_sums(acb_ode_sol_t sol, slong base, slong n,
          acb_srcptr pows, const char * shifted, slong npts,
          const fmpz * binom_n,
          slong nder, slong prec)
{
    acb_t tmp;
    acb_init(tmp);

    for (slong j = 0; j < npts; j++)
    {
        for (slong k = 0; k < sol->nlogs; k++)
        {
            acb_ptr c = (sol->series + k)->coeffs + n - base;
            acb_ptr s = acb_ode_sol_sum_ptr(sol, j, k, 0);

            if (shifted[j])
            {
                slong pows_offset = (n % nder) * npts;
                acb_srcptr ptpow = pows + pows_offset + j;
                acb_mul(tmp, c, ptpow, prec);
                for (slong i = 0; i < nder; i++)
                    acb_addmul_fmpz(s + i, tmp, binom_n + i, prec);
            }
            else
            {
                for (slong i = 0; i < nder; i++)
                {
                    slong pows_offset = ((n + nder - i) % nder) * npts;
                    acb_srcptr ptpow = pows + pows_offset + j;
                    acb_mul(tmp, c, ptpow, prec);
                    acb_addmul_fmpz(s + i, tmp, binom_n + i, prec);
                }
            }
        }
    }

    acb_clear(tmp);
}

void
_acb_ode_sum_forward_1(acb_ode_sum_struct * sum)
{
    slong ini_i, mult, len;
    acb_poly_t ind_n;

    acb_poly_init(ind_n);

    slong n = sum->n;

    /* find the position and multiplicity of the current shift as an indicial
       root */

    ini_i = -1;
    mult = 0;
    for (slong i = 0; i < sum->dop_len - 1; i++)
    {
        /* flint_printf("i=%wd shift=%wd mult=%wd\n", i, sum->group->shifts[i].n, sum->group->shifts[i].mult); */
        if (sum->group->shifts[i].n == n)
        {
            ini_i = i;
            mult = sum->group->shifts[i].mult;
            break;
        }
    }

    /* find the maximum log(x)-length of the corresponding coefficient of x on
       the rhs */

    len = acb_ode_sum_max_nlogs_xn(sum, n - sum->n0);

    /* evaluate the indicial polynomial */

    acb_ode_poly_taylor_shift_aps_trunc(ind_n, sum->ind, sum->group->leader, n,
                                        len + mult, sum->wp);

    /* flint_printf("fwd1: n=%wd mult=%wd len=%wd ind_n=%{acb_poly}\n", n, mult, len, ind_n); */

    /* advance solutions */

    next_binom_n(sum->binom_n, n, sum->nder);

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_struct * sol = sum->sol + m;
        /* XXX when rewriting with a better data structure, pass n0 - n instead
           of (n0, n) and an rhs pointer/offset separate from the sum
           pointer/offset */
        next_coeff(sol, sum->n0, n, sol->series, sol->nlogs, mult, ini_i,
                   ind_n, sum->flags & ACB_ODE_APPROX, sum->wp);
        next_sums(sol, sum->n0, n, sum->pows, sum->shifted_sums, sum->npts,
                  sum->binom_n, sum->nder, sum->sum_wp);
    }

    /* compute the next power of each evaluation point */
    /* XXX wasteful near the end */

    for (slong j = 0; j < sum->npts; j++)
    {
        slong prec = FLINT_MIN(sum->wp, sum->sum_wp);
        acb_mul(sum->pows + ((n + 1) % sum->nder) * sum->npts + j,
                sum->pows + (n       % sum->nder) * sum->npts + j,
                sum->pts + j,
                prec);
    }

    mag_mul(sum->magpow, sum->magpow, sum->mag);

    sum->n++;

    acb_poly_clear(ind_n);
}
