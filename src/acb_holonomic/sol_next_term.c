#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_holonomic.h"


typedef enum
{
    LC_UNKNOWN,
    LC_ACB,
    LC_FMPZ
}
lc_status_t;


/* XXX take off = n - base instead of (base, n)?
 * XXX separate rhs_off from off? */


static void
next_coeff(acb_holonomic_sol_struct * sol, slong base, slong n,
           acb_poly_struct * rhs, slong rhs_nlogs,
           slong mult, slong ini_i, const acb_poly_struct * ind_n,
           slong prec)
{
    fmpz_t lc_fmpz;
    acb_t invlc;
    slong k;

    lc_status_t lc_status = LC_UNKNOWN;
    fmpz_init(lc_fmpz);
    acb_init(invlc);

    while (rhs_nlogs > 0
           && acb_is_zero((rhs + rhs_nlogs - 1)->coeffs + n - base))
         rhs_nlogs--;

   /* Compute the high-log part of the new term from the coefficient of x^{λ+n}
    * on the rhs and the indicial polynomial */

    acb_struct * new_term = _acb_vec_init(rhs_nlogs);

    for (k = rhs_nlogs - 1; k >= 0; k--)
    {
        /* combin = rhs[k][0] + sum(ind_n[mult+1+j] * new_term[k+1+j], j≥0) */
        acb_dot(new_term + k,
                (rhs + k)->coeffs + n - base,
                0,
                ind_n->coeffs + mult + 1, 1,
                new_term + k + 1, 1,
                rhs_nlogs - k - 1,
                prec);

        if (acb_is_zero(new_term + k))
            continue;

        /* only compute invlc if needed */
        if (lc_status == 0)
        {
            if (acb_is_exact(ind_n->coeffs + mult)
                && acb_get_unique_fmpz(lc_fmpz, ind_n->coeffs + mult))
            {
                lc_status = LC_FMPZ;
                acb_indeterminate(invlc);
            }
            else
            {
                lc_status = LC_ACB;
                acb_inv(invlc, ind_n->coeffs + mult, prec);
            }
        }

        if (lc_status == LC_FMPZ)
            acb_div_fmpz(new_term + k, new_term + k, lc_fmpz, prec);
        else
            acb_mul(new_term + k, new_term + k, invlc, prec);

        acb_neg(new_term + k, new_term + k);
    }

    /* Write the new term to sol, store new extini */

    for (k = 0; k < mult; k++)  /* initial conditions */
    {
        acb_set((sol->series + k)->coeffs +  n - base,
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

    _acb_vec_clear(new_term, rhs_nlogs);
    acb_clear(invlc);
    fmpz_clear(lc_fmpz);
}


static void
next_sums(acb_holonomic_sol_struct * sol, slong base, slong n,
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
            acb_ptr s = acb_holonomic_sol_sum_ptr(sol, j, k)->coeffs;

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
_acb_holonomic_sol_add_term(
        acb_holonomic_sol_struct * sol, slong base, slong n,
        acb_poly_struct * rhs, slong rhs_nlogs,
        slong mult, slong ini_i, const acb_poly_struct * ind_n,
        acb_srcptr pows, const char * shifted, slong npts,
        const fmpz * binom_n,
        slong nder, slong prec, slong sums_prec)
{
    next_coeff(sol, base, n, rhs, rhs_nlogs, mult, ini_i, ind_n, prec);
    next_sums(sol, base, n, pows, shifted, npts, binom_n, nder, prec);
}
