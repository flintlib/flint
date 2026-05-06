#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb_poly.h"
#include "fmpz_mat.h"
#include "gr_mat.h"


/* Write dop ∈ K[x][θ] as rcoeff[0] + θ·rcoeff[1] + ··· */

static void
dop_rcoeffs(acb_poly_struct * rcoeff, const acb_poly_struct * dop, slong len,
            slong prec)
{
    fmpz_t c, pow;
    fmpz_mat_t binom;

    fmpz_init(c);
    fmpz_init(pow);
    fmpz_mat_init(binom, len, len);

    gr_ctx_t fmpz_ctx;
    gr_ctx_init_fmpz(fmpz_ctx);
    GR_MUST_SUCCEED(gr_mat_pascal((gr_mat_struct *) binom, -1, fmpz_ctx));
    gr_ctx_clear(fmpz_ctx);

    slong xlen = 0;
    for (slong i = 0; i < len; i++)
        xlen = FLINT_MAX(xlen, (dop + i)->length);

    for (slong k = 0; k < len; k++)
    {
        acb_poly_struct * pol = rcoeff + k;
        acb_poly_zero(pol);
        acb_poly_fit_length(pol, xlen);

        for (slong j = 0; j < xlen; j++)
        {
            fmpz_one(pow);
            for (slong i = 0; i < len - k; i++)
            {
                if ((dop + k + i)->length > j)
                {
                    fmpz_mul(c, pow, fmpz_mat_entry(binom, k + i, i));
                    acb_addmul_fmpz(pol->coeffs + j, (dop + k + i)->coeffs + j,
                                    c, prec);
                }
                if (j == 0)
                    break;
                fmpz_mul_si(pow, pow, -j);
            }
        }

        _acb_poly_set_length(pol, xlen);
        _acb_poly_normalise(pol);
    }

    fmpz_mat_clear(binom);
    fmpz_clear(c);
    fmpz_clear(pow);
}

/* Write (rcoeffs[0] + θ·rcoeffs[1] + ··· + θ^{dop_len-1}·lc)·lc⁻¹ as
       initial[0] + θ·initial[1] + ···
         + (rem_num[0] + θ·rem_num[1] + ···)·x^trunc·lc⁻¹
   where the polynomials initial[i] have degree < trunc.
   We assume trunc >= 1, so that rem_num has length < dop_len. */

static void
split_dop(acb_poly_struct * initial, acb_poly_struct * rem_num,
          const acb_poly_struct * rcoeffs, slong dop_len,
          slong trunc, slong prec)
{
    acb_poly_t inv, rev;

    FLINT_ASSERT(trunc >= 1);

    acb_poly_init(inv);
    acb_poly_init(rev);

    /* todo: when deg(lc) < trunc, we should only invert to order ~deg(lc)... */
    const acb_poly_struct * lc = rcoeffs + dop_len - 1;
    acb_poly_inv_series(inv, lc, trunc, prec);

    for (slong i = 0; i < dop_len - 1; i++)
    {
        acb_poly_mullow(initial + i, rcoeffs + i, inv, trunc, prec);
        acb_poly_mul(rem_num + i, initial + i, lc, prec);
        acb_poly_sub(rem_num + i, rcoeffs + i, rem_num + i, prec);
        /* drop known zeroes (todo: avoid computing them) */
        acb_poly_shift_right(rem_num + i, rem_num + i, trunc);
    }
    acb_poly_one(initial + dop_len - 1);

    acb_poly_clear(inv);
    acb_poly_clear(rev);
}

static void
swap_vars(acb_poly_struct * b, acb_poly_struct * a, slong len)
{
    for (slong i = 0; i < len; i++)
    {
        for (slong j = 0; j < (a + i)->length; j++)
        {
            acb_poly_fit_length(b + j, i + 1);
            acb_swap((b + j)->coeffs + i, (a + i)->coeffs + j);
            _acb_poly_set_length(b + j, i + 1);
            _acb_poly_normalise(b + j);
        }
    }
}

void
acb_ode_bound_precompute(acb_ode_bound_t bound,
                         const acb_poly_struct * dop, slong dop_len,
                         acb_srcptr lcrt,
                         slong pol_part_len, slong prec)
{
    const acb_poly_struct * den = dop + dop_len - 1;

    /* Denominator bound: 1/lc_Θ(dop) << cst/prod(rt-x) = cst/den_lbound */

    if (prec > bound->prec)
    {
        bound->prec = prec;

        acb_abs(bound->cst, den->coeffs + den->length - 1, prec);
        arb_inv(bound->cst, bound->cst, prec);

        slong rtlen = den->length - 1;
        if (bound->den_rt != NULL)
            _arb_vec_clear(bound->den_rt, bound->den_rt_len);
        bound->den_rt = _arb_vec_init(rtlen);
        bound->den_rt_len = rtlen;

        for (slong j = 0; j < rtlen; j++)
            acb_get_abs_lbound_arf(arb_midref(bound->den_rt + j),
                                   lcrt + j, prec);

        arb_poly_product_roots(bound->den_lbound,
                               bound->den_rt, den->length - 1, prec);
    }

    /* Write dop with θ on the left, right-divide by lc, isolate pol_part_len
       initial terms in the series expansion wrt x of the resulting fractions.

       We have order(initial) = order(dop), order(rem_num) < order(dop) (because
       the leading coefficient of dop/lc is 1 and pol_part_len >= 1).

       todo: Avoid recomputations when refining a bound computed for a different
       pol_part_len. */

    if (pol_part_len <= 0)
        pol_part_len = 1;

    if (prec > bound->prec || pol_part_len > bound->pol_part_len)
    {
        acb_poly_struct * rcoeffs = _acb_poly_vec_init(3*dop_len - 1);
        acb_poly_struct * initial = rcoeffs + dop_len;
        acb_poly_struct * rem_num = rcoeffs + 2*dop_len;

        dop_rcoeffs(rcoeffs, dop, dop_len, prec);

        /* Our pol_part_len = (that of ore_algebra) + 1, see comment below */
        split_dop(initial, rem_num, rcoeffs, dop_len, pol_part_len, prec);

        /* Rewrite the truncated part and the numerator of the remainder as
           elements of K[θ][x], which may also be viewed as elements of
           K[n][S⁻¹].

           Note: here initial_n[0] = indicial polynomial, whereas in ore_algebra
           this first element is dropped (so that only those that contribute to
           the integrand in Algo 6.11, step 3 of [M19] remain) and all indices
           are shifted by one. */

        slong rem_num_clen = 0;
        for (slong i = 0; i < dop_len - 1; i++)
            rem_num_clen = FLINT_MAX(rem_num_clen, (rem_num + i)->length);

        if (bound->all_nums != NULL)
            _acb_poly_vec_clear(bound->all_nums, bound->all_nums_len);
        bound->pol_part_len = pol_part_len;
        bound->all_nums_len = pol_part_len + rem_num_clen;
        bound->all_nums = _acb_poly_vec_init(bound->all_nums_len);

        acb_poly_struct * initial_n = bound->all_nums;
        acb_poly_struct * rem_num_n = bound->all_nums + pol_part_len;
        swap_vars(initial_n, initial, dop_len);
        swap_vars(rem_num_n, rem_num, dop_len - 1);

        _acb_poly_vec_clear(rcoeffs, 3*dop_len - 1);
    }
}
