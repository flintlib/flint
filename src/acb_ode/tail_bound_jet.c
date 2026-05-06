#include "acb_ode.h"
#include "arb_poly.h"


static void
arb_poly_taylor_shift_trunc(arb_poly_t res, const arb_poly_t pol,
                            const arb_t a, slong ord, slong prec)
{
    /* todo: ord is typically small, evaluate derivatives... */
    arb_poly_taylor_shift(res, pol, a, prec);
    arb_poly_truncate(res, ord);
}

static void
_arb_poly_product_roots_trunc(arb_ptr res, arb_ptr rts, slong rts_len,
                              slong trunc, slong prec)
{
    FLINT_ASSERT(trunc >= 1);

    if (rts_len == 1)
        arb_neg(res, rts);
    else if (rts_len > 1)
    {
        slong mid = (rts_len + 1) / 2;
        slong len1 = FLINT_MIN(mid + 1, trunc);
        slong len2 = FLINT_MIN(rts_len - mid + 1, trunc);

        arb_ptr scratch = _arb_vec_init(len1 + len2);

        _arb_poly_product_roots_trunc(scratch, rts, mid, trunc, prec);
        _arb_poly_product_roots_trunc(scratch + len1, rts + mid, rts_len - mid,
                                      trunc, prec);
        _arb_poly_mullow(res, scratch, len1, scratch + len1, len2,
                         FLINT_MIN(rts_len, trunc), prec);

        _arb_vec_clear(scratch, len1 + len2);
    }

    if (trunc > rts_len)
        arb_one(res + rts_len);
}

static void
hexp_series(arb_poly_t res, const arb_poly_t itg_pol, slong shift,
            const arb_poly_t itg_num, const arb_poly_t itg_den,
            slong ord, slong prec)
{
    if (ord > shift + 1) {
        /* XXX Not yet tested(?)
         * 1/itg_den could be cached.
         * And couldn't we just choose pol_part_len so that this case never
         * occurs? */
        arb_poly_div_series(res, itg_num, itg_den, ord - 1 - shift, prec);
        arb_poly_shift_left(res, res, shift);
    }
    arb_poly_add_series(res, res, itg_pol, ord - 1, prec);

    arb_poly_neg(res, res);
    arb_poly_integral(res, res, prec);
    arb_poly_exp_series(res, res, ord, prec);
}

/* [M19], Algorithm 6.11, step 5,
 * specialized to w⁻¹·â(w) = itg_pol + w^itg_rat_shift·itg_num/itg_den,
 * pre = g with an implicit x^n factor.
 *
 * (To reiterate: itg_rat_shift is about the expression of \hat{f} and has
 * nothing to do with the w^{n-1} factor under the outer integral.) */

static void
maj_pre_num(arb_poly_t pre, const arb_poly_t nres_maj, slong n,
            const arb_poly_t itg_pol, slong itg_rat_shift,
            const arb_poly_t itg_num, const arb_poly_t itg_den,
            slong prec)
{
    arb_poly_t inv_exp;
    arb_poly_init(inv_exp);

    arb_poly_fit_length(pre, nres_maj->length);
    _arb_poly_set_length(pre, nres_maj->length);
    for (slong j = 0; j < pre->length; j++)
        arb_mul_si(pre->coeffs + j, nres_maj->coeffs + j, n + j, prec);

    hexp_series(inv_exp, itg_pol, itg_rat_shift, itg_num,
                itg_den, nres_maj->length, prec);

    /* polynomial factor in rhs of majorizing equation */
    arb_poly_mullow(pre, pre, inv_exp, nres_maj->length, prec);
    /* flint_printf("rhs=%{arb_poly}\n", pre); */

    for (slong j = 0; j < pre->length; j++) {
        arb_ptr c = pre->coeffs + j;
        if (arb_contains_positive(c))
            arb_div_si(c, c, n + j, prec);
        else
            arb_zero(c);
    }
    _arb_poly_normalise(pre);

    arb_poly_clear(inv_exp);
}

/* Bound the series expansion at rad of
 *
 * x^n·pre_num/den × exp ∫(itg_pol + t^itg_rat_shift·itg_num/den).
 *
 * [M19], Algorithm 8.1, specialized to this case and with a slightly different
 * output format.
 */

static void
maj_jet(arb_poly_t res, slong n, const arb_poly_t pre_num,
        const arb_poly_t itg_pol, slong itg_rat_shift, const arb_poly_t itg_num,
        const arb_t cst, arb_ptr den_rt, slong den_rt_len,
        const arb_t rad, slong ord, slong prec)
{
    /* flint_printf("== maj_jet n=%ld ord=%ld rad=%{arb} ==\n", n, ord, rad); */

    arb_poly_t pre_ser, invden_ser, shx_ser, int_pol_ser, int_rat_ser, int_ser;
    arb_ptr rt_rad = _arb_vec_init(den_rt_len);

    arb_poly_init(pre_ser);
    arb_poly_init(invden_ser);
    arb_poly_init(shx_ser);
    arb_poly_init(int_pol_ser);
    arb_poly_init(int_rat_ser);
    arb_poly_init(int_ser);

    /* Denominator: (1/den)(rad+ε) */

    for (slong i = 0; i < den_rt_len; i++)
        arb_sub(rt_rad + i, den_rt + i, rad, prec);
    arb_poly_fit_length(invden_ser, ord);
    _arb_poly_product_roots_trunc(invden_ser->coeffs,
                                  rt_rad, den_rt_len, ord, prec);
    _arb_poly_set_length(invden_ser, ord);
    _arb_poly_normalise(invden_ser);

    if (arb_contains_zero(invden_ser->coeffs))
    {
        arb_poly_fit_length(res, ord);
        _arb_vec_indeterminate(res->coeffs, ord);
        _arb_poly_set_length(res, ord);
        goto cleanup;
    }
    arb_poly_inv_series(invden_ser, invden_ser, ord, prec);
    // XXX not 100% sure this cst should be there
    arb_poly_scalar_mul(invden_ser, invden_ser, cst, prec);

    /* flint_printf("invden_ser=%{arb_poly}\n", invden_ser); */

    /* Rational part */

    arb_poly_set_coeff_si(shx_ser, 1, 1);
    arb_poly_set_coeff_arb(shx_ser, 0, rad);
    arb_poly_pow_ui_trunc_binexp(shx_ser, shx_ser, n, ord, prec);

    arb_poly_taylor_shift_trunc(pre_ser, pre_num, rad, ord, prec);
    arb_poly_mullow(pre_ser, pre_ser, shx_ser, ord, prec);
    /* pre_ser corresponds to rat_ser in ore_algebra */
    /* subtract den' from the numerator of itg_rat below instead? */
    arb_poly_mullow(pre_ser, pre_ser, invden_ser, ord, prec);

    /* flint_printf("pre_ser=%{arb_poly}\n", pre_ser); */

    /* Exponential part. Up to the sign in the exponential, this is the
     * same expression as in hexp_series, but now we want a bound on the series
     * expansion at rad instead of a series expansion at zero. For this we bound
     * ∫(f+p/q) by ∫f + (∫p)/q, then compose with rad+ε. */

    arb_poly_integral(int_pol_ser, itg_pol, prec);
    arb_poly_taylor_shift_trunc(int_pol_ser, int_pol_ser, rad, ord, prec);

    arb_poly_fit_length(int_rat_ser, itg_num->length);
    _arb_poly_set_length(int_rat_ser, itg_num->length);
    itg_rat_shift++;
    for (slong i = 0; i < itg_num->length; i++)
        arb_div_si(int_rat_ser->coeffs + i, itg_num->coeffs + i,
                   itg_rat_shift + i, prec);
    arb_poly_taylor_shift_trunc(int_rat_ser, int_rat_ser, rad, ord, prec);

    arb_poly_set_arb(shx_ser, rad);
    arb_poly_set_coeff_si(shx_ser, 1, 1);
    arb_poly_pow_ui_trunc_binexp(shx_ser, shx_ser, itg_rat_shift, ord, prec);

    arb_poly_mullow(int_rat_ser, int_rat_ser, shx_ser, ord, prec);
    arb_poly_mullow(int_rat_ser, int_rat_ser, invden_ser, ord, prec);

    /* flint_printf("int_rat_ser=%{arb_poly}\n", int_ser); */

    arb_poly_add(int_ser, int_pol_ser, int_rat_ser, prec);

    /* flint_printf("int_ser=%{arb_poly}\n", int_ser); */

    arb_poly_exp_series(res, int_ser, ord, prec);

    /* flint_printf("exp_ser=%{arb_poly}\n", res); */

    /* Final series */

    arb_poly_mullow(res, res, pre_ser, ord, prec);

cleanup:
    arb_poly_clear(pre_ser);
    arb_poly_clear(invden_ser);
    arb_poly_clear(shx_ser);
    arb_poly_clear(int_pol_ser);
    arb_poly_clear(int_rat_ser);
    arb_poly_clear(int_ser);
    _arb_vec_clear(rt_rad, den_rt_len);
}

void
acb_ode_tail_bound_jet_precomp(arb_poly_t res,
            const acb_ode_bound_t bound, slong n,
            const arb_poly_t itg_pol, const arb_poly_t itg_num,  /* n0 <= n */
            const arb_poly_t nres_maj,  /* implicit x^n factor */
            const arb_t rad, slong ord, slong prec)
{
    // flint_printf("== acb_ode_tail_bound_jet_precomp n=%ld ==\n", n);

    arb_poly_t pre;

    // flint_printf("nres_maj=%{arb_poly}\n", nres_maj);

    if (arb_poly_is_zero(nres_maj))
    {
        arb_poly_zero(res);
        return;
    }

    arb_poly_init(pre);

    /* The solution of the homogeneous part of the majorizing equation obtained
     * after specialization at n is exp(∫w⁻¹·â(w)) where
     * w⁻¹·â(w) = itg_pol + w^{pol_part_len-1} · itg_num/itg_den.
     * In other words, the itg_pol term already incorporates the w⁻¹ factor, but
     * itg_num does not. */
    slong rat_shift = bound->pol_part_len - 1;

    /* Variation of constants */
    maj_pre_num(pre, nres_maj, n, itg_pol, rat_shift, itg_num,
                bound->den_lbound, prec);

    /* Evaluate the resulting majorant */
    maj_jet(res, n, pre, itg_pol, rat_shift, itg_num, bound->cst,
            bound->den_rt, bound->den_rt_len, rad, ord, prec);

    // flint_printf("n=%ld pre=%{arb_poly}\n", n, pre);
    // flint_printf("n=%ld maj_jet=%{arb_poly}\n", n, res);

    arb_poly_clear(pre);
}

void acb_ode_tail_bound_jet(arb_poly_t res,
            const acb_ode_group_t group, const acb_ode_bound_t bound, slong n,
            slong nlogs, const arb_poly_t nres_maj, const arb_t rad,
            slong ord, slong prec)
{
    arb_poly_t itg_pol, itg_num;

    arb_poly_init(itg_pol);
    arb_poly_init(itg_num);

    acb_ode_bound_precompute_integrand(itg_pol, itg_num, group, bound, n, nlogs,
                                       prec);
    acb_ode_tail_bound_jet_precomp(res, bound, n, itg_pol, itg_num, nres_maj,
                                   rad, ord, prec);

    arb_poly_clear(itg_pol);
    arb_poly_clear(itg_num);
}
