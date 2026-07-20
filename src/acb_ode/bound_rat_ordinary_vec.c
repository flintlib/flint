#include "arf.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb_poly.h"


/* [M19], Algorithm 7.1.

   This is quite wasteful for large n, where the high-order coefficients
   contribute very little. */

void
acb_ode_bound_rat_ordinary_vec(mag_ptr res,
                               const acb_poly_struct * num, slong len,
                               const acb_poly_t ind,
                               const acb_ode_ind_lbound_t ind_lbound,
                               slong n0, slong ord, slong prec)
{
    FLINT_ASSERT(n0 >= 0);

    acb_t invcst, range_acb;
    mag_t a, invcstmag;
    arb_poly_t range_jet0, range_tail;
    acb_poly_t range_jet0_acb, range_tail_acb;
    acb_poly_t rat_jet, ind_jet, common_jet;

    acb_init(invcst);
    acb_init(range_acb);
    mag_init(a);
    mag_init(invcstmag);
    arb_poly_init(range_jet0);
    arb_poly_init(range_tail);
    acb_poly_init(range_jet0_acb);
    acb_poly_init(range_tail_acb);
    acb_poly_init(rat_jet);
    acb_poly_init(ind_jet);
    acb_poly_init(common_jet);

    arb_ptr range = acb_realref(range_acb);

    /* [n0,∞]⁻¹ */
    mag_set_ui(arb_radref(range), 2*n0);
    mag_inv(arb_radref(range), arb_radref(range));
    arf_set_mag(arb_midref(range), arb_radref(range));

    /* 1/(1+[n0,∞]⁻¹ε) */
    arb_poly_one(range_jet0);
    arb_poly_set_coeff_arb(range_jet0, 1, range);
    arb_poly_inv_series(range_jet0, range_jet0, ord, prec);
    acb_poly_set_arb_poly(range_jet0_acb, range_jet0);

    /* 1/([n0,∞]+ε) = range + range_tail + O(x^ord) */
    arb_poly_fit_length(range_tail, range_jet0->length);
    _arb_vec_scalar_mul(range_tail->coeffs + 1, range_jet0->coeffs + 1,
                        range_jet0->length - 1, range, prec);
    _arb_poly_set_length(range_tail, range_jet0->length);
    _arb_poly_normalise(range_tail);
    acb_poly_set_arb_poly(range_tail_acb, range_tail);

    /* flint_printf("n0=%ld ord=%ld range_acb=%{acb} range_jet0_acb=%{acb_poly} range_tail_acb=%{acb_poly}\n0", n0, ord, range_acb, range_jet0_acb, range_tail_acb); */

    /* 1/([n0,∞]+ε)^r · ind([n0,∞]+ε) */
    // could precompute the reciprocal polynomial (or use a shallow copy)
    slong ind_jet_len0 = FLINT_MAX(ind->length, ord);
    acb_poly_fit_length(ind_jet, ind_jet_len0);
    _acb_poly_reverse(ind_jet->coeffs, ind->coeffs, ind->length, ind->length);
    _acb_poly_set_length(ind_jet, ind->length);
    _acb_poly_normalise(ind_jet);
    /* flint_printf("rcpq_den=%{acb_poly}\n0", ind_jet); */
    /* Use the factored form of ind? But ind_lbound currently stores only a
     * subset of the roots. */
    _acb_poly_taylor_shift(ind_jet->coeffs, range_acb, ind_jet_len0, prec);
    acb_poly_compose_series(ind_jet, ind_jet, range_tail_acb, ord, prec);

    /* flint_printf("ind_jet=%{acb_poly}\n0", ind_jet); */

    /* A lower bound on |ind(k)|/k^deg(ind), where ind is the shifted indicial
       polynomial, valid for all k ≥ n with n, k ∈ ℕ ∖ exponents.
       Cf. [M19], Lemma 7.2 */
    acb_ode_ind_lbound_eval(invcstmag, ind_lbound, n0, prec);

    mag_inv(invcstmag, invcstmag);
    acb_add_error_mag(invcst, invcstmag);  /* XXX overkill? */
    _acb_vec_scalar_mul(ind_jet->coeffs + 1, ind_jet->coeffs + 1,
                        ind_jet->length - 1, invcst, prec);
    acb_one(ind_jet->coeffs);

    /* flint_printf("invcstmag=%{mag} invcst=%{acb} improved ind_jet=%{acb_poly}\n0", invcstmag, invcst, ind_jet); */

    /* 1/[n0,∞]·([n0,∞]+ε)^{r-1}/ind([n0,∞]+ε) */
    acb_poly_inv_series(common_jet, ind_jet, ord, prec);
    acb_poly_mullow(common_jet, common_jet, range_jet0_acb, ord, prec);

    for (slong i = 0; i < len; i++)
    {
        acb_poly_fit_length(rat_jet, FLINT_MAX(ind->length, ord));
        _acb_poly_reverse(rat_jet->coeffs, (num + i)->coeffs,
                          FLINT_MIN(ind->length - 1, (num + i)->length),
                          ind->length - 1);
        _acb_poly_set_length(rat_jet, ind->length - 1);
        _acb_poly_normalise(rat_jet);
        /* flint_printf("i=%ld rcpqnum=%{acb_poly} range_tail_acb=%{acb_poly}\n0", i, rat_jet, range_tail_acb); */
        _acb_poly_taylor_shift(rat_jet->coeffs, range_acb,
                               FLINT_MAX(rat_jet->length, ord), prec);
        acb_poly_compose_series(rat_jet, rat_jet, range_tail_acb, ord, prec);
        /* flint_printf("num=%{acb_poly}\n0", rat_jet); */
        acb_poly_mullow(rat_jet, rat_jet, common_jet, ord, prec);
        mag_zero(res + i);
        for (slong j = 0; j < rat_jet->length; j++)
        {
            acb_get_mag(a, rat_jet->coeffs + j);
            mag_add(res + i, res + i, a);
        }
        mag_mul(res + i, res + i, invcstmag);
    }

    acb_clear(invcst);
    acb_clear(range_acb);
    mag_clear(a);
    mag_clear(invcstmag);
    arb_poly_clear(range_jet0);
    arb_poly_clear(range_tail);
    acb_poly_clear(range_jet0_acb);
    acb_poly_clear(range_tail_acb);
    acb_poly_clear(rat_jet);
    acb_poly_clear(ind_jet);
    acb_poly_clear(common_jet);
}

