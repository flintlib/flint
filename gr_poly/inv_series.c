#include "arb_poly.h"
#include "gr_poly.h"

int
_gr_poly_inv_series(gr_ptr Qinv,
    gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx)
{
    int status;
    slong sz = ctx->sizeof_elem;

    Qlen = FLINT_MIN(Qlen, len);

    status = gr_inv(Qinv, Q, ctx);

    if (status != GR_SUCCESS)
        return status;

    if (Qlen == 1)
    {
        status = _gr_vec_zero(GR_ENTRY(Qinv, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status = gr_mul(GR_ENTRY(Qinv, 1, sz), Qinv, Qinv, ctx);
        status |= gr_mul(GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Q, 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(Qinv, 1, sz), GR_ENTRY(Qinv, 1, sz), ctx);
    }
    else
    {
        int is_one;
        slong i, clen, blen;

        if (Qlen <= 8)
        {
            blen = len;
        }
        else
        {
            /* todo: only use Newton when we have fast multiplication */
            blen = FLINT_MIN(10, len / 2);
        }

        /* todo: special case for -1 as well */
        is_one = (gr_is_one(Qinv, ctx) == T_TRUE);

        for (i = 1; i < blen; i++)
        {
            clen = FLINT_MIN(i, Qlen - 1);

            status |= _gr_vec_dot_rev(GR_ENTRY(Qinv, i, sz), NULL, 1,
                GR_ENTRY(Q, 1, sz), GR_ENTRY(Qinv, i - clen, sz), clen, ctx);

            if (!is_one)
                status |= gr_mul(GR_ENTRY(Qinv, i, sz), GR_ENTRY(Qinv, i, sz), Qinv, ctx);
        }

        if (len > blen)
        {
            slong Qnlen, Wlen, W2len;
            gr_ptr W;

            GR_TMP_INIT_VEC(W, len, ctx);

            NEWTON_INIT(blen, len)
            NEWTON_LOOP(m, n)

            Qnlen = FLINT_MIN(Qlen, n);
            Wlen = FLINT_MIN(Qnlen + m - 1, n);
            W2len = Wlen - m;
            status |= _gr_poly_mullow(W, Q, Qnlen, Qinv, m, Wlen, ctx);
            status |= _gr_poly_mullow(GR_ENTRY(Qinv, m, sz), Qinv, m, GR_ENTRY(W, m, sz), W2len, n - m, ctx);
            status |= _gr_vec_neg(GR_ENTRY(Qinv, m, sz), GR_ENTRY(Qinv, m, sz), n - m, ctx);

            NEWTON_END_LOOP
            NEWTON_END

            GR_TMP_CLEAR_VEC(W, len, ctx);
        }
    }

    return status;
}

int
gr_poly_inv_series(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx)
{
    int status;

    if (len == 0)
        return gr_poly_zero(Qinv, ctx);

    if (Q->length == 0)
        return GR_DOMAIN;

    if (Qinv == Q)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series(t, Q, len, ctx);
        gr_poly_swap(Qinv, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Qinv, len, ctx);
    status = _gr_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, len, ctx);
    _gr_poly_set_length(Qinv, len, ctx);
    _gr_poly_normalise(Qinv, ctx);
    return status;
}
