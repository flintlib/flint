/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_types.h"
#include "acb_poly.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "acb_ode.h"

/* Partly written using Claude Code (Sonnet 5.5). */

TEST_FUNCTION_START(acb_ode_apply_diffop, state)
{
    for (slong iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        const char * fail = NULL;

        slong prec = 2 + n_randint(state, 100);

        gr_ctx_t CC, Pol, Ore;
        gr_ctx_init_complex_acb(CC, prec);
        gr_ctx_init_gr_poly(Pol, CC);
        gr_ore_poly_ctx_init(Ore, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

        gr_ore_poly_t p, q, pq;
        gr_ore_poly_init(p, Ore);
        gr_ore_poly_init(q, Ore);
        gr_ore_poly_init(pq, Ore);
        status |= gr_ore_poly_randtest(p, state, n_randint(state, 5), Ore);
        status |= gr_ore_poly_randtest(q, state, n_randint(state, 4), Ore);
        status |= gr_ore_poly_mul(pq, p, q, Ore);
        GR_MUST_SUCCEED(status);

        acb_poly_struct * dop = p->coeffs;
        slong dop_len = p->length;

        slong offset = n_randint(state, 2) ? n_randint(state, 10)
                                           : n_randint(state, 1000);

        /* iexpo = value of expo if it is a known integer, -1 otherwise */
        acb_t expo_buf;
        acb_init(expo_buf);
        acb_ptr expo = expo_buf;
        slong iexpo = -1;
        switch (n_randint(state, 4))
        {
            case 0:
                expo = NULL;
                iexpo = 0;
                break;
            case 1:
                iexpo = n_randint(state, 5);
                acb_set_si(expo, iexpo);
                break;
            default:
                acb_randtest(expo, state, 1 + n_randint(state, 200), 5);
        }

        slong logs = n_randint(state, 4);
        acb_poly_struct * f = _acb_poly_vec_init(logs);
        for (slong k = 0; k < logs; k++)
            acb_poly_randtest(f + k, state, 1 + n_randint(state, 20),
                              1 + n_randint(state, 200), 5);

        slong len1 = n_randint(state, 5);
        slong len2 = n_randint(state, 5);
        slong len = len1 + len2;
        slong start = n_randint(state, 10);
        slong flen = start + len;
        slong foff = n_randint(state, 3);
        slong goff = foff + flen + n_randint(state, 3);

        acb_poly_struct * g1 = _acb_poly_vec_init(logs);
        acb_poly_struct * g2 = _acb_poly_vec_init(logs);
        acb_poly_struct * g3 = _acb_poly_vec_init(logs);
        acb_poly_struct * buf0 = _acb_poly_vec_init(logs);
        acb_poly_struct * buf = _acb_poly_vec_init(logs);

        /* Non-underscore versions, with outputs not initially zero */

        for (slong k = 0; k < logs; k++)
        {
            acb_poly_randtest(g1 + k, state, n_randint(state, 30), 10, 5);
            acb_poly_randtest(g2 + k, state, n_randint(state, 30), 10, 5);
        }

        acb_ode_apply_diffop_polmul(g1, dop, dop_len, expo, offset, f, logs,
                                    start, len, prec);
        acb_ode_apply_diffop_basecase(g2, dop, dop_len, expo, offset, f, logs,
                                      start, len, prec);

        if (!_acb_poly_vec_overlaps(g1, g2, logs))
        {
            fail = "polmul vs basecase";
            goto end;
        }

        /* (p·q)(f) = p(q(f)), with aliasing */

        _acb_poly_vec_set(g2, f, logs);
        acb_ode_apply_diffop_polmul(g2, q->coeffs, q->length, expo, offset,
                                    g2, logs, 0, flen, prec);
        acb_ode_apply_diffop_basecase(g2, dop, dop_len, expo, offset,
                                      g2, logs, start, len, prec);
        if (n_randint(state, 2))
            acb_ode_apply_diffop_polmul(g3, pq->coeffs, pq->length, expo,
                                        offset, f, logs, start, len, prec);
        else
            acb_ode_apply_diffop_basecase(g3, pq->coeffs, pq->length, expo,
                                          offset, f, logs, start, len, prec);

        if (!_acb_poly_vec_overlaps(g2, g3, logs))
        {
            fail = "(p·q)(f) = p(q(f))";
            goto end;
        }

        /* Comparison with gr_ore_poly_apply */

        if (logs == 1 && iexpo >= 0)
        {
            acb_poly_shift_left(g2, f, iexpo + offset);
            status |= gr_ore_poly_apply(g2, p, g2, Ore);
            acb_poly_shift_right(g2, g2, iexpo + offset + start);
            acb_poly_truncate(g2, len);

            if (status != GR_SUCCESS || !acb_poly_overlaps(g1, g2))
            {
                fail = "gr_ore_poly_apply";
                goto end;
            }
        }

        /* Underscore versions, in place: buf0 holds f at offset foff and some
           initial values at offset goff */

        for (slong k = 0; k < logs; k++)
        {
            acb_poly_fit_length(buf0 + k, goff + len);
            _acb_poly_set_length(buf0 + k, goff + len);
            for (slong i = 0; i < goff + len; i++)
                acb_randtest((buf0 + k)->coeffs + i, state, prec, 5);
            for (slong i = 0; i < flen; i++)
                acb_poly_get_coeff_acb((buf0 + k)->coeffs + foff + i, f + k, i);

            acb_poly_shift_left(g2 + k, g1 + k, goff);
            acb_poly_add(g2 + k, g2 + k, buf0 + k, prec);
        }

        _acb_poly_vec_set(buf, buf0, logs);
        _acb_ode_apply_diffop_polmul(buf, goff, dop, dop_len, expo, offset,
                                     buf, foff, start + len1, logs,
                                     start, len1, prec);
        _acb_ode_apply_diffop_polmul(buf, goff + len1, dop, dop_len, expo,
                                     offset, buf, foff, flen, logs,
                                     start + len1, len2, prec);

        _acb_poly_vec_set(g3, buf0, logs);
        slong wlogs = logs + n_randint(state, 2);
        acb_ptr weights = _acb_vec_init(len * wlogs * flen);
        _acb_ode_apply_diffop_basecase_weights(weights, dop, dop_len, expo,
                                               offset, flen, wlogs,
                                               start, len, prec);
        _acb_ode_apply_diffop_basecase_precomp(g3, goff, weights, wlogs,
                                               g3, foff, flen, logs,
                                               start, len, prec);
        _acb_vec_clear(weights, len * wlogs * flen);

        if (!_acb_poly_vec_overlaps(buf, g2, logs)
            || !_acb_poly_vec_overlaps(g3, g2, logs))
            fail = "underscore versions";

end:
        if (fail != NULL)
        {
            flint_printf("FAIL: %s\n\n", fail);
            flint_printf("p = ");
            status |= gr_ore_poly_print(p, Ore);
            flint_printf("\nq = ");
            status |= gr_ore_poly_print(q, Ore);
            flint_printf("\nexpo = %{acb}\n", expo_buf);
            flint_printf("offset = %wd, start = %wd, len1 = %wd, len2 = %wd, "
                         "foff = %wd, goff = %wd\n",
                         offset, start, len1, len2, foff, goff);
            for (slong k = 0; k < logs; k++)
                flint_printf("f[%wd] = %{acb_poly}\n\n"
                             "g1[%wd] = %{acb_poly}\n\n"
                             "g2[%wd] = %{acb_poly}\n\n"
                             "g3[%wd] = %{acb_poly}\n\n"
                             "buf[%wd] = %{acb_poly}\n\n",
                             k, f + k, k, g1 + k, k, g2 + k, k, g3 + k,
                             k, buf + k);
            flint_abort();
        }

        _acb_poly_vec_clear(f, logs);
        _acb_poly_vec_clear(g1, logs);
        _acb_poly_vec_clear(g2, logs);
        _acb_poly_vec_clear(g3, logs);
        _acb_poly_vec_clear(buf0, logs);
        _acb_poly_vec_clear(buf, logs);
        acb_clear(expo_buf);
        gr_ore_poly_clear(p, Ore);
        gr_ore_poly_clear(q, Ore);
        gr_ore_poly_clear(pq, Ore);
        gr_ctx_clear(Ore);
        gr_ctx_clear(Pol);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
