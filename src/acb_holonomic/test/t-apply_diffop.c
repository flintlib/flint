#include "test_helpers.h"
#include "acb_types.h"
#include "acb_poly.h"
#include "fmpq_types.h"
#include "fmpq_poly.h"
#include "acb_holonomic.h"

TEST_FUNCTION_START(acb_holonomic_apply_diffop, state)
{
    fmpq_poly_t tmp;
    fmpq_poly_init(tmp);

    for (slong iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        slong dop_len = n_randint(state, 3);
        acb_poly_struct * dop = _acb_poly_vec_init(dop_len);
        for (slong i = 0; i < dop_len; i++)
        {
            acb_poly_randtest(dop + i, state, 1 + n_randint(state, 10), 1 +
                              n_randint(state, 200), 5);
            /* fmpq_poly_randtest(tmp, state, 1 + n_randint(state, 5), 1 + n_randint(state, 5)); */
            /* flint_printf("dop[%wd] = %{fmpq_poly}\n", i, tmp); */
            /* acb_poly_set_fmpq_poly(dop + i, tmp, 100); */
        }

        acb_t expo;
        acb_init(expo);
        acb_randtest(expo, state, 1 + n_randint(state, 200), 5);
        slong offset = n_randint(state, 1000);

        slong logs = n_randint(state, 4);
        acb_poly_struct * f = _acb_poly_vec_init(logs);
        for (slong k = 0; k < logs; k++)
        {
            acb_poly_randtest(f + k, state, 1 + n_randint(state, 20), 1 +
                              n_randint(state, 200), 5);
            /* fmpq_poly_randtest(tmp, state, 1 + n_randint(state, 10), 1 + n_randint(state, 5)); */
            /* flint_printf("f[%wd] = %{fmpq_poly}\n", k, tmp); */
            /* acb_poly_set_fmpq_poly(f + k, tmp, 100); */
        }

        slong len1 = 1 + n_randint(state, 5);
        slong len2 = n_randint(state, 5);
        slong len = len1 + len2;
        slong start = n_randint(state, len);

        slong prec = 2 + n_randint(state, 100);

        acb_poly_struct * g1 =_acb_poly_vec_init(logs);
        acb_poly_struct * g2 =_acb_poly_vec_init(logs);
        acb_poly_struct * g3 =_acb_poly_vec_init(logs);

        acb_holonomic_apply_diffop_polmul(g1, dop, dop_len, expo, offset, f, logs, start, len, prec);

        _acb_poly_vec_fit_length(g2, logs, len1 + len2);
        /* TODO play with foff and offset too */
        _acb_holonomic_apply_diffop_polmul(g2, 0, dop, dop_len, expo, offset, f, 0, start + len1, logs, start, len1, prec);
        _acb_holonomic_apply_diffop_polmul(g2, len1, dop, dop_len, expo, offset, f, 0, start + len, logs, start + len1, len2, prec);
        _acb_poly_vec_set_length(g2, logs, len);
        _acb_poly_vec_normalise(g2, logs);

        acb_holonomic_apply_diffop_basecase(g3, dop, dop_len, expo, offset, f, logs, start, len, prec);

        if (!_acb_poly_vec_overlaps(g1, g2, logs)
                || !_acb_poly_vec_overlaps(g1, g3, logs)
                || !_acb_poly_vec_overlaps(g2, g3, logs))
        {
            flint_printf("FAIL\n\n");
            for (slong i = 0; i < dop_len; i++)
                flint_printf("dop[%wd] = %{acb_poly}\n", i, dop + i);
            for (slong k = 0; k < logs; k++)
                flint_printf("f[%wd] = %{acb_poly}\n", k, f + k);
            flint_printf("expo = %{acb}\n", expo);
            flint_printf("offset = %wd\n", offset);
            flint_printf("start = %wd\n", start);
            flint_printf("len1 = %wd\n", len1);
            flint_printf("len2 = %wd\n", len2);
            for (slong k = 0; k < logs; k++)
            {
                flint_printf("g1[%wd] = %{acb_poly}\n", k, g1 + k);
                flint_printf("g2[%wd] = %{acb_poly}\n", k, g2 + k);
                flint_printf("g3[%wd] = %{acb_poly}\n", k, g3 + k);
            }
            flint_abort();
        }

        /* multiplier à g/d par des trucs */
        /* quelques séries connues ? */

        // TODO aliasing tests

        _acb_poly_vec_clear(dop, dop_len);
        acb_clear(expo);
        _acb_poly_vec_clear(g1, logs);
        _acb_poly_vec_clear(g2, logs);
    }

    TEST_FUNCTION_END(state);
}
