/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_airy_bound, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t ai, aip, bi, bip, z1, z2;
        slong prec;
        mag_t aib, aipb, bib, bipb, aim, aipm, bim, bipm;

        acb_init(ai);
        acb_init(aip);
        acb_init(bi);
        acb_init(bip);
        acb_init(z1);
        acb_init(z2);

        mag_init(aib);
        mag_init(aipb);
        mag_init(bib);
        mag_init(bipb);

        mag_init(aim);
        mag_init(aipm);
        mag_init(bim);
        mag_init(bipm);

        acb_randtest(z1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        arb_mul_ui(acb_realref(z1), acb_realref(z1), n_randint(state, 300), 1 + n_randint(state, 200));
        arb_mul_ui(acb_imagref(z1), acb_imagref(z1), n_randint(state, 300), 1 + n_randint(state, 200));

        acb_zero(z2);

        arf_set_mag(arb_midref(acb_realref(z2)), arb_radref(acb_realref(z1)));
        arf_set_mag(arb_midref(acb_imagref(z2)), arb_radref(acb_imagref(z1)));

        switch (n_randint(state, 5))
        {
            case 0:
                arf_add(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_add(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                break;
            case 1:
                arf_add(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_sub(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                break;
            case 2:
                arf_sub(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_add(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                break;
            case 3:
                arf_sub(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_sub(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
                break;
            default:
                arf_set(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)));
                arf_set(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)));
        }

        acb_hypgeom_airy_bound(aib, aipb, bib, bipb, z1);

        prec = MAG_BITS + 10;

        do {
            acb_hypgeom_airy(ai, aip, bi, bip, z2, prec);

            if (acb_rel_accuracy_bits(ai) >= MAG_BITS &&
                acb_rel_accuracy_bits(aip) >= MAG_BITS &&
                acb_rel_accuracy_bits(bi) >= MAG_BITS &&
                acb_rel_accuracy_bits(bip) >= MAG_BITS)
                break;

            prec *= 2;
        } while (1);

        acb_get_mag(aim, ai);
        acb_get_mag(aipm, aip);
        acb_get_mag(bim, bi);
        acb_get_mag(bipm, bip);

        if (mag_cmp(aim, aib) > 0 || mag_cmp(aipm, aipb) > 0 ||
            mag_cmp(bim, bib) > 0 || mag_cmp(aipm, bipb) > 0)
        {
            printf("FAIL\n");
            flint_printf("z1 = "); acb_printd(z1, 20); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 20); flint_printf("\n\n");
            flint_printf("ai = "); acb_printd(ai, 20); flint_printf("\n");
            flint_printf("aim = "); mag_printd(aim, 10); printf("\n");
            flint_printf("aib = "); mag_printd(aib, 10); printf("\n\n");
            flint_printf("api = "); acb_printd(aip, 20); flint_printf("\n");
            flint_printf("aipm = "); mag_printd(aipm, 10); printf("\n");
            flint_printf("aipb = "); mag_printd(aipb, 10); printf("\n\n");
            flint_printf("bi = "); acb_printd(bi, 20); flint_printf("\n");
            flint_printf("bim = "); mag_printd(bim, 10); printf("\n");
            flint_printf("bib = "); mag_printd(bib, 10); printf("\n\n");
            flint_printf("bpi = "); acb_printd(bip, 20); flint_printf("\n");
            flint_printf("bipm = "); mag_printd(bipm, 10); printf("\n");
            flint_printf("bipb = "); mag_printd(bipb, 10); printf("\n\n");
            flint_abort();
        }

        acb_clear(ai);
        acb_clear(aip);
        acb_clear(bi);
        acb_clear(bip);
        acb_clear(z1);
        acb_clear(z2);

        mag_clear(aib);
        mag_clear(aipb);
        mag_clear(bib);
        mag_clear(bipb);

        mag_clear(aim);
        mag_clear(aipm);
        mag_clear(bim);
        mag_clear(bipm);
    }

    TEST_FUNCTION_END(state);
}
