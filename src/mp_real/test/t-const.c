/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* All constants through all four paths: the ball from scratch and from
   the cache, and the limb floors from scratch and from the cache.  The
   two floors must be identical and exact, *err = 1, the floor must lie
   in the from-scratch ball, and both balls must contain the 160-digit
   references (independent of FLINT: computed with mpmath) as far as
   those reach.  The cache is grown in random orders and cleared at
   random, and the tests run on up to four threads, each with its own
   cache. */

typedef void (* const_limb_func)(nn_ptr, ulong *, slong, int);

static const mp_real_const_func const_ball[] = {
    mp_real_const_pi4, mp_real_const_log2, mp_real_const_euler,
    mp_real_const_e, mp_real_const_log10, mp_real_const_catalan,
    mp_real_const_zeta3, mp_real_const_zeta5, mp_real_const_gamma_1_3,
    mp_real_const_gamma_1_4, mp_real_const_2_div_pi,
};

static const const_limb_func const_limb[] = {
    _mp_real_const_pi4, _mp_real_const_log2, _mp_real_const_euler,
    _mp_real_const_e, _mp_real_const_log10, _mp_real_const_catalan,
    _mp_real_const_zeta3, _mp_real_const_zeta5, _mp_real_const_gamma_1_3,
    _mp_real_const_gamma_1_4, _mp_real_const_2_div_pi,
};

/* a units limb for the constants in [1, B) */
static const int const_units[] = { 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0 };

static const char * const_name[] = { "pi/4", "log(2)", "euler", "e",
    "log(10)", "catalan", "zeta(3)", "zeta(5)", "gamma(1/3)",
    "gamma(1/4)", "2/pi" };

static const char * const_ref[] = {
    "0.7853981633974483096156608458198757210492923498437764552437361480769541015715522496570087063355292669955370216283205766617734611523876455579313398520321202793626",
    "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516350759619307",
    "0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495146314472498070824809605040144865428362241739976449235362535",
    "2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738",
    "2.302585092994045684017991454684364207601101488628772976033327900967572609677352480235997205089598298341967784042286248633409525465082806756666287369098781689483",
    "0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694793565129261151062485744226191961995790358988033258590594315947374811584069953",
    "1.202056903159594285399738161511449990764986292340498881792271555341838205786313090186455873609335258146199157795260719418491995998673283213776396837207900161454",
    "1.036927755143369926331365486457034168057080919501912811974192677903803589786281484560043106557133336379620341466556609042800961779155970841835110721800876448663",
    "2.678938534707747633655692940974677644128689377957301100950428327590417610167743819540982889041188789419159049200072263335719084569504472259977713367708469768167",
    "3.625609908221908311930685155867672002995167682880065467433377999569919243538729121618360136723384300361471751392420719965891524094022559977426458890361450606414",
    "0.6366197723675813430755350534900574481378385829618257949906693762355871905369061403604552110650123438242913709070318321475716473844583146115118696429267993569169",
};

#define NUM_CONST 11

TEST_FUNCTION_START(mp_real_const, state)
{
    slong iter, i;
    arb_t ref[NUM_CONST];

    for (i = 0; i < NUM_CONST; i++)
    {
        arb_init(ref[i]);
        arb_set_str(ref[i], const_ref[i], 560);
        /* 160 digits, the last one rounded */
        arb_add_error_2exp_si(ref[i], -525);
    }

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong n, units, len;
        int c;
        ulong err0, err1;
        nn_ptr y0, y1;
        mp_real_t b0, b1;
        arb_t a0, a1, fl;
        fmpz_t f;

        c = n_randint(state, NUM_CONST);
        n = 1 + n_randint(state, (iter % 20 == 0) ? 400 : 24);
        units = const_units[c];
        len = n + units;

        if (n_randint(state, 8) == 0)
            _mp_real_const_clear_cache();
        if (n_randint(state, 4) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));

        y0 = flint_malloc(len * sizeof(ulong));
        y1 = flint_malloc(len * sizeof(ulong));
        mp_real_init(b0);
        mp_real_init(b1);
        arb_init(a0);
        arb_init(a1);
        arb_init(fl);
        fmpz_init(f);

        const_ball[c](b0, n, 0);
        const_ball[c](b1, n, 1);
        err0 = err1 = 0;
        const_limb[c](y0, &err0, n, 0);
        const_limb[c](y1, &err1, n, 1);
        const_limb[c](y1, NULL, n, 1);   /* err may be NULL */

        mp_real_get_arb(a0, b0);
        mp_real_get_arb(a1, b1);

        if (err0 != 1 || err1 != 1)
            TEST_FUNCTION_FAIL("%s: err = %wu, %wu\n", const_name[c], err0, err1);

        if (mpn_cmp(y0, y1, len) != 0)
            TEST_FUNCTION_FAIL("%s: floors differ, n = %wd\n", const_name[c], n);

        /* the floor f B^-n satisfies f B^-n <= c < (f + 1) B^-n: the
           interval [f, f + 1] B^-n must meet the from-scratch ball and
           lie within one ulp of the reference */
        fmpz_set_ui_array(f, y0, len);
        arb_set_fmpz(fl, f);
        arb_add_error_2exp_si(fl, 0);
        arb_mul_2exp_si(fl, fl, -FLINT_BITS * n);
        if (!arb_overlaps(fl, a0))
            TEST_FUNCTION_FAIL("%s: floor outside the ball, n = %wd\n", const_name[c], n);
        if (!arb_overlaps(fl, ref[c]) || !arb_overlaps(a0, ref[c]) || !arb_overlaps(a1, ref[c]))
            TEST_FUNCTION_FAIL("%s: reference mismatch, n = %wd\n", const_name[c], n);
        if (!arb_overlaps(a0, a1))
            TEST_FUNCTION_FAIL("%s: the balls differ, n = %wd\n", const_name[c], n);

        /* the cached ball is the floor with one ulp of radius */
        if (arb_rel_accuracy_bits(a1) < FLINT_BITS * n - 8)
            TEST_FUNCTION_FAIL("%s: cached ball accuracy %wd, n = %wd\n",
                const_name[c], arb_rel_accuracy_bits(a1), n);

        flint_free(y0);
        flint_free(y1);
        mp_real_clear(b0);
        mp_real_clear(b1);
        arb_clear(a0);
        arb_clear(a1);
        arb_clear(fl);
        fmpz_clear(f);
    }

    for (i = 0; i < NUM_CONST; i++)
        arb_clear(ref[i]);

    _mp_real_const_clear_cache();
    flint_set_num_threads(1);

    TEST_FUNCTION_END(state);
}
