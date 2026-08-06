/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "padic_radix.h"
#include "gr.h"
#include "nmod_poly.h"

/* Compressed table of lcm(1, ..., n) / i for small n. */

#if FLINT_BITS == 64
#define LOG_LCMTAB_MAXN 46
static const ulong log_lcm_data[449] = {
    UWORD(0), UWORD(1), UWORD(0), UWORD(2), UWORD(1), UWORD(0),
    UWORD(6), UWORD(3), UWORD(2), UWORD(0), UWORD(12), UWORD(6),
    UWORD(4), UWORD(3), UWORD(0), UWORD(60), UWORD(30), UWORD(20),
    UWORD(15), UWORD(12), UWORD(10), UWORD(0), UWORD(420), UWORD(210),
    UWORD(140), UWORD(105), UWORD(84), UWORD(70), UWORD(60), UWORD(0),
    UWORD(840), UWORD(420), UWORD(280), UWORD(210), UWORD(168), UWORD(140),
    UWORD(120), UWORD(105), UWORD(0), UWORD(2520), UWORD(1260), UWORD(840),
    UWORD(630), UWORD(504), UWORD(420), UWORD(360), UWORD(315), UWORD(280),
    UWORD(252), UWORD(0), UWORD(27720), UWORD(13860), UWORD(9240), UWORD(6930),
    UWORD(5544), UWORD(4620), UWORD(3960), UWORD(3465), UWORD(3080), UWORD(2772),
    UWORD(2520), UWORD(2310), UWORD(0), UWORD(360360), UWORD(180180), UWORD(120120),
    UWORD(90090), UWORD(72072), UWORD(60060), UWORD(51480), UWORD(45045), UWORD(40040),
    UWORD(36036), UWORD(32760), UWORD(30030), UWORD(27720), UWORD(25740), UWORD(24024),
    UWORD(0), UWORD(720720), UWORD(360360), UWORD(240240), UWORD(180180), UWORD(144144),
    UWORD(120120), UWORD(102960), UWORD(90090), UWORD(80080), UWORD(72072), UWORD(65520),
    UWORD(60060), UWORD(55440), UWORD(51480), UWORD(48048), UWORD(45045), UWORD(0),
    UWORD(12252240), UWORD(6126120), UWORD(4084080), UWORD(3063060), UWORD(2450448), UWORD(2042040),
    UWORD(1750320), UWORD(1531530), UWORD(1361360), UWORD(1225224), UWORD(1113840), UWORD(1021020),
    UWORD(942480), UWORD(875160), UWORD(816816), UWORD(765765), UWORD(720720), UWORD(680680),
    UWORD(0), UWORD(232792560), UWORD(116396280), UWORD(77597520), UWORD(58198140), UWORD(46558512),
    UWORD(38798760), UWORD(33256080), UWORD(29099070), UWORD(25865840), UWORD(23279256), UWORD(21162960),
    UWORD(19399380), UWORD(17907120), UWORD(16628040), UWORD(15519504), UWORD(14549535), UWORD(13693680),
    UWORD(12932920), UWORD(12252240), UWORD(11639628), UWORD(11085360), UWORD(10581480), UWORD(0),
    UWORD(5354228880), UWORD(2677114440), UWORD(1784742960), UWORD(1338557220), UWORD(1070845776), UWORD(892371480),
    UWORD(764889840), UWORD(669278610), UWORD(594914320), UWORD(535422888), UWORD(486748080), UWORD(446185740),
    UWORD(411863760), UWORD(382444920), UWORD(356948592), UWORD(334639305), UWORD(314954640), UWORD(297457160),
    UWORD(281801520), UWORD(267711444), UWORD(254963280), UWORD(243374040), UWORD(232792560), UWORD(223092870),
    UWORD(0), UWORD(26771144400), UWORD(13385572200), UWORD(8923714800), UWORD(6692786100), UWORD(5354228880),
    UWORD(4461857400), UWORD(3824449200), UWORD(3346393050), UWORD(2974571600), UWORD(2677114440), UWORD(2433740400),
    UWORD(2230928700), UWORD(2059318800), UWORD(1912224600), UWORD(1784742960), UWORD(1673196525), UWORD(1574773200),
    UWORD(1487285800), UWORD(1409007600), UWORD(1338557220), UWORD(1274816400), UWORD(1216870200), UWORD(1163962800),
    UWORD(1115464350), UWORD(1070845776), UWORD(1029659400), UWORD(0), UWORD(80313433200), UWORD(40156716600),
    UWORD(26771144400), UWORD(20078358300), UWORD(16062686640), UWORD(13385572200), UWORD(11473347600), UWORD(10039179150),
    UWORD(8923714800), UWORD(8031343320), UWORD(7301221200), UWORD(6692786100), UWORD(6177956400), UWORD(5736673800),
    UWORD(5354228880), UWORD(5019589575), UWORD(4724319600), UWORD(4461857400), UWORD(4227022800), UWORD(4015671660),
    UWORD(3824449200), UWORD(3650610600), UWORD(3491888400), UWORD(3346393050), UWORD(3212537328), UWORD(3088978200),
    UWORD(2974571600), UWORD(2868336900), UWORD(0), UWORD(2329089562800), UWORD(1164544781400), UWORD(776363187600),
    UWORD(582272390700), UWORD(465817912560), UWORD(388181593800), UWORD(332727080400), UWORD(291136195350), UWORD(258787729200),
    UWORD(232908956280), UWORD(211735414800), UWORD(194090796900), UWORD(179160735600), UWORD(166363540200), UWORD(155272637520),
    UWORD(145568097675), UWORD(137005268400), UWORD(129393864600), UWORD(122583661200), UWORD(116454478140), UWORD(110909026800),
    UWORD(105867707400), UWORD(101264763600), UWORD(97045398450), UWORD(93163582512), UWORD(89580367800), UWORD(86262576400),
    UWORD(83181770100), UWORD(80313433200), UWORD(77636318760), UWORD(0), UWORD(72201776446800), UWORD(36100888223400),
    UWORD(24067258815600), UWORD(18050444111700), UWORD(14440355289360), UWORD(12033629407800), UWORD(10314539492400), UWORD(9025222055850),
    UWORD(8022419605200), UWORD(7220177644680), UWORD(6563797858800), UWORD(6016814703900), UWORD(5553982803600), UWORD(5157269746200),
    UWORD(4813451763120), UWORD(4512611027925), UWORD(4247163320400), UWORD(4011209802600), UWORD(3800093497200), UWORD(3610088822340),
    UWORD(3438179830800), UWORD(3281898929400), UWORD(3139207671600), UWORD(3008407351950), UWORD(2888071057872), UWORD(2776991401800),
    UWORD(2674139868400), UWORD(2578634873100), UWORD(2489716429200), UWORD(2406725881560), UWORD(2329089562800), UWORD(0),
    UWORD(144403552893600), UWORD(72201776446800), UWORD(48134517631200), UWORD(36100888223400), UWORD(28880710578720), UWORD(24067258815600),
    UWORD(20629078984800), UWORD(18050444111700), UWORD(16044839210400), UWORD(14440355289360), UWORD(13127595717600), UWORD(12033629407800),
    UWORD(11107965607200), UWORD(10314539492400), UWORD(9626903526240), UWORD(9025222055850), UWORD(8494326640800), UWORD(8022419605200),
    UWORD(7600186994400), UWORD(7220177644680), UWORD(6876359661600), UWORD(6563797858800), UWORD(6278415343200), UWORD(6016814703900),
    UWORD(5776142115744), UWORD(5553982803600), UWORD(5348279736800), UWORD(5157269746200), UWORD(4979432858400), UWORD(4813451763120),
    UWORD(4658179125600), UWORD(4512611027925), UWORD(4375865239200), UWORD(4247163320400), UWORD(4125815796960), UWORD(4011209802600),
    UWORD(0), UWORD(5342931457063200), UWORD(2671465728531600), UWORD(1780977152354400), UWORD(1335732864265800), UWORD(1068586291412640),
    UWORD(890488576177200), UWORD(763275922437600), UWORD(667866432132900), UWORD(593659050784800), UWORD(534293145706320), UWORD(485721041551200),
    UWORD(445244288088600), UWORD(410994727466400), UWORD(381637961218800), UWORD(356195430470880), UWORD(333933216066450), UWORD(314290085709600),
    UWORD(296829525392400), UWORD(281206918792800), UWORD(267146572853160), UWORD(254425307479200), UWORD(242860520775600), UWORD(232301367698400),
    UWORD(222622144044300), UWORD(213717258282528), UWORD(205497363733200), UWORD(197886350261600), UWORD(190818980609400), UWORD(184239015760800),
    UWORD(178097715235440), UWORD(172352627647200), UWORD(166966608033225), UWORD(161907013850400), UWORD(157145042854800), UWORD(152655184487520),
    UWORD(148414762696200), UWORD(144403552893600), UWORD(140603459396400), UWORD(136998242488800), UWORD(133573286426580), UWORD(0),
    UWORD(219060189739591200), UWORD(109530094869795600), UWORD(73020063246530400), UWORD(54765047434897800), UWORD(43812037947918240), UWORD(36510031623265200),
    UWORD(31294312819941600), UWORD(27382523717448900), UWORD(24340021082176800), UWORD(21906018973959120), UWORD(19914562703599200), UWORD(18255015811632600),
    UWORD(16850783826122400), UWORD(15647156409970800), UWORD(14604012649306080), UWORD(13691261858724450), UWORD(12885893514093600), UWORD(12170010541088400),
    UWORD(11529483670504800), UWORD(10953009486979560), UWORD(10431437606647200), UWORD(9957281351799600), UWORD(9524356075634400), UWORD(9127507905816300),
    UWORD(8762407589583648), UWORD(8425391913061200), UWORD(8113340360725600), UWORD(7823578204985400), UWORD(7553799646192800), UWORD(7302006324653040),
    UWORD(7066457733535200), UWORD(6845630929362225), UWORD(6638187567866400), UWORD(6442946757046800), UWORD(6258862563988320), UWORD(6085005270544200),
    UWORD(5920545668637600), UWORD(5764741835252400), UWORD(5616927942040800), UWORD(5476504743489780), UWORD(5342931457063200), UWORD(5215718803323600),
    UWORD(0), UWORD(9419588158802421600), UWORD(4709794079401210800), UWORD(3139862719600807200), UWORD(2354897039700605400), UWORD(1883917631760484320),
    UWORD(1569931359800403600), UWORD(1345655451257488800), UWORD(1177448519850302700), UWORD(1046620906533602400), UWORD(941958815880242160), UWORD(856326196254765600),
    UWORD(784965679900201800), UWORD(724583704523263200), UWORD(672827725628744400), UWORD(627972543920161440), UWORD(588724259925151350), UWORD(554093421106024800),
    UWORD(523310453266801200), UWORD(495767797831706400), UWORD(470979407940121080), UWORD(448551817085829600), UWORD(428163098127382800), UWORD(409547311252279200),
    UWORD(392482839950100900), UWORD(376783526352096864), UWORD(362291852261631600), UWORD(348873635511200800), UWORD(336413862814372200), UWORD(324813384786290400),
    UWORD(313986271960080720), UWORD(303857682542013600), UWORD(294362129962575675), UWORD(285442065418255200), UWORD(277046710553012400), UWORD(269131090251497760),
    UWORD(261655226633400600), UWORD(254583463751416800), UWORD(247883898915853200), UWORD(241527901507754400), UWORD(235489703970060540), UWORD(229746052653717600),
    UWORD(224275908542914800), UWORD(219060189739591200), UWORD(214081549063691400), UWORD(209324181306720480), UWORD(204773655626139600),
};
static const unsigned short log_lcm_off[LOG_LCMTAB_MAXN + 1] = { 0, 0, 2, 5, 9, 14, 14, 21, 29, 38, 38, 49, 49, 62, 62, 62, 78, 95, 95, 114, 114, 114, 114, 137, 137, 162, 162, 189, 189, 218, 218, 249, 281, 281, 281, 281, 281, 318, 318, 318, 318, 359, 359, 402, 402, 402, 402 };
#endif

#if FLINT_BITS == 32
#define LOG_LCMTAB_MAXN 22
static const ulong log_lcm_data[137] = {
    UWORD(0), UWORD(1), UWORD(0), UWORD(2), UWORD(1), UWORD(0),
    UWORD(6), UWORD(3), UWORD(2), UWORD(0), UWORD(12), UWORD(6),
    UWORD(4), UWORD(3), UWORD(0), UWORD(60), UWORD(30), UWORD(20),
    UWORD(15), UWORD(12), UWORD(10), UWORD(0), UWORD(420), UWORD(210),
    UWORD(140), UWORD(105), UWORD(84), UWORD(70), UWORD(60), UWORD(0),
    UWORD(840), UWORD(420), UWORD(280), UWORD(210), UWORD(168), UWORD(140),
    UWORD(120), UWORD(105), UWORD(0), UWORD(2520), UWORD(1260), UWORD(840),
    UWORD(630), UWORD(504), UWORD(420), UWORD(360), UWORD(315), UWORD(280),
    UWORD(252), UWORD(0), UWORD(27720), UWORD(13860), UWORD(9240), UWORD(6930),
    UWORD(5544), UWORD(4620), UWORD(3960), UWORD(3465), UWORD(3080), UWORD(2772),
    UWORD(2520), UWORD(2310), UWORD(0), UWORD(360360), UWORD(180180), UWORD(120120),
    UWORD(90090), UWORD(72072), UWORD(60060), UWORD(51480), UWORD(45045), UWORD(40040),
    UWORD(36036), UWORD(32760), UWORD(30030), UWORD(27720), UWORD(25740), UWORD(24024),
    UWORD(0), UWORD(720720), UWORD(360360), UWORD(240240), UWORD(180180), UWORD(144144),
    UWORD(120120), UWORD(102960), UWORD(90090), UWORD(80080), UWORD(72072), UWORD(65520),
    UWORD(60060), UWORD(55440), UWORD(51480), UWORD(48048), UWORD(45045), UWORD(0),
    UWORD(12252240), UWORD(6126120), UWORD(4084080), UWORD(3063060), UWORD(2450448), UWORD(2042040),
    UWORD(1750320), UWORD(1531530), UWORD(1361360), UWORD(1225224), UWORD(1113840), UWORD(1021020),
    UWORD(942480), UWORD(875160), UWORD(816816), UWORD(765765), UWORD(720720), UWORD(680680),
    UWORD(0), UWORD(232792560), UWORD(116396280), UWORD(77597520), UWORD(58198140), UWORD(46558512),
    UWORD(38798760), UWORD(33256080), UWORD(29099070), UWORD(25865840), UWORD(23279256), UWORD(21162960),
    UWORD(19399380), UWORD(17907120), UWORD(16628040), UWORD(15519504), UWORD(14549535), UWORD(13693680),
    UWORD(12932920), UWORD(12252240), UWORD(11639628), UWORD(11085360), UWORD(10581480),
};
static const unsigned short log_lcm_off[LOG_LCMTAB_MAXN + 1] = { 0, 0, 2, 5, 9, 14, 14, 21, 29, 38, 38, 49, 49, 62, 62, 62, 78, 95, 95, 114, 114, 114, 114 };
#endif

/* a <- a / p^d  (a divisible by p^d), in place on an l-limb array mod B^l. */
static void
_arr_rshift_digits(nn_ptr a, slong l, slong d, const radix_t radix)
{
    slong e = radix->exp;
    slong dl = d / e, t;
    unsigned int dr = (unsigned int) (d % e);

    if (dl > 0)
    {
        for (t = 0; t + dl < l; t++)
            a[t] = a[t + dl];
        for (; t < l; t++)
            a[t] = 0;
    }
    if (dr > 0)
        radix_rshift_digits(a, a, l, dr, radix);
}

/* a <- (a * p^d) mod B^l, in place on an l-limb array. */
static void
_arr_lshift_digits(nn_ptr a, slong l, slong d, const radix_t radix)
{
    slong e = radix->exp;
    slong dl = d / e, t;
    unsigned int dr = (unsigned int) (d % e);

    if (dl > 0)
    {
        for (t = l - 1; t >= dl; t--)
            a[t] = a[t - dl];
        for (; t >= 0; t--)
            a[t] = 0;
    }
    if (dr > 0)
        radix_lshift_digits(a, a, l, dr, radix);
}

/* pows[0] = 1, pows[1] = y mod B^l, pows[i] = y^i mod B^l for i = 2, ..., b. */
static void
_build_pows(nn_ptr * pows, slong b, const radix_integer_t y, slong l,
    const radix_t radix)
{
    slong e = radix->exp, i, t;
    radix_integer_t xi;

    radix_integer_init(xi, radix);
    radix_integer_mod_digits(xi, y, e * l, radix);
    for (t = 0; t < l; t++)
        pows[1][t] = (t < xi->size) ? xi->d[t] : 0;
    radix_integer_clear(xi, radix);

    for (t = 0; t < l; t++)
        pows[0][t] = 0;
    pows[0][0] = 1;
    for (i = 2; i <= b; i++)
        radix_mulmid(pows[i], pows[i - 1], l, pows[1], l, 0, l, radix);
}

/* cc = sum_{i=1}^{hi} coef[i] * pows[i]  (mod B^l), with single-limb
   coefficients coef[1..hi], accumulated in plain limbs through a three-word
   delayed-carry accumulator acc (3*l words); one mpn reduction per limb. */
static void
_block_inner_sum(nn_ptr cc, nn_ptr acc, nn_ptr * pows, slong hi,
    const ulong * coef, slong l, nmod_t Bmod)
{
    slong i, t;
    ulong cy[3];

    for (t = 0; t < 3 * l; t++)
        acc[t] = 0;
    for (i = 1; i <= hi; i++)
    {
        ulong c = coef[i];
        nn_srcptr pi = pows[i];

        for (t = 0; t < l; t++)
        {
            ulong phi, plo;
            umul_ppmm(phi, plo, pi[t], c);
            add_sssaaaaaa(acc[3 * t + 2], acc[3 * t + 1], acc[3 * t],
                          acc[3 * t + 2], acc[3 * t + 1], acc[3 * t],
                          (ulong) 0, phi, plo);
        }
    }

    cy[0] = cy[1] = cy[2] = 0;
    if (Bmod.norm == 0)
        for (t = 0; t < l; t++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                          acc[3 * t + 2], acc[3 * t + 1], acc[3 * t]);
            cc[t] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, Bmod.n, Bmod.ninv);
        }
    else
        for (t = 0; t < l; t++)
        {
            add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                          acc[3 * t + 2], acc[3 * t + 1], acc[3 * t]);
            cc[t] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, Bmod.n, Bmod.ninv, Bmod.norm);
        }
}

/*
    z = sum_{i=1}^{n} y^i / i  modulo p^N, by rectangular splitting

        sum_{i=1}^{n} y^i/i
            = sum_{j} ( sum_{i=1}^{b} y^i/(i+jb) ) y^{jb},

    evaluated by Horner in y^b.  Inner sums are single-limb-coefficient
    combinations of the precomputed powers y^i, accumulated with a three-word
    delayed-carry accumulator.  There are three paths.

    When n <= LOG_LCMTAB_MAXN the common denominator L = lcm(1, ..., n) fits one
    radix limb, so the numerator sum_{i=1}^{n} (L/i) y^i can be evaluated
    directly from the table (coefficients L/i) and divided by L exactly once at
    the end -- no per-block inversion, and b is free of the single-limb block
    constraint.  If in addition the working precision N + v_p(L) fits one limb,
    the numerator is a single Horner evaluation (the fast path).  Otherwise it
    is the multi-limb rectangular evaluation with one final division.

    When n > LOG_LCMTAB_MAXN the lcm exceeds the limb radix, so the sum is
    carried scaled by p^k (k = floor(log_p n) >= v_p(i) for i <= n) with a block
    denominator f = (1+jb)...(hi+jb) that fits a limb: each block needs one
    Hensel inversion and a digit shift by k - v_p(f), and p^k is stripped at the
    end.  The block width is b = min(floor(sqrt(n)), largest m with n^m < B) and
    l covers N + k + g digits, where g = b/(p-1) bounds a block's k - v_p(f)
    loss.  Does not support aliasing of z with y.
*/
static void
_padic_radix_log_rectangular_series(radix_integer_t z, const radix_integer_t y,
    slong n, slong N, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    ulong B = LIMB_RADIX(radix);
    nmod_t Bmod = radix->B;
    slong e = radix->exp;
    slong b, bfit, k, g, lN, l, i, j, t, nsums;
    nn_ptr scratch, zarr, cc, tmp, acc, cf;
    nn_ptr * pows;
    radix_integer_struct zv;

    if (n <= 0)
    {
        radix_integer_zero(z, radix);
        return;
    }

    if (n == 1)
    {
        radix_integer_mod_digits(z, y, N, radix);
        return;
    }

    if (n <= LOG_LCMTAB_MAXN)
    {
        const ulong * row = log_lcm_data + log_lcm_off[n]; /* row[i] = lcm(1..n)/i */
        ulong L = row[1];                                  /* L = lcm(1..n)        */

        if (L < B)
        {
            ulong Lu;                                      /* L = p^w * Lu, Lu unit */
            slong w;

            if (p == 2)
            {
                w = flint_ctz(L);
                Lu = L >> w;
            }
            else
            {
                w = (slong) _radix_valuation_digits_1(L, radix);
                Lu = L * radix->bpow_oddinv[w].a;
            }

            /* Single-limb fast path: N + v_p(L) fits one limb. */
            if (N + w <= e)
            {
                ulong yl, num, num_u, res, mag0;

                mag0 = (y->size == 0) ? 0 : y->d[0];       /* |y| mod B */
                yl = (y->size >= 0) ? mag0 : (mag0 == 0 ? 0 : B - mag0);

                num = _nmod_poly_evaluate_nmod(row, n + 1, yl, Bmod);
                num_u = (p == 2) ? (num >> w) : (num * radix->bpow_oddinv[w].a);

                if (!radix_divmod_bn_1(&res, NULL, &num_u, 1, Lu, 1, radix))
                    flint_throw(FLINT_ERROR, "_padic_radix_log_rectangular_series: "
                        "denominator is not invertible after stripping\n");

                res = res % radix->bpow[N];                /* reduce to N digits */
                radix_integer_set_ui(z, res, radix);
                return;
            }

            /* Multi-limb path with the lcm as a single common denominator:
               evaluate the numerator sum_{i=1}^{n} (L/i) y^i mod B^l, then strip
               p^w and divide by Lu once. */
            b = (slong) n_sqrt((ulong) n);
            if (b < 1)
                b = 1;
            lN = (N + e - 1) / e;
            l = lN + (w + e - 1) / e;                      /* covers N + w digits */

            /* pows[0..b] (b+1), zarr, cc, tmp (3), acc (3) -> b+7 slots */
            scratch = flint_malloc(sizeof(ulong) * (size_t) (b + 7) * (size_t) l);
            pows = flint_malloc(sizeof(nn_ptr) * (b + 1));
            for (i = 0; i <= b; i++)
                pows[i] = scratch + (slong) i * l;
            zarr = scratch + (slong) (b + 1) * l;
            cc   = scratch + (slong) (b + 2) * l;
            tmp  = scratch + (slong) (b + 3) * l;
            acc  = scratch + (slong) (b + 4) * l;

            _build_pows(pows, b, y, l, radix);

            for (t = 0; t < l; t++)
                zarr[t] = 0;
            nsums = (n + b - 1) / b;

            for (j = nsums - 1; j >= 0; j--)
            {
                slong base = j * b;
                slong hi = FLINT_MIN(b, n - base);

                /* cc = sum_{i=1}^{hi} (L/(i+base)) y^i = sum row[base+i] y^i */
                _block_inner_sum(cc, acc, pows, hi, row + base, l, Bmod);

                /* zarr <- zarr * y^b + cc */
                radix_mulmid(tmp, zarr, l, pows[b], l, 0, l, radix);
                radix_add(zarr, tmp, l, cc, l, radix);
            }

            flint_free(pows);

            /* zarr = sum (L/i) y^i = p^w Lu * sum y^i/i; strip p^w, divide by Lu */
            if (w > 0)
                _arr_rshift_digits(zarr, l, w, radix);
            if (!radix_divmod_bn_1(cc, NULL, zarr, l, Lu, lN, radix))
                flint_throw(FLINT_ERROR, "_padic_radix_log_rectangular_series: "
                    "denominator is not invertible after stripping\n");

            zv.d = cc;
            zv.alloc = lN;
            zv.size = lN;
            while (zv.size > 0 && zv.d[zv.size - 1] == 0)
                zv.size--;
            radix_integer_mod_digits(z, &zv, N, radix);

            flint_free(scratch);
            return;
        }
    }

    /* General per-block loop: n > LOG_LCMTAB_MAXN (or lcm >= B). */
    {
        ulong prod = 1, hh, ll;
        bfit = 0;
        for (;;)
        {
            umul_ppmm(hh, ll, prod, (ulong) n);
            if (hh != 0 || ll >= B)
                break;
            prod = ll;
            bfit++;
        }
    }

    b = (slong) n_sqrt((ulong) n);
    if (b < 1)
        b = 1;
    if (b > bfit)
        b = bfit;

    k = (slong) n_flog((ulong) n, p);             /* >= v_p(i) for all i <= n */
    g = b / (slong) (p - 1) + 1;                  /* bound on a block's w - k  */
    lN = (N + e - 1) / e;
    l = lN + (k + g + e - 1) / e;                 /* limbs covering N + k + g  */

    /* pows[0..b] (b+1), zarr, cc, tmp (3), acc (3) -> b+7 slots; cf separate. */
    scratch = flint_malloc(sizeof(ulong) * (size_t) (b + 7) * (size_t) l);
    pows = flint_malloc(sizeof(nn_ptr) * (b + 1));
    cf = flint_malloc(sizeof(ulong) * (b + 1));
    for (i = 0; i <= b; i++)
        pows[i] = scratch + (slong) i * l;
    zarr = scratch + (slong) (b + 1) * l;
    cc   = scratch + (slong) (b + 2) * l;
    tmp  = scratch + (slong) (b + 3) * l;
    acc  = scratch + (slong) (b + 4) * l;

    _build_pows(pows, b, y, l, radix);

    for (t = 0; t < l; t++)
        zarr[t] = 0;
    nsums = (n + b - 1) / b;

    for (j = nsums - 1; j >= 0; j--)
    {
        slong base = j * b;
        slong hi = FLINT_MIN(b, n - base);
        slong w;
        ulong f, fu;

        /* f = (1+base)...(hi+base); cf[i] = f/(i+base), all single limbs */
        f = 1;
        for (i = 1; i <= hi; i++)
            f *= (ulong) (i + base);
        for (i = 1; i <= hi; i++)
            cf[i] = f / (ulong) (i + base);

        _block_inner_sum(cc, acc, pows, hi, cf, l, Bmod);

        /* f = p^w * fu (fu a unit) */
        w = 0;
        fu = f;
        while (fu % p == 0)
        {
            fu /= p;
            w++;
        }

        /* cc <- cc * p^{k-w} */
        if (w > k)
            _arr_rshift_digits(cc, l, w - k, radix);
        else if (w < k)
            _arr_lshift_digits(cc, l, k - w, radix);

        /* block contribution tmp = cc / fu  (mod B^l), fu a single-limb unit */
        if (!radix_divmod_bn_1(tmp, NULL, cc, l, fu, l, radix))
            flint_throw(FLINT_ERROR, "_padic_radix_log_rectangular_series: "
                "denominator is not invertible after stripping\n");

        /* zarr <- zarr * y^b + tmp */
        radix_mulmid(cc, pows[b], l, zarr, l, 0, l, radix);
        radix_add(zarr, cc, l, tmp, l, radix);
    }

    flint_free(pows);
    flint_free(cf);

    /* zarr holds p^k * sum y^i/i; strip the guard and reduce to p^N. */
    if (k > 0)
        _arr_rshift_digits(zarr, l, k, radix);

    zv.d = zarr;
    zv.alloc = l;
    zv.size = l;
    while (zv.size > 0 && zv.d[zv.size - 1] == 0)
        zv.size--;

    radix_integer_mod_digits(z, &zv, N, radix);

    flint_free(scratch);
}
/*
    rop = log(1 - y) modulo p^N by rectangular splitting, with y an exact
    integer of p-adic valuation at least the convergence threshold and < N.
    The nonnegative residue is returned.
*/
void
_padic_radix_log_rectangular(radix_integer_t rop, const radix_integer_t y,
    slong N, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    ulong B = LIMB_RADIX(radix);
    slong vr, n;
    radix_integer_t s, zero;

    if (N <= 0)
    {
        radix_integer_zero(rop, radix);
        return;
    }

    vr = radix_integer_valuation_digits(y, radix);
    if (vr >= N)
    {
        radix_integer_zero(rop, radix);
        return;
    }

    n = _padic_radix_log_bound(vr, N, p) - 1;

    /* Probably an artificially small limb radix, not supported by the implementation */
    if (n >= 2 && (ulong) n >= B)
    {
        _padic_radix_log_balanced(rop, y, N, radix);
        return;
    }

    radix_integer_init(s, radix);
    radix_integer_init(zero, radix);

    _padic_radix_log_rectangular_series(s, y, n, N, radix);

    /* log(1 - y) = - sum_{i>=1} y^i/i */
    radix_integer_zero(zero, radix);
    radix_integer_sub(rop, zero, s, radix);
    radix_integer_mod_digits(rop, rop, N, radix);

    radix_integer_clear(s, radix);
    radix_integer_clear(zero, radix);
}
