/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "padic_polyxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

padicxx_ctx ctx(fmpzxx(5), 10, 20, PADIC_TERSE);

void
test_init()
{
    padic_polyxx a(ctx, 20);
    tassert(a.get_ctx()._ctx() == ctx._ctx());
    tassert(a.prec() == 20);

    padic_polyxx c(a);
    tassert(a == c);

    padic_polyxx b(ctx, 30);
    tassert((a + b).estimate_ctx()._ctx() == ctx._ctx());
    tassert((a + b).prec() == 30);

    tassert((a + b).create_temporary().prec() == 30);
    padic_polyxx d((a + b).create_temporary());
    tassert(d.get_ctx()._ctx() == ctx._ctx());
    tassert(d.prec() == 30);

    padic_polyxx e(a + b);
    tassert(e.prec() == 30);

    tassert(padic_polyxx::zero(ctx).is_zero()
            && padic_polyxx::zero(ctx, 25).is_zero());
    tassert(padic_polyxx::one(ctx).is_one()
            && padic_polyxx::one(ctx, 25).is_one());

    a.prec() = 25;
    tassert(a.prec() == 25);
    a.val() = 17;
    tassert(a.val() == 17);
}

void
test_assignment()
{
    padic_polyxx a(ctx, 20), b(ctx, 20);
    ulong o1 = 1;
    slong o2 = 1;
    fmpzxx o3(1);
    fmpqxx o4 = fmpqxx::integer(o3);
    fmpz_polyxx o5 = fmpz_polyxx::from_ground(o3);
    fmpq_polyxx o6;o6.set_coeff(0, o4);

    a = o1;tassert(a.is_one());a = 0;
    a = o2;tassert(a.is_one());a = 0;
    a = o3;tassert(a.is_one());a = 0;
    a = o4;tassert(a.is_one());a = 0;
    a = o5;tassert(a.is_one());a = 0;
    a = o6;tassert(a.is_one());

    tassert(a != b);
    a = b;
    tassert(a == b);

    a = padicxx::from_QQ(1, ctx);tassert(a.is_one());
}

void
test_conversion()
{
    tassert(padic_polyxx::from_QQ(1, ctx).is_one());
    tassert(padic_polyxx::from_QQX(fmpz_polyxx::from_ground(1), ctx).is_one());
    tassert(padic_polyxx::from_ground(padicxx::one(ctx)).is_one());

    frandxx state;
    fmpz_polyxx f = fmpz_polyxx::randtest_unsigned(state, 5, 10);
    tassert(padic_polyxx::from_QQX(f, ctx, 100).to<fmpz_polyxx>() == f);
    fmpq_polyxx fq;fq = f;
    tassert(padic_polyxx::from_QQX(f, ctx, 100).to<fmpq_polyxx>() == fq);
}

void
test_arithmetic()
{
    frandxx state;
    padicxx a = padicxx::randtest(state, ctx);
    padicxx b = padicxx::randtest(state, ctx);
    padic_polyxx ap = padic_polyxx::from_ground(a);
    padic_polyxx bp = padic_polyxx::from_ground(b);

    tassert(ap + bp == padic_polyxx::from_ground(a + b));
    tassert(ap - bp == padic_polyxx::from_ground(a - b));
    tassert(ap * bp == padic_polyxx::from_ground(a * b));
    tassert(-ap == padic_polyxx::from_ground(-a));
    tassert(ap * b == padic_polyxx::from_ground(a * b));
}

void
test_functions()
{
    frandxx state;

    padic_polyxx p = padic_polyxx::randtest(state, 4, ctx, 25);
    tassert(p.pow(3u) == p*p*p);
    tassert(p.degree() == 3);
    tassert(p.length() == 4);
    tassert(p.val() + 1 == (p*padicxx::from_QQ(5, ctx, 25)).val());

    fmpz_polyxx F1 = fmpz_polyxx::randtest(state, 5, 10);
    fmpz_polyxx F2 = fmpz_polyxx::randtest(state, 5, 10);
    padic_polyxx f1 = padic_polyxx::from_QQX(F1, ctx);
    padic_polyxx f2 = padic_polyxx::from_QQX(F2, ctx);
    tassert(f1.derivative() == padic_polyxx::from_QQX(F1.derivative(), ctx));
    tassert(compose(f1, f2) == padic_polyxx::from_QQX(F1(F2), ctx));
    tassert(f1(f2) == compose(f1, f2));

    fmpzxx X = fmpzxx::randtest(state, 10);
    padicxx x = padicxx::from_QQ(X, ctx);
    tassert(f1(x) == padicxx::from_QQ(F1(X), ctx));

    tassert(f1.shift_left(5) == padic_polyxx::from_QQX(F1.shift_left(5), ctx));
    tassert(f1.shift_right(2)
            == padic_polyxx::from_QQX(F1.shift_right(2), ctx));

    F1.set_coeff(0, 1);
    f1 = padic_polyxx::from_QQX(F1, ctx);
    tassert(f1.inv_series(10) == padic_polyxx::from_QQX(F1.inv_series(10), ctx));

    f2.set_zero();
    f2.set_coeff(1, padicxx::one(ctx));
    tassert(p.compose_pow(5) == p(f2.pow(5u)));

    tassert(!padic_polyxx::randtest_not_zero(state, 5, ctx).is_zero());
    tassert(padic_polyxx::randtest_val(state, 5, 5, ctx).val() == 5);
}

// test stuff which we should get automatically - references etc
void
test_extras()
{
    padic_polyxx a = padic_polyxx::from_QQ(fmpqxx(3, 5u), ctx);
    padic_polyxx b = padic_polyxx::from_QQ(fmpqxx(3, 1u), ctx);

    padic_polyxx_ref ar(a);
    padic_polyxx_srcref asr(a);
    tassert(a == ar && ar == asr);
    ar = 3;
    tassert(a == b && asr == b);

    tassert(ar + asr == a + a);
    tassert(padic_polyxx(ar) == ar);
    tassert(padic_polyxx(asr) == ar);
}

void
test_prec()
{
    padic_polyxx a(ctx, 5), b(ctx, 7);
    padic_polyxx_ref ar(a);
    padic_polyxx_srcref br(b);

    tassert((a + a).prec() == 5);
    tassert((a + ar).prec() == 5);
    tassert((a + b).prec() == 7);
    tassert((a + br).prec() == 7);
    tassert((a.toN(15) + br.toN(10)).prec() == 15);
}

void
test_printing()
{
    padic_polyxx p(ctx);
    p.set_coeff(2, padicxx::one(ctx));
    tassert_fprint(p, "3  0 0 1");
    tassert_fprint_pretty(p, "x", "x^2");
}

int
main()
{
    std::cout << "padic_polyxx....";

    test_init();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_extras();
    test_prec();
    test_printing();

    std::cout << "PASS" << std::endl;
    return 0;
}
