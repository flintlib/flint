/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "qadicxx.h"
#include "flintxx/test/helpers.h"

// This is a test in itself: used to not compile.
#include "padic_polyxx.h"

using namespace flint;

void
test_init()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);

    qadicxx a(ctx, 20);
    tassert(&a.get_qctx() == &ctx);
    tassert(a.prec() == 20);

    qadicxx c(a);
    tassert(a == c);

    qadicxx b(ctx, 30);
    tassert(&(a + b).estimate_ctx() == &ctx);
    tassert((a + b).prec() == 30);

    tassert((a + b).create_temporary().prec() == 30);
    qadicxx d((a + b).create_temporary());
    tassert(&d.get_qctx() == &ctx);
    tassert(d.prec() == 30);

    qadicxx e(a + b);
    tassert(e.prec() == 30);

    tassert(qadicxx::zero(ctx).is_zero() && qadicxx::zero(ctx, 25).is_zero());
    tassert(qadicxx::one(ctx).is_one() && qadicxx::one(ctx, 25).is_one());

    qadicxx zero = qadicxx::zero(ctx);
    tassert((zero + zero).val() == 0);
}

void
test_assignment()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    qadicxx a(ctx, 20), b(ctx, 20);

    a = 17u; tassert(a != b);
    b = 17u; tassert(a == b);
}

void
test_conversion()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);

    frandxx state;
    padicxx p = padicxx::randtest(state, ctx.pctx());
    qadicxx a = qadicxx::from_ground(ctx, p);
    tassert(a.to<padicxx>() == p);
    assert_exception(qadicxx::gen(ctx).to<padicxx>());
}

void
test_arithmetic()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    frandxx state;
    padicxx pa = padicxx::randtest(state, ctx.pctx());
    padicxx pb = padicxx::randtest(state, ctx.pctx());
    qadicxx qa = qadicxx::from_ground(ctx, pa);
    qadicxx qb = qadicxx::from_ground(ctx, pb);

    tassert((qa + qb).to<padicxx>() == pa + pb);
    tassert((qa - qb).to<padicxx>() == pa - pb);
    tassert((qa * qb).to<padicxx>() == pa * pb);
    tassert((-qa).to<padicxx>() == -pa);
}

void
test_functions()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    frandxx state;
    padicxx pa = padicxx::randtest_not_zero(state, ctx.pctx());
    padicxx pb(padicxx::randtest_int(state, ctx.pctx())
                 * padicxx::from_QQ(5*5*5, ctx.pctx()));
    qadicxx qa = qadicxx::from_ground(ctx, pa);
    qadicxx qb = qadicxx::from_ground(ctx, pb);
    qadicxx qc = qadicxx::randtest_val(state, 5, ctx);

    tassert(qa.inv().to<padicxx>() == pa.inv());
    tassert(qc.pow(fmpzxx(3)) == qc*qc*qc);

    tassert(qb.exp().to<padicxx>() == pb.exp());
    tassert(qc.exp() == qc.exp_rectangular());
    tassert(qc.exp() == qc.exp_balanced());

    tassert((qadicxx::one(ctx) + qb).log().to<padicxx>()
            == (padicxx::one(ctx.pctx()) + pb).log());
    qc += qadicxx::one(ctx);
    tassert(qc.log() == qc.log_balanced());

    qadicxx res(ctx, PADIC_DEFAULT_PREC);
    qadic_frobenius(res._qadic(), qc._qadic(), 3, ctx._ctx());
    tassert(res == qc.frobenius(3));

    qadic_teichmuller(res._qadic(), qc._qadic(), ctx._ctx());
    tassert(res == qc.teichmuller());

    padicxx resp(ctx.pctx(), PADIC_DEFAULT_PREC);
    qadic_trace(resp._padic(), qc._qadic(), ctx._ctx());
    tassert(resp == qc.trace());

    qadic_norm(resp._padic(), qc._qadic(), ctx._ctx());
    tassert(resp == qc.norm());
    tassert(qc.norm() == qc.norm_analytic());
    tassert(qc.norm() == qc.norm_resultant());
}

// test stuff which we should get automatically - references etc
void
test_extras()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    qadicxx a(ctx);
    qadicxx b(ctx);b = 3u;

    qadicxx_ref ar(a);
    qadicxx_srcref asr(a);
    tassert(a == ar && ar == asr);
    ar = 3u;
    tassert(a == b && asr == b);

    tassert(ar + asr == a + a);
    tassert(qadicxx(ar) == ar);
    tassert(qadicxx(asr) == ar);
}

void
test_prec()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    qadicxx a(ctx, 5), b(ctx, 7);
    qadicxx_ref ar(a);
    qadicxx_srcref br(b);

    tassert((a + a).prec() == 5);
    tassert((a + ar).prec() == 5);
    tassert((a + b).prec() == 7);
    tassert((a + br).prec() == 7);
    tassert((a.toN(15) + br.toN(10)).prec() == 15);
}

void
test_printing()
{
    qadicxx_ctx ctx(fmpzxx(5), 3, 10, 20, PADIC_TERSE);
    if(0)
        print(ctx); // make sure this compiles
    tassert_fprint_pretty(qadicxx::gen(ctx).pow(fmpzxx(2)), "x^2");
}

int
main()
{
    std::cout << "qadicxx....";

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
