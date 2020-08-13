/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "padic_matxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

padicxx_ctx ctx(fmpzxx(5), 10, 20, PADIC_TERSE);

void
test_init()
{
    padic_matxx a(ctx, 5, 7, 20);
    tassert(a.get_ctx()._ctx() == &ctx._ctx()[0]);
    tassert(a.prec() == 20);
    tassert(a.rows() == 5 && a.cols() == 7);

    padic_matxx c(a);
    tassert(a == c);

    padic_matxx b(ctx, 3, 4, 30);
    tassert((a + b).estimate_ctx()._ctx() == &ctx._ctx()[0]);
    tassert((a + b).prec() == 30);

    tassert((a + b).create_temporary().prec() == 30);
    padic_matxx d((a + b).create_temporary());
    tassert(d.get_ctx()._ctx() == &ctx._ctx()[0]);
    tassert(d.prec() == 30);

    padic_matxx e(a + b);
    tassert(e.prec() == 30);

    tassert(padic_matxx::zero(ctx, 5, 7).is_zero()
            && padic_matxx::zero(ctx, 5, 7, 25).is_zero());

    a.prec() = 25;
    tassert(a.prec() == 25);
    a.val() = 17;
    tassert(a.val() == 17);

    padic_matxx eye = padic_matxx::one(ctx, 3, 3);
    for(slong i = 0;i < eye.rows();++i)
        for(slong j = 0;j < eye.cols();++j)
            if(i == j)
                tassert(eye.get_entry(i, j).is_one());
            else
                tassert(eye.get_entry(i, j).is_zero());
}

void
test_manipulation()
{
    padic_matxx a(ctx, 5, 7);
    tassert(a.at(0, 0) == 0);
    padic_matxx b(a);
    b.set_entry(0, 0, padicxx::one(ctx));
    tassert(a != b);
    b /= fmpzxx(5);
    tassert(b.at(0, 0) == 1);
    b.at(0, 0) = 2;
    tassert((b / fmpzxx(5)).val() == -2);
    tassert(b.get_entry(0, 0) == padicxx::from_QQ(fmpqxx::frac(2, 5), ctx));
    tassert(b.get_entry(0, 0)
            == (b + padic_matxx::zero(ctx, 5, 7)).get_entry(0, 0));
}

void
test_assignment()
{
    padic_matxx a(ctx, 3, 7), b(ctx, 3, 7);
    b.set_entry(0, 0, padicxx::from_QQ(3, ctx));
    tassert(a != b);

    fmpq_matxx M(3, 7);
    b = M;
    tassert(a == b);
    tassert(b.to<fmpq_matxx>() == M);
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
void
test_arithmetic()
{
    frandxx state;

    padic_matxx v(ctx, 1, 10);
    padic_matxx A(ctx, 10, 10);

    tassert(!has_explicit_temporaries(transpose(transpose(v))));
    tassert(!has_explicit_temporaries(transpose(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(A + v*transpose(v)));
    tassert(!has_explicit_temporaries(transpose((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(transpose(v*transpose(v) + v*transpose(v))));
    tassert(!has_explicit_temporaries(v*transpose(v) + v*transpose(v)));

    fmpz_matxx a = fmpz_matxx::randtest(4, 4, state, 10);
    fmpz_matxx b = fmpz_matxx::randtest(4, 4, state, 10);
    fmpq_matxx aq(4, 4);aq = a;
    fmpq_matxx bq(4, 4);bq = b;
    padic_matxx ap = padic_matxx::from_QQ(aq, ctx);
    padic_matxx bp = padic_matxx::from_QQ(bq, ctx);

    tassert(ap + bp == padic_matxx::from_QQ(aq + bq, ctx));
    tassert(ap - bp == padic_matxx::from_QQ(aq - bq, ctx));
    tassert(ap * bp == padic_matxx::from_QQ(aq * bq, ctx));
    tassert(-ap == padic_matxx::from_QQ(-aq, ctx));

    fmpzxx s = fmpzxx::randtest(state, 15);
    tassert(ap * s == padic_matxx::from_QQ(aq * s, ctx));
    ltupleref(_, s) = s.remove(ctx.get_p());
    tassert(ap / s == padic_matxx::from_QQ(aq / s, ctx));
}

void
test_functions()
{
    frandxx state;
    padic_matxx A = padic_matxx::randtest(3, 4, state, ctx);
    padic_matxx B(A.transpose());
    for(slong i = 0;i < B.rows();++i)
        for(slong j = 0;j < B.cols();++j)
            tassert(B.get_entry(i, j) == A.get_entry(j, i));

    tassert(!A.is_square() && !A.is_empty());
}

// test stuff which we should get automatically - references etc
void
test_extras()
{
    padic_matxx a(ctx, 2, 2);
    padic_matxx b(ctx, 2, 2);
    a.set_entry(0, 0, padicxx::from_QQ(3, ctx));
    b.set_entry(0, 1, padicxx::from_QQ(2, ctx));

    padic_matxx_ref ar(a);
    padic_matxx_srcref asr(a);
    tassert(a == ar && ar == asr);

    tassert(ar + asr == a + a);
    tassert(padic_matxx(ar) == ar);
    tassert(padic_matxx(asr) == ar);

    ar.set_entry(0, 0, padicxx::from_QQ(1, ctx));
    tassert(a.get_entry(0, 0).is_one());
}

void
test_prec()
{
    padic_matxx a(ctx, 2, 2, 5), b(ctx, 2, 2, 7);
    padic_matxx_ref ar(a);
    padic_matxx_srcref br(b);

    tassert((a + a).prec() == 5);
    tassert((a + ar).prec() == 5);
    tassert((a + b).prec() == 7);
    tassert((a + br).prec() == 7);
    tassert((a.toN(15) + br.toN(10)).prec() == 15);
}

void
test_printing()
{
    tassert_fprint(padic_matxx::one(ctx, 2, 2), "2 2  1 0 0 1");
    tassert_fprint_pretty(padic_matxx::one(ctx, 2, 2), "[[1 0]\n[0 1]]");
}

int
main()
{
    std::cout << "padic_matxx....";

    test_init();
    test_manipulation();
    test_assignment();
    test_arithmetic();
    test_functions();
    test_extras();
    test_prec();
    test_printing();

    std::cout << "PASS" << std::endl;
    return 0;
}
