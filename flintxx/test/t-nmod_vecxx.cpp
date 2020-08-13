/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "nmod_vecxx.h"
#include "flintxx/ltuple.h"
#include "flintxx/test/helpers.h"

using namespace flint;

bool test_conv(nmodxx_ctx_srcref) {return true;}
void
test_init()
{
    nmodxx_ctx c(31);
    tassert(test_conv(c));
    tassert(c.n() == 31);

    nmodxx a = nmodxx::make_nored(65, c);
    tassert(a._limb() == 65);
    a.reduce();
    tassert(a._limb() == (65 % 31));

    nmodxx b = nmodxx::red(65, c);
    tassert(b._limb() == a._limb());

    tassert(a == b);
    nmodxx zero(c);
    tassert(zero != a);

    a = zero;
    tassert(a == zero);

    tassert(nmodxx::red(65, c) == nmodxx::red(fmpzxx(65), c));
    tassert(nmodxx::red(3, c)/nmodxx::red(2, c)
            == nmodxx::red(fmpqxx(3, 2u), c));
}

void
test_conversion()
{
    nmodxx_ctx c(31);
    nmodxx a = nmodxx::red(65, c);
    tassert(a.to<mp_limb_t>() == a._limb());
    tassert(a.to_string() == "3 mod 31");
}

void
test_arithmetic()
{
    nmodxx_ctx c(31);
    nmodxx a = nmodxx::red(65, c);
    nmodxx b = nmodxx::red(21, c);

    tassert(a + b == nmodxx::red(65 + 21, c));
    tassert((a + b).estimate_ctx() == a._ctx());
    tassert((-a) + a == nmodxx::red(0, c));
    tassert(b + (-a) == b - a);
    tassert(nmodxx::red(3, c) * a == a + a + a);
    tassert((a / b) * b == a);
    tassert(inv(b) == nmodxx::red(1, c) / b);
    tassert(pow(a, 4u) == a*a*a*a);

    tassert(b.inv() == inv(b));
    tassert(a.pow(4u) == pow(a, 4u));
}

void
test_references()
{
    nmodxx_ctx c(31);
    nmodxx a = nmodxx::red(17, c);
    nmodxx b = nmodxx::red(19, c);
    nmodxx_ref ar(a);
    nmodxx_srcref bsr(b);
    tassert(a == ar);
    tassert(bsr != ar);
    ar = nmodxx::red(19, c);
    tassert(a == bsr);
    tassert(ar + bsr + a == nmodxx::red(3, c) * b);

    nmodxx_ref ar2 = nmodxx_ref::make(a._limb(), a.estimate_ctx());
    nmodxx_srcref asr2 = nmodxx_srcref::make(a._limb(), a.estimate_ctx());
    ar2 = nmodxx::red(4, c);
    tassert(asr2 == nmodxx::red(4, c));
}

void
test_vec()
{
    nmodxx_ctx c(31);
    nmod_vecxx v1(10, c), v2(10, c);
    for(slong i = 0;i < (slong) v1.size();++i)
    {
        v1[i] = nmodxx::red(0, c);
        v2[i] = nmodxx::red(0, c);
    }
    tassert(v1 == v2);
    v1[0] = nmodxx::red(15, c);
    tassert(v1 != v2);
    tassert(v1 + v2 == v1);
}

namespace flint {
typedef make_ltuple<mp::make_tuple<nmodxx>::type>::type ltup_t;
FLINT_DEFINE_UNOP(make_lazy_test)
namespace rules {
FLINT_DEFINE_UNARY_EXPR_COND(make_lazy_test_op, ltup_t, NMOD_VECXX_COND_S,
        to.template get<0>() = from[0])
}
}
void
test_temporaries()
{
    nmodxx_ctx ctx(17);
    nmod_vecxx v(1, ctx);v[0] = nmodxx::red(7, ctx);
    tassert(make_lazy_test(v).get<0>() == v[0]);
}

int
main()
{
    std::cout << "nmod_vecxx....";

    test_init();
    test_conversion();
    test_arithmetic();
    test_references();
    test_vec();
    test_temporaries();

    std::cout << "PASS" << std::endl;
    return 0;
}
