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

#include "fmpz_mod_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init(fmpz_modxx_ctx& M)
{
    M.set_modulus(2003);
    fmpz_mod_polyxx p(M);
    tassert(p.length() == 0);
    tassert(p.modulus() == 2003);
    tassert(fmpz_mod_polyxx::zero(M).is_zero());
}

void
test_manipulation(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx p(M), q(M);
    p.set_coeff(5, 17u + M.modulus());
    tassert(p.degree() == 5);
    q.set_coeff(5, fmpzxx(16) + fmpzxx(1));
    tassert((q + fmpz_mod_polyxx(M)).get_coeff(5) == 17);
    p.set_coeff(0, 1);
    tassert(p != q);
    p.set_coeff(0, 0);
    tassert(p == q);

    tassert((p + p).lead() == 2*p.lead());

    q.lead() = 0;
    q._normalise();
    tassert(q.is_zero());

    p.realloc(0);
    tassert(p.is_zero());
}

void
test_assignment(fmpz_modxx_ctx& M)
{
    M.set_modulus(31);
    fmpz_mod_polyxx p(M), q(M);
    p.set_coeff(0, 1);
    tassert(p != q);
    p = q;
    tassert(p == q);
}

void
test_conversion(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx p(M);

    p = 4u + 1031;
    tassert(p.length() == 1 && p.get_coeff(0) == 4);

    p = fmpzxx(5) + M.modulus();
    tassert(p.length() == 1 && p.get_coeff(0) == 5);

    frandxx rand;
    fmpz_polyxx P = fmpz_polyxx::randtest(rand, 10, 20);
    p = P;
    for(slong i = 0;i < P.length();++i)
        tassert(P.get_coeff(i) % M.modulus() == p.get_coeff(i));
    fmpz_polyxx Pp = p.to<fmpz_polyxx>();
    for(slong i = 0;i < P.length();++i)
        tassert(P.get_coeff(i) % M.modulus() == Pp.get_coeff(i));
}

void
test_arithmetic(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx g(M), h(M);
    g.set_coeff(0, 17); h.set_coeff(0, 15u + M.modulus());
    tassert((g + h).get_coeff(0) == 15 + 17);

    frandxx state;
    g.set_randtest(state, 10);
    h.set_randtest(state, 10);

    tassert(((-g) + g).is_zero());
    tassert(g - h == g + (-h));

    tassert(g*fmpzxx(3) == g + g + g);
    tassert(g.make_monic() == g*g.lead().invmod(M.modulus()));

    fmpz_mod_polyxx f(M);f = 15u;
    tassert(f*g == fmpzxx(15)*g);

    f = h*g;f.truncate(7);
    tassert(f == mullow(h, g, 7));

    f = h / g;
    tassert(f*g + (h % g) == h);
    tassert(((h*g) % h).is_zero());

    f.set_randtest(state, 10);
    tassert(h.mulmod(g, f) == ((h*g) % f));

    fmpz_mod_polyxx X(M);X.set_coeff(1, 1);
    fmpz_mod_polyxx one(M);one.set_coeff(0, 1);
    f = X*X + one;
    fmpzxx x(7);
    tassert(evaluate(f, x) == x*x + 1u);
    tassert(f(x) == evaluate(f, x));

    fmpz_mod_polyxx seven(M);
    seven.set_coeff(0, x);
    tassert(compose(f, seven).get_coeff(0) == f(x));
    tassert(f(seven).length() == 1);
}

void
test_functions(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx g(M), res(M);

    g.set_coeff(5, 15);

    g.truncate(3);
    tassert(g.is_zero());

    g.set_coeff(15, 1);
    fmpz_mod_polyxx one(M);one = 1u;
    tassert(g.shift_right(15) == one);
    tassert(g.shift_right(15).shift_left(15) == g);

    frandxx rand;
    g.set_randtest(rand, 15);
    tassert(g.length() <= 15);
    g.set_randtest_irreducible(rand, 15);
    tassert(g.length() <= 15);
    g.set_randtest_not_zero(rand, 15);
    tassert(g.length() <= 15 && !g.is_zero());

    g.set_coeff(15, 1);
    g.zero_coeffs(14, 15);
    tassert(g.get_coeff(14) == 0);

    // multiplication, division, modulo tested in arithmetic

    tassert(g.pow(3u) == g*g*g);

    res = g.pow(15u);res.truncate(12);
    tassert(res == g.pow_trunc(15u, 12));
    tassert(res == g.pow_trunc_binexp(15u, 12));

    fmpz_mod_polyxx f(M);f.set_randtest(rand, 10);
    res = g.pow(10u) % f;
    tassert(res == g.powmod_binexp(10u, f));
    tassert(res == g.powmod_binexp(fmpzxx(10), f));

    fmpz_mod_polyxx tmp(M);
    ltupleref(res, tmp) = f.gcdinv(g);
    tassert(res == gcd(f, g) && tmp*f % g == res);

    g.set_randtest_irreducible(rand, 5);
    tassert(f.invmod(g)*f % g == one);
    assert_exception((f*g).invmod(g).evaluate());

    res = g*f;
    res.remove(f);
    tassert(res == g);

    fmpz_polyxx lift;
    lift = "5  1 1 1 1 1";
    res = lift;
    tassert(res.derivative().to<fmpz_polyxx>().to_string() == "4  1 2 3 4");

    tassert(f.divrem(g) == ltuple(f / g, f % g));
    tassert(f.divrem_basecase(g) == f.divrem(g));
    tassert(f.divrem_divconquer(g) == f.divrem(g));
    tassert(f.divrem_f(g) == ltuple(1, f / g, f % g));

    tassert(f.div_basecase(g) == f / g);
    tassert(f.rem_basecase(g) == f % g);

    f.set_coeff(0, 17);
    res = f*f.inv_series_newton(15);res.truncate(15);
    tassert(res == one);

    tassert(f(g) == f.compose_divconquer(g));
    tassert(f(g) == f.compose_horner(g));

    fmpz_mod_polyxx h(M);
    h.set_randtest(rand, 15);
    tassert(f.compose_mod(g, h) == f(g) % h);
    tassert(f.compose_mod(g, h) == f.compose_mod_horner(g, h));
    tassert(f.compose_mod(g, h) == f.compose_mod_brent_kung(g, h));

    h.set_randtest_irreducible(rand, 12);
    tassert(h.gcd(f) == one);
    tassert(f.gcd_euclidean(f) == f.make_monic());
    tassert(f.gcd_f(g) == ltuple(1, f.gcd(g)));
    tassert(f.gcd_euclidean_f(g) == ltuple(1, f.gcd(g)));

    fmpz_mod_polyxx R(M), S(M);
    ltupleref(res, R, S) = f.xgcd(g);
    tassert(res == R*f + S*g && res == gcd(f, g));
    tassert(f.xgcd(g) == f.xgcd_euclidean(g));
}

bool equiv_fac(const fmpz_mod_poly_factorxx& fac1,
        const fmpz_mod_poly_factorxx& fac2)
{
    tassert(fac1.size() == 2);
    if(fac1.exp(0) == fac1.exp(1))
    {
        if(fac2.exp(0) != fac1.exp(0) || fac2.exp(1) != fac1.exp(0))
            return false;
        return (fac1.p(0) == fac2.p(0) && fac1.p(1) == fac2.p(1))
            || (fac1.p(1) == fac2.p(0) && fac1.p(0) == fac2.p(1));
    }
    if(fac1.size() != fac2.size())
        return false;
    if(fac1.exp(0) == fac2.exp(0))
        return fac1.exp(1) == fac2.exp(1)
            && fac1.p(0) == fac2.p(0)
            && fac1.p(1) == fac2.p(1);
    else
        return fac1.exp(0) == fac2.exp(1)
            && fac1.exp(1) == fac2.exp(0)
            && fac1.p(0) == fac2.p(1)
            && fac1.p(1) == fac2.p(0);
}
void
test_factoring(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx f(M), g(M);
    frandxx state;
    f.set_randtest_irreducible(state, 4); f = f.make_monic();
    g.set_randtest_irreducible(state, 5); g = g.make_monic();

    fmpz_mod_poly_factorxx fac(M);
    fac = factor(f*f*g);
    tassert(fac.size() == 2);
    if(fac.exp(0) == 1)
    {
        tassert(fac.p(0) == g);
        tassert(fac.p(1) == f && fac.exp(1) == 2);
    }
    else
    {
        tassert(fac.p(0) == f && fac.exp(0) == 2);
        tassert(fac.p(1) == g && fac.exp(1) == 1);
    }

    fmpz_mod_poly_factorxx fac2(M);fac2 = fac;fac2.pow(2);
    fac.insert(g, 1);
    fac.insert(f, 2);
    tassert(fac == fac2);

    fmpz_mod_polyxx prod(f*f*f*g*g);
    fac = factor(prod);
    tassert(equiv_fac(fac, factor_cantor_zassenhaus(prod)));
    tassert(equiv_fac(factor(f*g), factor_berlekamp(f*g)));
    tassert(equiv_fac(fac, factor_kaltofen_shoup(prod)));

    std::vector<slong> degs(2);
    fac.realloc(0);fac.set_factor_distinct_deg(f*g, degs);
    tassert(degs.size() == 2);
    tassert((degs[0] == f.degree() && degs[1] == g.degree())
            || (degs[1] == f.degree() && degs[0] == g.degree()));

    tassert(f.is_irreducible() && f.is_irreducible_ddf()
            && f.is_irreducible_rabin());
    tassert(f.is_squarefree());
    // TODO test set_factor_equal_deg*

    if(0)
        print(fac); // test this compiles
}

void
test_randomisation(fmpz_modxx_ctx& M)
{
    frandxx state, state2;
    M.set_modulus(1031);
    fmpz_mod_polyxx p(M);

    p.set_randtest(state, 10);
    tassert(p == fmpz_mod_polyxx::randtest(M, state2, 10));
    p.set_randtest_irreducible(state, 10);
    tassert(p == fmpz_mod_polyxx::randtest_irreducible(M, state2, 10));
    p.set_randtest_not_zero(state, 10);
    tassert(p == fmpz_mod_polyxx::randtest_not_zero(M, state2, 10));
}

void
test_radix(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_poly_vecxx v1(10, M), v2(10, M);
    v1[0].set_coeff(7, 1);
    tassert(v1 != v2);
    fmpz_mod_poly_vecxx v3(v1);
    tassert(v3 == v1);
    v3[0].set_coeff(1, 1);
    tassert(v3 != v1);
    v2[0].set_coeff(7, 1);
    tassert(v1 == v2);

    frandxx rand;
    fmpz_mod_polyxx F = fmpz_mod_polyxx::randtest(M, rand, 10);
    fmpz_mod_polyxx R = fmpz_mod_polyxx::randtest(M, rand, 3);
    fmpz_mod_poly_vecxx b(F.degree() / R.degree() + 1, M);
    fmpz_mod_poly_radixxx rad(R, 15);
    b = F.radix(rad);
    tassert(b == F.radix(rad));

    fmpz_mod_polyxx f(M);
    for(slong i = 0;i < b.size();++i)
        f += b[i]*R.pow(static_cast<ulong>(i));
    tassert(f == F);
}

void
test_printing(fmpz_modxx_ctx& M)
{
    frandxx state;
    M.set_modulus(7);
    fmpz_mod_polyxx f = fmpz_mod_polyxx::randtest(M, state, 4);
    test_print_read(f);
    f.set_zero();
    f.set_coeff(0, 3);
    f.set_coeff(1, 1);
    tassert_fprint_pretty(f, "x", "x+3");
}

void
test_unified_access(fmpz_modxx_ctx& M)
{
    M.set_modulus(1031);
    fmpz_mod_polyxx p(M);
    p.set_coeff(0, 1);
    const fmpz_mod_polyxx& q = p;
    tassert(q.lead() == 1);
}

int
main()
{
    std::cout << "fmpz_mod_polyxx....";
    fmpz_modxx_ctx M(2);

    test_init(M);
    test_manipulation(M);
    test_assignment(M);
    test_conversion(M);
    test_arithmetic(M);
    test_functions(M);
    test_factoring(M);
    test_randomisation(M);
    test_radix(M);
    test_printing(M);
    test_unified_access(M);

    std::cout << "PASS" << std::endl;
    return 0;
}

