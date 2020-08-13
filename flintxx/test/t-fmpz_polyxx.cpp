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

#include "fmpz_polyxx.h"
#include "nmod_polyxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpz_polyxx p(10);
    tassert(p.length() == 0);

    tassert(fmpz_polyxx::zero().is_zero());
    tassert(fmpz_polyxx::one().is_one());
}

void
test_manipulation()
{
    fmpz_polyxx p, q;
    tassert(p.get_coeff(5) == 0);
    p.set_coeff(5, 17);
    tassert(p.degree() == 5);
    q.set_coeff(5, 17u);
    p.set_coeff(0, fmpzxx(1));
    q.coeff(0) = 1;
    tassert(p == q);

    tassert(p.lead() == 17);
    tassert(p.length() == 6);

    p.zero_coeffs(0, 6);
    tassert(p.is_zero() && !p.is_one() && !p.is_unit());
    p.set_coeff(0, 1);
    tassert(p.is_one() && p.is_unit());
    p.set_coeff(0, -1);
    tassert(p.is_unit());

    const fmpz_polyxx pc = p;
    tassert(p.get_coeff(0) == pc.get_coeff(0));
    tassert((p*q).coeff(0) == p.get_coeff(0)*q.get_coeff(0));
    tassert((p*q).lead() == p.lead()*q.lead());
}

void
test_assignment()
{
    fmpz_polyxx p, q;
    p = 1;
    tassert(p.is_one());
    q = UWORD(0);
    tassert(q.is_zero());
    tassert(p != q);
    p = q;
    tassert(p == q);
    p = "4  0 0 0 1";
    q.set_coeff(3, 1);
    tassert(p == q);

    // TODO XXX this does not always fail?
    //assert_exception(p = "4  1 2");
    assert_exception(p = "2 1 2");
}

void
test_conversion()
{
    fmpz_polyxx p;
    p.set_coeff(3, 1);
    tassert(p.to_string() == "4  0 0 0 1");
    tassert(p.pretty("x") == "x^3");

    tassert(p == fmpz_polyxx("4  0 0 0 1"));
}

void
test_arithmetic()
{
    fmpz_polyxx p, q;
    p = 1;
    q = "4  0 0 0 1";
    tassert((p + q).to_string() == "4  1 0 0 1");
    tassert((p - q).to_string() == "4  1 0 0 -1");
    tassert((-p).to_string() == "1  -1");

    fmpzxx two(2);
    tassert(two * q == 2 * q && 2u * q == q * 2
            && (two * q).to_string() == "4  0 0 0 2");
    q *= 2;
    tassert(q / two == q / 2u && q / 2 == q / two &&
            (q / two).to_string() == "4  0 0 0 1");
    // q == "4  0 0 0 2"

    q.set_coeff(1, 17); // q == "4 0 17 0 2"
    tassert((q % fmpzxx(5)).to_string() == "4  0 2 0 2");

    p = "3  1 0 1";
    tassert((p*q).to_string() == "6  0 17 0 19 0 2");

    tassert((p*q) / p == q);
    tassert(p + q % q == p);

    tassert(p(fmpzxx(1) + fmpzxx(1)) == 5);
    q = "3  0 0 1";
    tassert(p(q).to_string() == "5  1 0 0 0 1");
}

// Won't compile if the expression is not done using addmul
template<class T>
bool is_ternary(const T&)
{
    return T::ev_traits_t::temp_rule_t::TERNARY_OP_MARKER + 1;
}

// test stuff which we should get automatically - addmul, references etc
void
test_extras()
{
    // TODO addmul when we have it
}

ulong pow(ulong base, ulong exp)
{
    ulong res = 1;
    while(exp-- > 0)
        res *= base;
    return base;
}

void
test_functions()
{
    // test swap
    fmpz_polyxx p, q;
    p = 1; q = 0;
    swap(p, q);
    tassert(p.is_zero() && q.is_one());
    // p = 0, q = 1

    fmpz_polyxx f, g, xp1;
    f = "4  2 0 0 1"; // f = x^3 + 2
    g = "5  1 2 3 4 5";
    xp1 = "2  1 1";
    fmpzxx two(2);

    // test lazy functions
    tassert(reverse(q, 4).pretty("x") == "x^3");
    tassert(mul_2exp(f, 3u) == f * 8);
    tassert(f == fdiv_2exp(mul_2exp(f, 3u), 3u));
    tassert(tdiv(-f, two) == tdiv(-f, 2) && tdiv(-f, 2) == tdiv(-f, 2u)
            && tdiv(-f, two).to_string() == "4  -1 0 0 0");
    tassert(f == divexact(2*f, two) && f == divexact(2*f, 2u)
            && f == divexact(2*f, 2));
    tassert(-f == tdiv_2exp(-8*f - q, 3u));
    tassert(smod(5*f, fmpzxx(3)).to_string() == "4  1 0 0 -1");
    tassert(f == fmpz_polyxx::bit_unpack(bit_pack(f, 10u), 10u));
    tassert(f == fmpz_polyxx::bit_unpack_unsigned(bit_pack(f, 10u), 10u));

    tassert(mul_classical(f, g) == f*g);
    tassert(mul_karatsuba(f, g) == f*g);
    tassert(mul_SS(f, g) == f*g);
    tassert(mul_KS(f, g) == f*g);

    fmpz_polyxx res; res = f*g; res.truncate(3);
    tassert(mullow_classical(f, g, 3) == res);
    tassert(mullow_karatsuba_n(f, g, 3) == res);
    tassert(mullow_KS(f, g, 3) == res);
    tassert(mullow_SS(f, g, 3) == res);
    tassert(mullow(f, g, 3) == res);

    res = f*g; res.zero_coeffs(0, 5);
    tassert(mulhigh_classical(f, g, 5) == res);
    tassert(mulhigh_karatsuba_n(f, g, 6) == res);
    res = mulhigh_n(f, g, 5);
    res.zero_coeffs(0, 5);
    tassert(res == mulhigh_classical(f, g, 5));

    tassert(mulmid_classical(g, f).to_string() == "2  9 12");

    tassert(sqr(f) == f*f);
    tassert(sqr_KS(f) == f*f);
    tassert(sqr_karatsuba(f) == f*f);
    tassert(sqr_classical(f) == f*f);

    res = sqr(f);res.truncate(4);
    tassert(sqrlow(f, 4) == res);
    tassert(sqrlow_KS(f, 4) == res);
    tassert(sqrlow_karatsuba_n(f, 4) == res);
    tassert(sqrlow_classical(f, 4) == res);

    tassert(pow(f, 3u) == f*f*f);
    tassert(pow_multinomial(f, 3u) == f*f*f);
    tassert(pow_binexp(f, 3u) == f*f*f);
    tassert(pow_addchains(f, 3u) == f*f*f);
    res = pow(f, 10u);
    res.truncate(10);
    tassert(pow_trunc(f, 10u, 10) == res);
    fmpz_polyxx binomial;
    binomial = "2  1 2";
    tassert(pow_binomial(binomial, 3u) == binomial*binomial*binomial);

    fmpz_polyxx x; x = "2  0 1";
    tassert(shift_left(f, 5) == f*pow(x, 5u));
    tassert(shift_right(shift_left(f, 5), 5) == f);

    tassert(height(g) == 5);
    tassert(twonorm(g) == 7);

    tassert(resultant(f, g) == 1797);
    tassert(gcd(f, g).is_one());
    tassert(gcd_subresultant(f, g).is_one());
    tassert(gcd_heuristic(f, g).is_one());
    tassert(gcd_modular(f, g).is_one());

    fmpz_polyxx r, s;fmpzxx number;
    res = 1797;
    ltupleref(number, r, s) = xgcd(f, g);
    tassert(r*f + s*g == res && number == 1797);
    r = 0; s = 0; number = 0;
    ltupleref(number, r, s) = xgcd_modular(f, g);
    tassert(r*f + s*g == res && number == 1797);

    tassert(lcm(f, g) == f*g);

    tassert(content(2*g) == 2);
    tassert(primitive_part(2*g) == g);

    tassert(div_basecase(f*g, g) == f);
    tassert(div_divconquer(f*g, g) == f);
    tassert(rem_basecase(f + g, g) == f);
    res = 1;
    tassert(div_root(f*(x - res), fmpzxx(1)) == f);

    tassert(inv_series(xp1, 5).to_string() == "5  1 -1 1 -1 1");
    tassert(inv_series(xp1, 10) == inv_series_newton(xp1, 10));

    tassert(derivative(xp1).is_one());

    tassert((compose(xp1, f) - f).is_one());
    tassert((compose_divconquer(xp1, f) - f).is_one());
    tassert((compose_horner(xp1, f) - f).is_one());
    tassert(evaluate(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_horner(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_divconquer(xp1, fmpzxx(1)) == 2);
    tassert(evaluate_mod(xp1, 1u, 10u) == 2);

    fmpz_vecxx xs(3), ys(3);
    xs[0] = 0; xs[1] = 1; xs[2] = 2;
    ys[0] = 1; ys[1] = 2; ys[2] = 3;
    tassert(evaluate(xp1, xs) == ys);

    tassert(f(xp1) == taylor_shift(f, fmpzxx(1)));
    tassert(f(xp1) == taylor_shift_horner(f, fmpzxx(1)));
    tassert(f(xp1) == taylor_shift_divconquer(f, fmpzxx(1)));

    fmpz_polyxx inner;inner = "3  0 -1 5";
    res = g(inner); res.truncate(10);
    tassert(compose_series(g, inner, 10) == res);
    tassert(compose_series_horner(g, inner, 10) == res);
    tassert(compose_series_brent_kung(g, inner, 10) == res);

    res = "2  0 1";
    tassert(compose_series(inner, revert_series(inner, 10), 10) == res);
    tassert(compose_series(inner, revert_series_lagrange(inner, 10), 10)
            == res);
    tassert(compose_series(inner, revert_series_lagrange_fast(inner, 10), 10)
            == res);
    tassert(compose_series(inner, revert_series_newton(inner, 10), 10) == res);

    tassert(sqrt(f*f) == f);
    tassert(sqrt_classical(f*f) == f);
    assert_exception(sqrt(f*f + xp1).evaluate());
    assert_exception(sqrt_classical(f*f + xp1).evaluate());

    res = "5  32 0 0 0 1";
    tassert(bound_roots(res) >= 2);

    // test immediate functions
    p.set_coeff(3, 1);
    p.truncate(2);
    tassert(p.is_zero());
    tassert(f.max_limbs() == 1);
    tassert(g.max_bits() == 3);

    r = 0; s = 0;
    ltupleref(r, s) = divrem(g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");
    r = 0; s = 0;
    ltupleref(r, s) = divrem_basecase(g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");
    r = 0; s = 0;
    ltupleref(r, s) = divrem_divconquer(g, f);
    tassert(r*f + s == g);
    tassert(r.to_string() == "2  4 5");

    bool does_divide;
    ltupleref(does_divide, res) = divides(f*g, g);
    tassert(does_divide && res == f);

    tassert(div_series(f, xp1, 10) * xp1 % pow(x, 10u) == f);

    f *= 2;
    ulong d = 0;
    ltupleref(r, s, d) = pseudo_divrem(g, f);
    tassert(r*f + s == g*fmpzxx(pow(2, d) + 0.5));
    r = 0; r = 0; d = 0;
    ltupleref(r, s, d) = pseudo_divrem_basecase(g, f);
    tassert(r*f + s == g*fmpzxx(pow(2, d) + 0.5));
    r = 0; r = 0; d = 0;
    ltupleref(r, s, d) = pseudo_divrem_divconquer(g, f);
    tassert(r*f + s == g*fmpzxx(pow(2, d) + 0.5));

    tassert(pseudo_div(g, f).get<0>() == r);
    tassert(pseudo_rem(g, f).get<0>() == s);

    r = 0; r = 0;
    ltupleref(r, s) = pseudo_divrem_cohen(g, f);
    tassert(r*f + s == g*4);
    tassert(pseudo_rem_cohen(g, f).to_string() == "3  -28 -32 12");

    f = "4  1 0 0 1";
    slong r1, r2;
    f.signature(r1, r2);
    tassert(r1 == 1 && r2 == 1);

    // test static functions
    frandxx state;
    tassert(fmpz_polyxx::randtest(state, 4, 10).length() <= 4);
    tassert(fmpz_polyxx::randtest_unsigned(state, 4, 10).get_coeff(0) >= 0);
    tassert(fmpz_polyxx::randtest_not_zero(state, 4, 10).is_zero() == false);
    tassert(fmpz_polyxx::interpolate(xs, ys) == xp1);

    xs[0] = 0;xs[1] = -1; xs[2] = -1;
    tassert(fmpz_polyxx::product_roots(xs) == x*xp1*xp1);
}

void
test_member_functions()
{
    // just a sample, since they all come from macros
    fmpz_polyxx f, g;
    f = "4  2 0 0 1";
    g = "5  1 2 3 4 5";

    tassert(f.bit_pack(17u) == bit_pack(f, 17u));
    tassert(f.divrem(g) == divrem(f, g));
    tassert(f.derivative() == derivative(f));
    tassert(f.bound_roots() == bound_roots(f));
    tassert(f.content() == content(f));
    tassert(f.pow_trunc(15u, 10) == pow_trunc(f, 15u, 10));
}

void
test_factoring()
{
    fmpz_polyxx f, g;
    // two irreducible polynomials
    f = "4  1 17 0 1";
    g = "6  2 0 2 0 0 1";

    // TODO are these deterministic?
    fmpz_poly_factorxx f1, f2;
    f1.insert(f, 1);
    f2.insert(g, 2);
    f1.content() = 7;
    f1.concat(f2);

    tassert(f1 == factor_zassenhaus(7*f*g*g));

    f1.realloc(0);
    f1.insert(f, 1);
    f1.insert(g, 2);
    f1.content() = 1;
    tassert(f1 == factor_squarefree(f*g*g));

    // TODO test set_factor_zassenhaus_recombination

    if(0)
        print(f1); // make sure this compiles
}

void
test_hensel()
{
    mp_limb_t pl = 1031;
    frandxx state;
    nmod_polyxx gl(nmod_polyxx::randtest_irreducible(pl, state, 5).make_monic());
    nmod_polyxx hl(nmod_polyxx::randtest_irreducible(pl, state, 6).make_monic());
    while(gl.length() != 5)
        gl = nmod_polyxx::randtest_irreducible(pl, state, 5).make_monic();
    while(hl.length() != 5)
        hl = nmod_polyxx::randtest_irreducible(pl, state, 5).make_monic();
    nmod_polyxx al(pl), bl(pl);
    ltupleref(_, al, bl) = xgcd(gl, hl);
    tassert((al*gl + bl*hl).is_one());

    fmpz_polyxx f = fmpz_polyxx::lift(gl*hl);
    fmpz_polyxx g = fmpz_polyxx::lift(gl);
    fmpz_polyxx h = fmpz_polyxx::lift(hl);
    fmpz_polyxx a = fmpz_polyxx::lift(al);
    fmpz_polyxx b = fmpz_polyxx::lift(bl);

    fmpz_polyxx G, H, A, B;
    fmpzxx p(pl), p1(1031);
    tassert(((f - g*h) % p).is_zero());

    //ltupleref(G, H, A, B) = hensel_lift(f, g, h, a, b, p, p1);
    fmpz_poly_hensel_lift(G._poly(), H._poly(), A._poly(), B._poly(), f._poly(), g._poly(), h._poly(), a._poly(), b._poly(), p._fmpz(), p1._fmpz());
    tassert(((f - G*H) % (p1*p)).is_zero());
    tassert(((A*G + B*H) % (p1*p)).is_one());

    tassert(ltupleref(G, H) == hensel_lift_without_inverse(f, g, h, a, b, p, p1));
    tassert(ltupleref(A, B) == hensel_lift_only_inverse(G, H, a, b, p, p1));

    nmod_poly_factorxx local_fac = factor(gl*hl);
    fmpz_poly_factorxx lifted_fac = hensel_lift_once(f, local_fac, 3);
    tassert(lifted_fac.size() == 2
            && lifted_fac.exp(0) == 1 && lifted_fac.exp(1) == 1);
    tassert(((lifted_fac.p(0)*lifted_fac.p(1) - f) % p.pow(3u)).is_zero());
}

void
test_printing()
{
    frandxx state;
    fmpz_polyxx f = fmpz_polyxx::randtest(state, 4, 10);
    test_print_read(f);
    test_print_read_pretty(f);
}

void
test_unified_access()
{
    fmpz_polyxx a = fmpz_polyxx::from_ground(1);
    const fmpz_polyxx& b = a;
    tassert(b.lead() == 1 && b.coeff(0) == 1);
}

int
main()
{
    std::cout << "fmpz_polyxx....";

    test_init();
    test_manipulation();
    test_assignment();
    test_conversion();
    test_arithmetic();
    test_functions();
    test_member_functions();
    test_extras();
    test_factoring();
    test_hensel();
    test_printing();
    test_unified_access();

    std::cout << "PASS" << std::endl;
    return 0;
}

