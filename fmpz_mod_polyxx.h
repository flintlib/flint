/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#ifndef FMPZ_MOD_POLYXX_H
#define FMPZ_MOD_POLYXX_H

#include "fmpz_mod_poly.h"

#include "fmpzxx.h"
#include "fmpz_polyxx.h"
#include "nmod_polyxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/flint_exception.h"
#include "flintxx/frandxx.h"
#include "flintxx/ltuple.h"
#include "flintxx/stdmath.h"
#include "flintxx/vector.h"

// TODO create/use fmpz_modxx class?

namespace flint {
FLINT_DEFINE_BINOP(divrem_f)
FLINT_DEFINE_BINOP(gcd_euclidean_f)
FLINT_DEFINE_BINOP(gcd_f)
FLINT_DEFINE_BINOP(radix)

FLINT_DEFINE_UNOP(fmpz_mod_polyxx_lead) // TODO standardise?
FLINT_DEFINE_BINOP(fmpz_mod_polyxx_get_coeff) // TODO standardise?

namespace detail {
template<class Poly>
struct fmpz_mod_poly_traits
{
    typedef FLINT_UNOP_BUILD_RETTYPE(
                fmpz_mod_polyxx_lead, fmpzxx, Poly) lead_ref_t;
    typedef lead_ref_t lead_srcref_t;
    static lead_ref_t lead(const Poly& p) {return fmpz_mod_polyxx_lead(p);}
};
}

template<class Operation, class Data>
class fmpz_mod_polyxx_expression
    : public expression<derived_wrapper<fmpz_mod_polyxx_expression>,
                            Operation, Data>
{
    typedef expression<derived_wrapper< ::flint::fmpz_mod_polyxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_mod_poly_traits<fmpz_mod_polyxx_expression>
        poly_traits_t;

    FLINTXX_DEFINE_BASICS(fmpz_mod_polyxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_mod_polyxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_mod_polyxx_expression, fmpz_mod_poly_struct, _poly)

    template<class Fmpz>
    static fmpz_mod_polyxx_expression randtest(const Fmpz& m,
            frandxx& state, slong len)
    {
        fmpz_mod_polyxx_expression res(m);
        res.set_randtest(state, len);
        return res;
    }
    template<class Fmpz>
    static fmpz_mod_polyxx_expression randtest_irreducible(const Fmpz& m,
            frandxx& state, slong len)
    {
        fmpz_mod_polyxx_expression res(m);
        res.set_randtest_irreducible(state, len);
        return res;
    }
    template<class Fmpz>
    static fmpz_mod_polyxx_expression randtest_not_zero(const Fmpz& m,
            frandxx& state, slong len)
    {
        fmpz_mod_polyxx_expression res(m);
        res.set_randtest_not_zero(state, len);
        return res;
    }

    template<class Fmpz>
    static fmpz_mod_polyxx_expression zero(const Fmpz& m)
        {return fmpz_mod_polyxx_expression(m);}

    // these only make sense with immediates
    fmpzxx_srcref _mod() const
        {return fmpzxx_srcref::make(fmpz_mod_poly_modulus(_poly()));}

    // These only make sense with target immediates
    void realloc(slong alloc) {fmpz_mod_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {fmpz_mod_poly_fit_length(_poly(), len);}
    void _normalise() {_fmpz_mod_poly_normalise(_poly());}
    void set_coeff(slong n, ulong c) {fmpz_mod_poly_set_coeff_ui(_poly(), n, c);}
    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_coeff(slong j, const Fmpz& c)
    {
        fmpz_mod_poly_set_coeff_fmpz(_poly(), j, c.evaluate()._fmpz());
    }
    void truncate(slong n) {fmpz_mod_poly_truncate(_poly(), n);}
    void zero_coeffs(slong i, slong j) {fmpz_mod_poly_zero_coeffs(_poly(), i, j);}

    void set_randtest(frandxx& state, slong len)
        {fmpz_mod_poly_randtest(_poly(), state._data(), len);}
    void set_randtest_irreducible(frandxx& state, slong len)
        {fmpz_mod_poly_randtest_irreducible(_poly(), state._data(), len);}
    void set_randtest_not_zero(frandxx& state, slong len)
        {fmpz_mod_poly_randtest_not_zero(_poly(), state._data(), len);}

    template<class Poly>
    slong remove(const Poly& p)
    {
        return fmpz_mod_poly_remove(_poly(), p.evaluate()._poly());
    }

    void set_zero() {fmpz_mod_poly_zero(_poly());}

    // unified coefficient access
    typename poly_traits_t::lead_ref_t lead()
    {
        return poly_traits_t::lead(*this);
    }
    typename poly_traits_t::lead_srcref_t lead() const
    {
        return poly_traits_t::lead(*this);
    }

    // this works without evaluation
    fmpzxx_srcref modulus() const;

    evaluated_t create_temporary() const
    {
        return evaluated_t(modulus());
    }

    // These cause evaluation
    slong length() const {return fmpz_mod_poly_length(this->evaluate()._poly());}
    slong degree() const {return fmpz_mod_poly_degree(this->evaluate()._poly());}
    bool is_zero() const {return fmpz_mod_poly_is_zero(this->evaluate()._poly());}
    bool is_squarefree() const
        {return fmpz_mod_poly_is_squarefree(this->evaluate()._poly());}
    bool is_irreducible() const
        {return fmpz_mod_poly_is_irreducible(this->evaluate()._poly());}
    bool is_irreducible_ddf() const
        {return fmpz_mod_poly_is_irreducible_ddf(this->evaluate()._poly());}
    bool is_irreducible_rabin() const
        {return fmpz_mod_poly_is_irreducible_rabin(this->evaluate()._poly());}

    // Lazy members
    FLINTXX_DEFINE_MEMBER_BINOP_(get_coeff, fmpz_mod_polyxx_get_coeff)
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)

    FLINTXX_DEFINE_MEMBER_UNOP(derivative)
    FLINTXX_DEFINE_MEMBER_UNOP(invmod)
    FLINTXX_DEFINE_MEMBER_UNOP(make_monic)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr)

    FLINTXX_DEFINE_MEMBER_BINOP(compose_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(compose_horner)
    FLINTXX_DEFINE_MEMBER_BINOP(div_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem_f)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_euclidean)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_euclidean_f)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_f)
    FLINTXX_DEFINE_MEMBER_BINOP(gcdinv)
    FLINTXX_DEFINE_MEMBER_BINOP(invmod)
    FLINTXX_DEFINE_MEMBER_BINOP(inv_series_newton)
    FLINTXX_DEFINE_MEMBER_BINOP(shift_left)
    FLINTXX_DEFINE_MEMBER_BINOP(shift_right)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_BINOP(radix)
    FLINTXX_DEFINE_MEMBER_BINOP(rem_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(xgcd)
    FLINTXX_DEFINE_MEMBER_BINOP(xgcd_euclidean)

    FLINTXX_DEFINE_MEMBER_3OP(compose_mod)
    FLINTXX_DEFINE_MEMBER_3OP(compose_mod_brent_kung)
    FLINTXX_DEFINE_MEMBER_3OP(compose_mod_horner)
    FLINTXX_DEFINE_MEMBER_3OP(mullow)
    FLINTXX_DEFINE_MEMBER_3OP(mulmod)
    FLINTXX_DEFINE_MEMBER_3OP(powmod_binexp)
    FLINTXX_DEFINE_MEMBER_3OP(pow_trunc)
    FLINTXX_DEFINE_MEMBER_3OP(pow_trunc_binexp)
};

namespace detail {
struct fmpz_mod_poly_data;
}

typedef fmpz_mod_polyxx_expression<operations::immediate, detail::fmpz_mod_poly_data>
           fmpz_mod_polyxx;
typedef fmpz_mod_polyxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_mod_polyxx, fmpz_mod_poly_struct> >
           fmpz_mod_polyxx_ref;
typedef fmpz_mod_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_mod_polyxx, fmpz_mod_polyxx_ref, fmpz_mod_poly_struct> >
           fmpz_mod_polyxx_srcref;

#define FMPZ_MOD_POLYXX_COND_S FLINTXX_COND_S(fmpz_mod_polyxx)
#define FMPZ_MOD_POLYXX_COND_T FLINTXX_COND_T(fmpz_mod_polyxx)

namespace detail {
template<>
struct fmpz_mod_poly_traits<fmpz_mod_polyxx_srcref>
{
    typedef fmpzxx_srcref lead_srcref_t;
    typedef fmpzxx_srcref lead_ref_t;

    template<class P>
    static lead_srcref_t lead(P p)
        {return lead_srcref_t::make(fmpz_mod_poly_lead(p._poly()));}
};
template<>
struct fmpz_mod_poly_traits<fmpz_mod_polyxx_ref>
{
    typedef fmpzxx_ref lead_ref_t;
    typedef fmpzxx_ref lead_srcref_t;

    template<class P>
    static lead_ref_t lead(P p)
        {return lead_ref_t::make(fmpz_mod_poly_lead(p._poly()));}
};
template<>
struct fmpz_mod_poly_traits<fmpz_mod_polyxx>
{
    typedef fmpzxx_ref lead_ref_t;
    typedef fmpzxx_srcref lead_srcref_t;

    template<class P>
    static lead_ref_t lead(P& p)
        {return lead_ref_t::make(fmpz_mod_poly_lead(p._poly()));}
    template<class P>
    static lead_srcref_t lead(const P& p)
        {return lead_srcref_t::make(fmpz_mod_poly_lead(p._poly()));}
};

struct fmpz_mod_poly_data
{
    fmpz_mod_poly_t inner;
    typedef fmpz_mod_poly_t& data_ref_t;
    typedef const fmpz_mod_poly_t& data_srcref_t;

    template<class Fmpz>
    fmpz_mod_poly_data(const Fmpz& n,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpz_mod_poly_init(inner, n.evaluate()._fmpz());
    }
    template<class Fmpz>
    fmpz_mod_poly_data(const Fmpz& n, slong alloc,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpz_mod_poly_init2(inner, n.evaluate()._fmpz(), alloc);
    }
    ~fmpz_mod_poly_data() {fmpz_mod_poly_clear(inner);}

    fmpz_mod_poly_data(const fmpz_mod_poly_data& o)
    {
        fmpz_mod_poly_init(inner, fmpz_mod_poly_modulus(o.inner));
        fmpz_mod_poly_set(inner, o.inner);
    }

    fmpz_mod_poly_data(fmpz_mod_polyxx_srcref r)
    {
        fmpz_mod_poly_init(inner, r.modulus()._fmpz());
        fmpz_mod_poly_set(inner, r._poly());
    }
};

struct is_fmpz_mod_polyxx_predicate
{
    template<class T> struct type : FMPZ_MOD_POLYXX_COND_S<T> { };
};
}
template<class Operation, class Data>
inline fmpzxx_srcref
fmpz_mod_polyxx_expression<Operation, Data>::modulus() const
{
    return tools::find_subexpr<detail::is_fmpz_mod_polyxx_predicate>(
            *this)._mod();
}

namespace traits {
template<class T> struct is_fmpz_mod_polyxx : mp::or_<
     traits::is_T_expr<T, fmpz_mod_polyxx>,
     flint_classes::is_source<fmpz_mod_polyxx, T> > { };
}

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MOD_POLYXX_COND_T,
        FMPZ_MOD_POLYXX_COND_S, fmpz_mod_poly_set(to._poly(), from._poly()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MOD_POLYXX_COND_T,
        traits::is_unsigned_integer, fmpz_mod_poly_set_ui(to._poly(), from))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MOD_POLYXX_COND_T,
        FMPZXX_COND_S, fmpz_mod_poly_set_fmpz(to._poly(), from._fmpz()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MOD_POLYXX_COND_T, FMPZ_POLYXX_COND_S,
        fmpz_mod_poly_set_fmpz_poly(to._poly(), from._poly()))
FLINTXX_DEFINE_CONVERSION_TMP(fmpz_polyxx, fmpz_mod_polyxx,
        fmpz_mod_poly_get_fmpz_poly(to._poly(), from._poly()))

FLINTXX_DEFINE_SWAP(fmpz_mod_polyxx, fmpz_mod_poly_swap(e1._poly(), e2._poly()))

FLINTXX_DEFINE_EQUALS(fmpz_mod_polyxx, fmpz_mod_poly_equal(e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpz_mod_polyxx_get_coeff_op, fmpzxx,
        FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_mod_poly_get_coeff_fmpz(to._fmpz(), e1._poly(), e2))

FLINT_DEFINE_PRINT_COND(FMPZ_MOD_POLYXX_COND_S, fmpz_mod_poly_fprint(to, from._poly()))
FLINT_DEFINE_PRINT_PRETTY_COND_2(FMPZ_MOD_POLYXX_COND_S, const char*,
        fmpz_mod_poly_fprint_pretty(to, from._poly(), extra))
FLINT_DEFINE_READ_COND(FMPZ_MOD_POLYXX_COND_T, fmpz_mod_poly_fread(from, to._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_add(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_sub(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_mod_poly_scalar_mul_fmpz(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_mul(to._poly(), e1._poly(), e2._poly()))

// TODO expose the temporary
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_polyxx tmp(to.modulus());
        fmpz_mod_poly_divrem(to._poly(), tmp._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_mod_polyxx, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_neg(to._poly(), from._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(fmpz_mod_polyxx_lead_op, fmpzxx,
        FMPZ_MOD_POLYXX_COND_S,
        fmpz_set(to._fmpz(), fmpz_mod_poly_lead(from._poly())))

FLINT_DEFINE_BINARY_EXPR_COND2(shift_left_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_mod_poly_shift_left(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(shift_right_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_mod_poly_shift_right(to._poly(), e1._poly(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(make_monic_op, fmpz_mod_polyxx, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_make_monic(to._poly(), from._poly()))

FLINT_DEFINE_THREEARY_EXPR_COND3(mullow_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
    fmpz_mod_poly_mullow(to._poly(), e1._poly(), e2._poly(), e3))
FLINT_DEFINE_THREEARY_EXPR_COND3(mulmod_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_mulmod(to._poly(), e1._poly(), e2._poly(), e3._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_op, fmpz_mod_polyxx, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_sqr(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_rem(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_THREEARY_EXPR_COND3(powmod_binexp_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, traits::is_unsigned_integer, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_powmod_ui_binexp(to._poly(), e1._poly(), e2, e3._poly()))
FLINT_DEFINE_THREEARY_EXPR_COND3(powmod_binexp_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, FMPZXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_powmod_fmpz_binexp(
        to._poly(), e1._poly(), e2._fmpz(), e3._poly()))

FLINT_DEFINE_THREEARY_EXPR_COND3(pow_trunc_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, traits::is_unsigned_integer, traits::fits_into_slong,
    fmpz_mod_poly_pow_trunc(to._poly(), e1._poly(), e2, e3))
FLINT_DEFINE_THREEARY_EXPR_COND3(pow_trunc_binexp_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, traits::is_unsigned_integer, traits::fits_into_slong,
    fmpz_mod_poly_pow_trunc_binexp(to._poly(), e1._poly(), e2, e3))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, traits::is_unsigned_integer,
    fmpz_mod_poly_pow(to._poly(), e1._poly(), e2))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<fmpz_mod_polyxx, fmpz_mod_polyxx>::type>::type
    fmpz_mod_polyxx_pair;
typedef make_ltuple<mp::make_tuple<
        fmpzxx, fmpz_mod_polyxx, fmpz_mod_polyxx>::type>::type
    fmpz_mod_poly_divrem_f_rt;
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(divrem_basecase_op, rdetail::fmpz_mod_polyxx_pair,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_divrem_basecase(
        to.template get<0>()._poly(), to.template get<1>()._poly(),
        e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(divrem_divconquer_op, rdetail::fmpz_mod_polyxx_pair,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_divrem_divconquer(
        to.template get<0>()._poly(), to.template get<1>()._poly(),
        e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(divrem_op, rdetail::fmpz_mod_polyxx_pair,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_divrem_divconquer(
        to.template get<0>()._poly(), to.template get<1>()._poly(),
        e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(divrem_f_op, rdetail::fmpz_mod_poly_divrem_f_rt,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_divrem_f(
        to.template get<0>()._fmpz(), to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(div_basecase_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_div_basecase(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(rem_basecase_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_rem_basecase(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(inv_series_newton_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_mod_poly_inv_series_newton(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(gcd_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_gcd(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(gcd_euclidean_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_gcd_euclidean(to._poly(), e1._poly(), e2._poly()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<
        fmpz_mod_polyxx, fmpz_mod_polyxx, fmpz_mod_polyxx>::type>::type
    fmpz_mod_polyxx_triple;
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(xgcd_op, rdetail::fmpz_mod_polyxx_triple,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_xgcd(to.template get<0>()._poly(), to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(xgcd_euclidean_op, rdetail::fmpz_mod_polyxx_triple,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_xgcd_euclidean(to.template get<0>()._poly(),
        to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<fmpzxx, fmpz_mod_polyxx>::type>::type
    fmpz_mod_gcd_f_rt;
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(gcd_f_op, rdetail::fmpz_mod_gcd_f_rt,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_gcd_f(to.template get<0>()._fmpz(),
            to.template get<1>()._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(gcd_euclidean_f_op, rdetail::fmpz_mod_gcd_f_rt,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_gcd_euclidean_f(to.template get<0>()._fmpz(),
            to.template get<1>()._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(gcdinv_op, rdetail::fmpz_mod_polyxx_pair,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    fmpz_mod_poly_gcdinv(
        to.template get<0>()._poly(), to.template get<1>()._poly(),
        e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(invmod_op, fmpz_mod_polyxx,
    FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
    execution_check(fmpz_mod_poly_invmod(to._poly(), e1._poly(), e2._poly()),
        "invmod", "fmpz_mod_polyxx"))

FLINT_DEFINE_UNARY_EXPR_COND(derivative_op, fmpz_mod_polyxx, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_derivative(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpzxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_mod_poly_evaluate_fmpz(to._fmpz(), e1._poly(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(compose_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(compose_divconquer_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose_divconquer(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(compose_horner_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose_horner(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_THREEARY_EXPR_COND3(compose_mod_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose_mod(to._poly(), e1._poly(), e2._poly(), e3._poly()))
FLINT_DEFINE_THREEARY_EXPR_COND3(compose_mod_horner_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose_mod_horner(
            to._poly(), e1._poly(), e2._poly(), e3._poly()))
FLINT_DEFINE_THREEARY_EXPR_COND3(compose_mod_brent_kung_op, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_compose_mod_brent_kung(
            to._poly(), e1._poly(), e2._poly(), e3._poly()))
} // rules


///////////////////////////////////////////////////////////////////////////////
// fmpz_mod_poly_vecxx (for radix conversion)
///////////////////////////////////////////////////////////////////////////////
namespace detail {
struct fmpz_mod_poly_vector_data
{
    slong size;
    fmpz_mod_poly_struct** array;

    template<class Fmpz>
    void init(slong n, const Fmpz& f)
    {
        size = n;
        array = new fmpz_mod_poly_struct*[n];
        for(slong i = 0;i < n;++i)
        {
            array[i] = new fmpz_mod_poly_struct();
            fmpz_mod_poly_init(array[i], f._fmpz());
        }
    }

    template<class Fmpz>
    fmpz_mod_poly_vector_data(slong n, const Fmpz& f,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        init(n, f.evaluate());
    }

    ~fmpz_mod_poly_vector_data()
    {
        for(slong i = 0;i < size;++i)
        {
            fmpz_mod_poly_clear(array[i]);
            delete array[i];
        }
        delete[] array;
    }

    fmpz_mod_poly_vector_data(const fmpz_mod_poly_vector_data& o)
        : size(o.size)
    {
        array = new fmpz_mod_poly_struct*[size];
        for(slong i = 0;i < size;++i)
        {
            array[i] = new fmpz_mod_poly_struct();
            fmpz_mod_poly_init(array[i], &o.array[0]->p);
            fmpz_mod_poly_set(array[i], o.array[i]);
        }
    }

    fmpz_mod_polyxx_ref at(slong i)
        {return fmpz_mod_polyxx_ref::make(array[i]);}
    fmpz_mod_polyxx_srcref at(slong i) const
        {return fmpz_mod_polyxx_srcref::make(array[i]);}

    bool equals(const fmpz_mod_poly_vector_data& o) const
    {
        if(size != o.size)
            return false;
        for(slong i = 0;i < size;++i)
            if(!fmpz_mod_poly_equal(array[i], o.array[i]))
                return false;
        return true;
    }
};

// TODO extend this somewhat similarly to the nmod* setup
template<class T>
fmpzxx_srcref find_fmpz_mod_polyxx_mod(const T& t)
{
    // this works because there are not actually any functions operating on
    // fmpz_mod_poly_vecxx
    return tools::find_subexpr<detail::is_fmpz_mod_polyxx_predicate>(
            t)._mod();
}

struct fmpz_mod_poly_vector_traits
    : wrapped_vector_traits<fmpz_mod_polyxx, slong, fmpz_mod_polyxx_ref,
          fmpz_mod_polyxx_srcref, fmpz_mod_poly_struct*>
{
    template<class Expr>
    static typename Expr::evaluated_t create_temporary(const Expr& e)
    {
        return typename Expr::evaluated_t(
                e.size(), find_fmpz_mod_polyxx_mod(e));
    }
};
} // detail

// TODO would it make more sense to have this have its own class?
typedef vector_expression<
    detail::fmpz_mod_poly_vector_traits, operations::immediate,
    detail::fmpz_mod_poly_vector_data> fmpz_mod_poly_vecxx;
// TODO references

template<>
struct enable_vector_rules<fmpz_mod_poly_vecxx> : mp::false_ { };

namespace rules {
// TODO hack to make code look like references are implemented
template<class T> struct FMPZ_MOD_POLY_VECXX_COND_S
    : mp::equal_types<T, fmpz_mod_poly_vecxx> { };
#define FMPZ_MOD_POLY_VECXX_COND_T FMPZ_MOD_POLY_VECXX_COND_S

// TODO references
FLINT_DEFINE_GET(equals, bool, fmpz_mod_poly_vecxx, e1._data().equals(e2._data()))
} // rules


///////////////////////////////////////////////////////////////////////////////
// radix conversion
///////////////////////////////////////////////////////////////////////////////

class fmpz_mod_poly_radixxx
{
private:
    fmpz_mod_poly_radix_t inner;

    // not copyable
    fmpz_mod_poly_radixxx(const fmpz_mod_poly_radixxx&);

public:
    template<class Fmpz_mod_poly>
    fmpz_mod_poly_radixxx(const Fmpz_mod_poly& r, slong deg,
            typename mp::enable_if<
                traits::is_fmpz_mod_polyxx<Fmpz_mod_poly> >::type* = 0)
        {fmpz_mod_poly_radix_init(inner, r.evaluate()._poly(), deg);}
    ~fmpz_mod_poly_radixxx() {fmpz_mod_poly_radix_clear(inner);}

    fmpz_mod_poly_radix_t& _data() {return inner;}
    const fmpz_mod_poly_radix_t& _data() const {return inner;}

    slong degR() const {return inner->degR;}
};

namespace traits {
template<class T> struct is_fmpz_mod_poly_radixxx
   : mp::equal_types<T, fmpz_mod_poly_radixxx> { };
} // traits

namespace vectors {
template<>
struct outsize<operations::radix_op>
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        return e._data().first().degree() / e._data().second().degR() + 1;
    }
};
}

namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(radix_op, fmpz_mod_poly_vecxx,
        FMPZ_MOD_POLYXX_COND_S, traits::is_fmpz_mod_poly_radixxx,
        fmpz_mod_poly_radix(to._array(), e1._poly(), e2._data()))
}
} // flint

#include "fmpz_mod_poly_factorxx.h"

#endif
