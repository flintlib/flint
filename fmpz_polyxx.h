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

#ifndef FMPZ_POLYXX_H
#define FMPZ_POLYXX_H

#include <cstdlib>
#include <string>

#include "fmpz_poly.h"

#include "fmpzxx.h"
#include "fmpz_vecxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/flint_exception.h"
#include "flintxx/frandxx.h"
#include "flintxx/ltuple.h"
#include "flintxx/traits.h"

// TODO exhibit this as a specialisation of a generic poly<fmpzxx>
// TODO newton basis?
// TODO power series class?
// TODO hensel lifting
// TODO input
// TODO modular reduction
// TODO addmul
// TODO rename poly_divides to divides?

namespace flint {
// function "declarations"
FLINT_DEFINE_UNOP(sqr_classical)
FLINT_DEFINE_UNOP(sqr_karatsuba)
FLINT_DEFINE_UNOP(sqr_KS)
FLINT_DEFINE_UNOP(sqrt_classical)

FLINT_DEFINE_BINOP(compose_divconquer)
FLINT_DEFINE_BINOP(compose_horner)
FLINT_DEFINE_BINOP(div_basecase)
FLINT_DEFINE_BINOP(div_divconquer)
FLINT_DEFINE_BINOP(divrem_basecase)
FLINT_DEFINE_BINOP(divrem_divconquer)
FLINT_DEFINE_BINOP(evaluate_divconquer)
FLINT_DEFINE_BINOP(evaluate_horner)
FLINT_DEFINE_BINOP(gcd_heuristic)
FLINT_DEFINE_BINOP(gcd_modular)
FLINT_DEFINE_BINOP(gcd_subresultant)
FLINT_DEFINE_BINOP(mul_karatsuba)
FLINT_DEFINE_BINOP(mulmid_classical)
FLINT_DEFINE_BINOP(mul_SS)
FLINT_DEFINE_BINOP(poly_divides)
FLINT_DEFINE_BINOP(pow_addchains)
FLINT_DEFINE_BINOP(pow_binomial)
FLINT_DEFINE_BINOP(pow_multinomial)
FLINT_DEFINE_BINOP(pseudo_div)
FLINT_DEFINE_BINOP(pseudo_divrem)
FLINT_DEFINE_BINOP(pseudo_divrem_basecase)
FLINT_DEFINE_BINOP(pseudo_divrem_cohen)
FLINT_DEFINE_BINOP(pseudo_divrem_divconquer)
FLINT_DEFINE_BINOP(pseudo_rem)
FLINT_DEFINE_BINOP(pseudo_rem_cohen)
FLINT_DEFINE_BINOP(rem_basecase)
FLINT_DEFINE_BINOP(sqrlow)
FLINT_DEFINE_BINOP(sqrlow_classical)
FLINT_DEFINE_BINOP(sqrlow_karatsuba_n)
FLINT_DEFINE_BINOP(sqrlow_KS)
FLINT_DEFINE_BINOP(taylor_shift_divconquer)
FLINT_DEFINE_BINOP(taylor_shift_horner)
FLINT_DEFINE_BINOP(xgcd_modular)

FLINT_DEFINE_THREEARY(compose_series)
FLINT_DEFINE_THREEARY(compose_series_brent_kung)
FLINT_DEFINE_THREEARY(compose_series_horner)
FLINT_DEFINE_THREEARY(div_series)
FLINT_DEFINE_THREEARY(mulhigh_karatsuba_n)
FLINT_DEFINE_THREEARY(mulhigh_n)
FLINT_DEFINE_THREEARY(mullow_karatsuba_n)
FLINT_DEFINE_THREEARY(mullow_SS)


FLINT_DEFINE_BINOP(fmpz_polyxx_interpolate)
FLINT_DEFINE_UNOP(fmpz_polyxx_product_roots)
FLINT_DEFINE_UNOP(fmpz_polyxx_lead)
FLINT_DEFINE_BINOP(fmpz_polyxx_get_coeff)

namespace traits {
template<class T> struct is_fmpz_polyxx;
} // traits

namespace detail {
template<class Poly>
struct fmpz_poly_traits
{
    typedef FLINT_UNOP_BUILD_RETTYPE(
                fmpz_polyxx_lead, fmpzxx, Poly) lead_ref_t;
    typedef lead_ref_t lead_srcref_t;
    static lead_ref_t lead(const Poly& p) {return fmpz_polyxx_lead(p);}

    template<class T>
    struct get_coeff
    {
        typedef FLINT_BINOP_ENABLE_RETTYPE(fmpz_polyxx_get_coeff, Poly, T) ref_t;
        typedef ref_t srcref_t;
        static ref_t get(const Poly& p, const T& t)
            {return fmpz_polyxx_get_coeff(p, t);}
    };
};
}

template<class Operation, class Data>
class fmpz_polyxx_expression
    : public expression<derived_wrapper<fmpz_polyxx_expression>,
                            Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_polyxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_poly_traits<fmpz_polyxx_expression> poly_traits_t;

    FLINTXX_DEFINE_BASICS(fmpz_polyxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_polyxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_polyxx_expression, fmpz_poly_struct, _poly)

    // static methods which only make sense with fmpz_polyxx
    static fmpz_polyxx_expression randtest(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpz_polyxx_expression res;
        fmpz_poly_randtest(res._poly(), state._data(), len, bits);
        return res;
    }
    static fmpz_polyxx_expression randtest_unsigned(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpz_polyxx_expression res;
        fmpz_poly_randtest_unsigned(res._poly(), state._data(), len, bits);
        return res;
    }
    static fmpz_polyxx_expression randtest_not_zero(frandxx& state, slong len,
            mp_bitcnt_t bits)
    {
        fmpz_polyxx_expression res;
        fmpz_poly_randtest_not_zero(res._poly(), state._data(), len, bits);
        return res;
    }

    template<class Fmpz_vec1, class Fmpz_vec2>
    static FLINT_BINOP_ENABLE_RETTYPE(fmpz_polyxx_interpolate,
        Fmpz_vec1, Fmpz_vec2)
    interpolate(const Fmpz_vec1& xs, const Fmpz_vec2& ys)
    {
        return fmpz_polyxx_interpolate(xs, ys);
    }

    template<class Fmpz_vec>
    static FLINT_UNOP_ENABLE_RETTYPE(fmpz_polyxx_product_roots, Fmpz_vec)
    product_roots(const Fmpz_vec& xs)
    {
        return fmpz_polyxx_product_roots(xs);
    }

    // These only make sense with immediates
    void realloc(slong alloc) {fmpz_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {fmpz_poly_fit_length(_poly(), len);}
    void _normalise() {_fmpz_poly_normalise(_poly());}
    void _set_length(slong len) {_fmpz_poly_set_length(_poly(), len);}
    void zero_coeffs(slong i, slong j) {fmpz_poly_zero_coeffs(_poly(), i, j);}

    // The result of these are undefined if n is >= length
    // You also may have to call _normalise().
    template<class T>
    typename poly_traits_t::template get_coeff<T>::ref_t get_coeff(const T& n)
    {
        return poly_traits_t::template get_coeff<T>::get(*this, n);
    }
    template<class T>
    typename poly_traits_t::template get_coeff<T>::srcref_t
        get_coeff(const T& n) const
    {
        return poly_traits_t::template get_coeff<T>::get(*this, n);
    }
    typename poly_traits_t::lead_ref_t lead()
    {
        return poly_traits_t::lead(*this);
    }
    typename poly_traits_t::lead_srcref_t lead() const
    {
        return poly_traits_t::lead(*this);
    }

    // These only make sense with target immediates
    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_coeff(slong n, const Fmpz& x)
    {
        fmpz_poly_set_coeff_fmpz(_poly(), n, x.evaluate()._fmpz());
    }
    template<class T>
    typename mp::enable_if<traits::is_signed_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpz_poly_set_coeff_si(_poly(), n, x);
    }
    template<class T>
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpz_poly_set_coeff_ui(_poly(), n, x);
    }

    void truncate(slong n) {fmpz_poly_truncate(_poly(), n);}

    // These cause evaluation
    slong length() const {return fmpz_poly_length(this->evaluate()._poly());}
    slong degree() const {return fmpz_poly_degree(this->evaluate()._poly());}
    bool is_one() const {return fmpz_poly_is_one(this->evaluate()._poly());}
    bool is_zero() const {return fmpz_poly_is_zero(this->evaluate()._poly());}
    bool is_unit() const {return fmpz_poly_is_unit(this->evaluate()._poly());}
    ulong max_limbs() const {return fmpz_poly_max_limbs(this->evaluate()._poly());}
    slong max_bits() const {return fmpz_poly_max_bits(this->evaluate()._poly());}

    std::string pretty(const char* x) const
    {
        char* str = fmpz_poly_get_str_pretty(this->evaluate()._poly(), x);
        std::string res(str);
        std::free(str);
        return res;
    }

    void signature(slong& r1, slong& r2) const
    {
        fmpz_poly_signature(&r1, &r2, this->evaluate()._poly());
    }

    // lazy member forwarding
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)
    FLINTXX_DEFINE_MEMBER_BINOP_(bit_pack, poly_bit_pack)
    FLINTXX_DEFINE_MEMBER_BINOP_(divides, poly_divides)

    FLINTXX_DEFINE_MEMBER_BINOP(compose_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(compose_horner)
    FLINTXX_DEFINE_MEMBER_BINOP(div_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(div_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(divexact)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(divrem_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(div_root)
    FLINTXX_DEFINE_MEMBER_BINOP(evaluate_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(evaluate_horner)
    FLINTXX_DEFINE_MEMBER_BINOP(fdiv_2exp)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_heuristic)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_modular)
    FLINTXX_DEFINE_MEMBER_BINOP(gcd_subresultant)
    FLINTXX_DEFINE_MEMBER_BINOP(inv_series)
    FLINTXX_DEFINE_MEMBER_BINOP(inv_series_newton)
    FLINTXX_DEFINE_MEMBER_BINOP(lcm)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_2exp)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_karatsuba)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_KS)
    FLINTXX_DEFINE_MEMBER_BINOP(mulmid_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_SS)
    FLINTXX_DEFINE_MEMBER_BINOP(poly_shift_left)
    FLINTXX_DEFINE_MEMBER_BINOP(poly_shift_right)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)
    FLINTXX_DEFINE_MEMBER_BINOP(pow_addchains)
    FLINTXX_DEFINE_MEMBER_BINOP(pow_binexp)
    FLINTXX_DEFINE_MEMBER_BINOP(pow_binomial)
    FLINTXX_DEFINE_MEMBER_BINOP(pow_multinomial)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_div)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_divrem)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_divrem_basecase)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_divrem_cohen)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_divrem_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_rem)
    FLINTXX_DEFINE_MEMBER_BINOP(pseudo_rem_cohen)
    FLINTXX_DEFINE_MEMBER_BINOP(resultant)
    FLINTXX_DEFINE_MEMBER_BINOP(reverse)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_lagrange)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_lagrange_fast)
    FLINTXX_DEFINE_MEMBER_BINOP(revert_series_newton)
    FLINTXX_DEFINE_MEMBER_BINOP(smod)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrlow)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrlow_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrlow_karatsuba_n)
    FLINTXX_DEFINE_MEMBER_BINOP(sqrlow_KS)
    FLINTXX_DEFINE_MEMBER_BINOP(taylor_shift)
    FLINTXX_DEFINE_MEMBER_BINOP(taylor_shift_divconquer)
    FLINTXX_DEFINE_MEMBER_BINOP(taylor_shift_horner)
    FLINTXX_DEFINE_MEMBER_BINOP(tdiv)
    FLINTXX_DEFINE_MEMBER_BINOP(tdiv_2exp)
    FLINTXX_DEFINE_MEMBER_BINOP(xgcd)
    FLINTXX_DEFINE_MEMBER_BINOP(xgcd_modular)

    FLINTXX_DEFINE_MEMBER_UNOP(derivative)
    FLINTXX_DEFINE_MEMBER_UNOP(primitive_part)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_classical)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_karatsuba)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_KS)
    FLINTXX_DEFINE_MEMBER_UNOP(sqrt)
    FLINTXX_DEFINE_MEMBER_UNOP(sqrt_classical)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(fmpzxx, bound_roots, poly_bound_roots)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(fmpzxx, norm, poly_2norm)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, content)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, height)

    FLINTXX_DEFINE_MEMBER_3OP(compose_series)
    FLINTXX_DEFINE_MEMBER_3OP(compose_series_brent_kung)
    FLINTXX_DEFINE_MEMBER_3OP(compose_series_horner)
    FLINTXX_DEFINE_MEMBER_3OP(div_series)
    FLINTXX_DEFINE_MEMBER_3OP(mulhigh_classical)
    FLINTXX_DEFINE_MEMBER_3OP(mulhigh_karatsuba_n)
    FLINTXX_DEFINE_MEMBER_3OP(mulhigh_n)
    FLINTXX_DEFINE_MEMBER_3OP(mullow)
    FLINTXX_DEFINE_MEMBER_3OP(mullow_classical)
    FLINTXX_DEFINE_MEMBER_3OP(mullow_karatsuba_n)
    FLINTXX_DEFINE_MEMBER_3OP(mullow_KS)
    FLINTXX_DEFINE_MEMBER_3OP(mullow_SS)
    FLINTXX_DEFINE_MEMBER_3OP(pow_trunc)
};

namespace detail {
struct fmpz_poly_data;
}

typedef fmpz_polyxx_expression<operations::immediate, detail::fmpz_poly_data>
           fmpz_polyxx;
typedef fmpz_polyxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_polyxx, fmpz_poly_struct> >
           fmpz_polyxx_ref;
typedef fmpz_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_polyxx, fmpz_polyxx_ref, fmpz_poly_struct> >
           fmpz_polyxx_srcref;

namespace detail {
template<>
struct fmpz_poly_traits<fmpz_polyxx_srcref>
{
    typedef fmpzxx_srcref lead_srcref_t;
    typedef fmpzxx_srcref lead_ref_t;

    template<class P>
    static lead_srcref_t lead(const P& p)
        {return lead_srcref_t::make(fmpz_poly_lead(p._poly()));}

    template<class T>
    struct get_coeff
    {
        typedef typename mp::enable_if<
            traits::is_integer<T>, fmpzxx_srcref>::type ref_t;
        typedef ref_t srcref_t;

        template<class P>
        static srcref_t get(const P& p, const T& n)
            {return srcref_t::make(fmpz_poly_get_coeff_ptr(p._poly(), n));}
    };
};
template<>
struct fmpz_poly_traits<fmpz_polyxx_ref>
    : fmpz_poly_traits<fmpz_polyxx_srcref>
{
    typedef fmpzxx_ref lead_ref_t;

    template<class P>
    static lead_ref_t lead(P& p)
        {return lead_ref_t::make(fmpz_poly_lead(p._poly()));}

    template<class T>
    struct get_coeff
        : fmpz_poly_traits<fmpz_polyxx_srcref>::template get_coeff<T>
    {
        typedef fmpzxx_ref ref_t;

        template<class P>
        static ref_t get(P& p, const T& n)
            {return ref_t::make(fmpz_poly_get_coeff_ptr(p._poly(), n));}
    };
};
template<>
struct fmpz_poly_traits<fmpz_polyxx>
    : fmpz_poly_traits<fmpz_polyxx_ref>
{ };

struct fmpz_poly_data
{
    fmpz_poly_t inner;
    typedef fmpz_poly_t& data_ref_t;
    typedef const fmpz_poly_t& data_srcref_t;

    fmpz_poly_data() {fmpz_poly_init(inner);}
    ~fmpz_poly_data() {fmpz_poly_clear(inner);}

    fmpz_poly_data(const fmpz_poly_data& o)
    {
        fmpz_poly_init(inner);
        fmpz_poly_set(inner, o.inner);
    }

    fmpz_poly_data(fmpz_polyxx_srcref r)
    {
        fmpz_poly_init(inner);
        fmpz_poly_set(inner, r._poly());
    }

    fmpz_poly_data(slong alloc)
    {
        fmpz_poly_init2(inner, alloc);
    }
};
} // detail

namespace traits {
template<class T> struct is_fmpz_polyxx : mp::or_<
     traits::is_T_expr<T, fmpz_polyxx>,
     flint_classes::is_source<fmpz_polyxx, T> > { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpz_polyxx : mp::and_<all_fmpz_polyxx<T1>, all_fmpz_polyxx<T2, T3, T4> > { };
template<class T>
struct all_fmpz_polyxx<T, void, void, void> : traits::is_fmpz_polyxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpz_polyxx
    : mp::enable_if<all_fmpz_polyxx<T1, T2, T3, T4>, Out> { };
} // mp

namespace rules {
#define FMPZ_POLYXX_COND_S FLINTXX_COND_S(fmpz_polyxx)
#define FMPZ_POLYXX_COND_T FLINTXX_COND_T(fmpz_polyxx)

FLINTXX_DEFINE_EQUALS(fmpz_polyxx, fmpz_poly_equal(e1._poly(), e2._poly()))

FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLYXX_COND_T, FMPZ_POLYXX_COND_S,
        fmpz_poly_set(to._poly(), from._poly()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLYXX_COND_T,
        traits::is_signed_integer,
        fmpz_poly_set_si(to._poly(), from))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLYXX_COND_T,
        traits::is_unsigned_integer,
        fmpz_poly_set_ui(to._poly(), from))
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_POLYXX_COND_T, FMPZXX_COND_S,
        fmpz_poly_set_fmpz(to._poly(), from._fmpz()))
FLINTXX_DEFINE_ASSIGN_STR(fmpz_polyxx, execution_check(
            !fmpz_poly_set_str(to._poly(), from), "assign string", "fmpz_polyxx"))

FLINTXX_DEFINE_TO_STR(fmpz_polyxx, fmpz_poly_get_str(from._poly()))
FLINTXX_DEFINE_SWAP(fmpz_polyxx, fmpz_poly_swap(e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(reverse_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_poly_reverse(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_add(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_sub(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_neg(to._poly(), from._poly()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_scalar_mul_fmpz(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_signed_integer,
        fmpz_poly_scalar_mul_si(to._poly(), e1._poly(), e2))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_scalar_mul_ui(to._poly(), e1._poly(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(mul_2exp_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_scalar_mul_2exp(to._poly(), e1._poly(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(divided_by, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_scalar_fdiv_fmpz(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_scalar_fdiv_ui(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_signed_integer,
        fmpz_poly_scalar_fdiv_si(to._poly(), e1._poly(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(fdiv_2exp_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_scalar_fdiv_2exp(to._poly(), e1._poly(), e2))

#define FMPZ_POLYXX_DEFINE_SCALAR_DIVFUNCS(name) \
FLINT_DEFINE_CBINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S, \
        fmpz_poly_scalar_##name##_fmpz(to._poly(), e1._poly(), e2._fmpz())) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer, \
        fmpz_poly_scalar_##name##_ui(to._poly(), e1._poly(), e2)) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, traits::is_signed_integer, \
        fmpz_poly_scalar_##name##_si(to._poly(), e1._poly(), e2))
FMPZ_POLYXX_DEFINE_SCALAR_DIVFUNCS(tdiv)
FMPZ_POLYXX_DEFINE_SCALAR_DIVFUNCS(divexact)

FLINT_DEFINE_CBINARY_EXPR_COND2(tdiv_2exp_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer,
        fmpz_poly_scalar_tdiv_2exp(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_scalar_mod_fmpz(to._poly(), e1._poly(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(smod_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_scalar_smod_fmpz(to._poly(), e1._poly(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_bit_pack_op, fmpzxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_mp_bitcnt_t,
        fmpz_poly_bit_pack(to._fmpz(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_bit_unpack_op, fmpz_polyxx,
        FMPZXX_COND_S, traits::fits_into_mp_bitcnt_t,
        fmpz_poly_bit_unpack(to._poly(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_bit_unpack_unsigned_op, fmpz_polyxx,
        FMPZXX_COND_S, traits::fits_into_mp_bitcnt_t,
        fmpz_poly_bit_unpack_unsigned(to._poly(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_mul(to._poly(), e1._poly(), e2._poly()))

#define FMPZ_POLYXX_DEFINE_MUL(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S, \
        fmpz_poly_##name(to._poly(), e1._poly(), e2._poly()))
FMPZ_POLYXX_DEFINE_MUL(mul_classical)
FMPZ_POLYXX_DEFINE_MUL(mulmid_classical)
FMPZ_POLYXX_DEFINE_MUL(mul_karatsuba)
FMPZ_POLYXX_DEFINE_MUL(mul_SS)
FMPZ_POLYXX_DEFINE_MUL(mul_KS)

FLINT_DEFINE_UNARY_EXPR_COND(sqr_KS_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_sqr_KS(to._poly(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_karatsuba_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_sqr_karatsuba(to._poly(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_classical_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_sqr_classical(to._poly(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_sqr(to._poly(), from._poly()))

#define FMPZ_POLYXX_DEFINE_SQRLOW(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, traits::fits_into_slong, \
        fmpz_poly_##name(to._poly(), e1._poly(), e2))
FMPZ_POLYXX_DEFINE_SQRLOW(sqrlow_KS)
FMPZ_POLYXX_DEFINE_SQRLOW(sqrlow_karatsuba_n)
FMPZ_POLYXX_DEFINE_SQRLOW(sqrlow_classical)
FMPZ_POLYXX_DEFINE_SQRLOW(sqrlow)

#define FMPZ_POLYXX_DEFINE_POW(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpz_polyxx, \
        FMPZ_POLYXX_COND_S, traits::is_unsigned_integer, \
        fmpz_poly_##name(to._poly(), e1._poly(), e2))
FMPZ_POLYXX_DEFINE_POW(pow_multinomial)
FMPZ_POLYXX_DEFINE_POW(pow_binomial)
FMPZ_POLYXX_DEFINE_POW(pow_addchains)
FMPZ_POLYXX_DEFINE_POW(pow_binexp)
FMPZ_POLYXX_DEFINE_POW(pow)

FLINT_DEFINE_BINARY_EXPR_COND2(poly_shift_left_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_poly_shift_left(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(poly_shift_right_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_poly_shift_right(to._poly(), e1._poly(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(height_op, fmpzxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_height(to._fmpz(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(poly_2norm_op, fmpzxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_2norm(to._fmpz(), from._poly()))

FMPZ_POLYXX_DEFINE_MUL(gcd)
FMPZ_POLYXX_DEFINE_MUL(gcd_subresultant)
FMPZ_POLYXX_DEFINE_MUL(gcd_heuristic)
FMPZ_POLYXX_DEFINE_MUL(gcd_modular)
FMPZ_POLYXX_DEFINE_MUL(lcm)

FLINT_DEFINE_BINARY_EXPR_COND2(resultant_op, fmpzxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_resultant(to._fmpz(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(content_op, fmpzxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_content(to._fmpz(), from._poly()))
FLINT_DEFINE_UNARY_EXPR_COND(primitive_part_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_primitive_part(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(div_basecase_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_div_basecase(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(div_divconquer_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_div_divconquer(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_div(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(rem_basecase_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_rem_basecase(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_rem(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(div_root_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_div_root(to._poly(), e1._poly(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(inv_series_newton_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_poly_inv_series_newton(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(inv_series_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_poly_inv_series(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(pseudo_rem_cohen_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
        fmpz_poly_pseudo_rem_cohen(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(derivative_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_derivative(to._poly(), from._poly()))

FMPZ_POLYXX_DEFINE_MUL(compose)
FMPZ_POLYXX_DEFINE_MUL(compose_horner)
FMPZ_POLYXX_DEFINE_MUL(compose_divconquer)

FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpzxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_evaluate_fmpz(to._fmpz(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_divconquer_op, fmpzxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_evaluate_divconquer_fmpz(to._fmpz(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_horner_op, fmpzxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_evaluate_horner_fmpz(to._fmpz(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, fmpz_vecxx,
        FMPZ_POLYXX_COND_S, FMPZ_VECXX_COND_S,
        fmpz_poly_evaluate_fmpz_vec(to._data().array, e1._poly(),
            e2._data().array, e2.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpz_polyxx_interpolate_op, fmpz_polyxx,
        FMPZ_VECXX_COND_S, FMPZ_VECXX_COND_S,
        fmpz_poly_interpolate_fmpz_vec(to._poly(), e1._data().array,
            e2._data().array, e2.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(taylor_shift_horner_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_taylor_shift_horner(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(taylor_shift_divconquer_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_taylor_shift_divconquer(to._poly(), e1._poly(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(taylor_shift_op, fmpz_polyxx,
        FMPZ_POLYXX_COND_S, FMPZXX_COND_S,
        fmpz_poly_taylor_shift(to._poly(), e1._poly(), e2._fmpz()))

FMPZ_POLYXX_DEFINE_SQRLOW(revert_series)
FMPZ_POLYXX_DEFINE_SQRLOW(revert_series_newton)
FMPZ_POLYXX_DEFINE_SQRLOW(revert_series_lagrange_fast)
FMPZ_POLYXX_DEFINE_SQRLOW(revert_series_lagrange)

FLINT_DEFINE_UNARY_EXPR_COND(sqrt_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        execution_check(fmpz_poly_sqrt(to._poly(), from._poly()),
            "sqrt", "fmpz_polyxx"))
FLINT_DEFINE_UNARY_EXPR_COND(sqrt_classical_op, fmpz_polyxx, FMPZ_POLYXX_COND_S,
        execution_check(fmpz_poly_sqrt_classical(to._poly(), from._poly()),
            "sqrt_classical", "fmpz_polyxx"))

FLINT_DEFINE_UNARY_EXPR_COND(fmpz_polyxx_product_roots_op, fmpz_polyxx,
        FMPZ_VECXX_COND_S,
        fmpz_poly_product_roots_fmpz_vec(to._poly(),
            from._data().array, from.size()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpz_polyxx_get_coeff_op, fmpzxx,
        FMPZ_POLYXX_COND_S, traits::is_integer,
        fmpz_poly_get_coeff_fmpz(to._fmpz(), e1._poly(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(fmpz_polyxx_lead_op, fmpzxx,
        FMPZ_POLYXX_COND_S,
        fmpz_set(to._fmpz(), fmpz_poly_lead(from._poly())))

FLINT_DEFINE_UNARY_EXPR_COND(poly_bound_roots_op, fmpzxx, FMPZ_POLYXX_COND_S,
        fmpz_poly_bound_roots(to._fmpz(), from._poly()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<fmpz_polyxx, fmpz_polyxx>::type>::type
    fmpz_polyxx_pair;
typedef make_ltuple<mp::make_tuple<fmpzxx, fmpz_polyxx, fmpz_polyxx>::type>::type
    fmpzxx_fmpz_polyxx_pair;
} // rdetail
#define FMPZ_POLYXX_DEFINE_DIVREM(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, rdetail::fmpz_polyxx_pair, \
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S, \
    fmpz_poly_##name(to.template get<0>()._poly(), to.template get<1>()._poly(), \
        e1._poly(), e2._poly()))
FMPZ_POLYXX_DEFINE_DIVREM(divrem_basecase)
FMPZ_POLYXX_DEFINE_DIVREM(divrem_divconquer)
FMPZ_POLYXX_DEFINE_DIVREM(divrem)
FMPZ_POLYXX_DEFINE_DIVREM(pseudo_divrem_cohen)

#define FMPZ_POLYXX_DEFINE_MULFUNC(name) \
FLINT_DEFINE_THREEARY_EXPR_COND3(name##_op, fmpz_polyxx, \
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S, traits::fits_into_slong, \
    fmpz_poly_##name(to._poly(), e1._poly(), e2._poly(), e3))

FMPZ_POLYXX_DEFINE_MULFUNC(mullow_classical)
FMPZ_POLYXX_DEFINE_MULFUNC(mulhigh_classical)
FMPZ_POLYXX_DEFINE_MULFUNC(mullow_karatsuba_n)
FMPZ_POLYXX_DEFINE_MULFUNC(mulhigh_karatsuba_n)
FMPZ_POLYXX_DEFINE_MULFUNC(mullow_KS)
FMPZ_POLYXX_DEFINE_MULFUNC(mullow_SS)
FMPZ_POLYXX_DEFINE_MULFUNC(mullow)
FMPZ_POLYXX_DEFINE_MULFUNC(mulhigh_n)

FLINT_DEFINE_BINARY_EXPR_COND2(xgcd_op, rdetail::fmpzxx_fmpz_polyxx_pair,
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
    fmpz_poly_xgcd(to.template get<0>()._fmpz(), to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(xgcd_modular_op, rdetail::fmpzxx_fmpz_polyxx_pair,
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
    fmpz_poly_xgcd_modular(to.template get<0>()._fmpz(),
        to.template get<1>()._poly(),
        to.template get<2>()._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_THREEARY_EXPR_COND3(pow_trunc_op, fmpz_polyxx,
    FMPZ_POLYXX_COND_S, traits::is_unsigned_integer, traits::fits_into_slong,
    fmpz_poly_pow_trunc(to._poly(), e1._poly(), e2, e3))

#define FMPZ_POLYXX_DEFINE_SERIES(name) \
FLINT_DEFINE_THREEARY_EXPR_COND3(name##_op, fmpz_polyxx, \
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S, traits::fits_into_slong, \
    fmpz_poly_##name(to._poly(), e1._poly(), e2._poly(), e3))
FMPZ_POLYXX_DEFINE_SERIES(div_series)
FMPZ_POLYXX_DEFINE_SERIES(compose_series_brent_kung)
FMPZ_POLYXX_DEFINE_SERIES(compose_series_horner)
FMPZ_POLYXX_DEFINE_SERIES(compose_series)

namespace rdetail {
typedef make_ltuple<mp::make_tuple<fmpz_polyxx, fmpz_polyxx, ulong>::type>::type
    fmpz_polyxx_pair_ulong;
typedef make_ltuple<mp::make_tuple<fmpz_polyxx, ulong>::type>::type
    fmpz_polyxx_ulong;
typedef make_ltuple<mp::make_tuple<bool, fmpz_polyxx>::type>::type
    bool_fmpz_polyxx;
}
#define FMPZ_POLYXX_DEFINE_PSEUDO_DIVREM(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, rdetail::fmpz_polyxx_pair_ulong, \
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S, \
    fmpz_poly_##name(to.template get<0>()._poly(), to.template get<1>()._poly(), \
        &to.template get<2>(), e1._poly(), e2._poly()))
FMPZ_POLYXX_DEFINE_PSEUDO_DIVREM(pseudo_divrem_basecase)
FMPZ_POLYXX_DEFINE_PSEUDO_DIVREM(pseudo_divrem_divconquer)
FMPZ_POLYXX_DEFINE_PSEUDO_DIVREM(pseudo_divrem)

FLINT_DEFINE_BINARY_EXPR_COND2(pseudo_div_op, rdetail::fmpz_polyxx_ulong,
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
    fmpz_poly_pseudo_div(to.template get<0>()._poly(), &to.template get<1>(),
        e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(pseudo_rem_op, rdetail::fmpz_polyxx_ulong,
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
    fmpz_poly_pseudo_rem(to.template get<0>()._poly(), &to.template get<1>(),
        e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_divides_op, rdetail::bool_fmpz_polyxx,
    FMPZ_POLYXX_COND_S, FMPZ_POLYXX_COND_S,
    to.template get<0> () =
        fmpz_poly_divides(to.template get<1>()._poly(), e1._poly(), e2._poly()))
} // rules

// immediate functions
// TODO make lazy when we have nmod class
template<class Poly>
inline typename mp::enable_all_fmpz_polyxx<mp_limb_t, Poly>::type
evaluate_mod(const Poly& p, mp_limb_t x, mp_limb_t n)
{
    return fmpz_poly_evaluate_mod(p.evaluate()._poly(), x, n);
}
} // flint

#endif
