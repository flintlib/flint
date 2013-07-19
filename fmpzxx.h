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

#ifndef CXX_FMPZXX_H
#define CXX_FMPZXX_H

#include <cstdlib>

#include "flintxx/evaluation_tools.h"
#include "flintxx/expression.h"
#include "flintxx/expression_traits.h"
#include "flintxx/flint_classes.h"
#include "flintxx/frandxx.h"
#include "flintxx/stdmath.h"

#include "fmpz.h"

// TODO swap

namespace flint {

template<class Operation, class Data>
class fmpzxx_expression
    : public expression<derived_wrapper<fmpzxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpzxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(fmpzxx_expression)
    FLINTXX_DEFINE_CTORS(fmpzxx_expression)
    FLINTXX_DEFINE_C_REF(fmpzxx_expression, fmpz, _fmpz)

    // these only make sense with fmpzxx
    FLINTXX_DEFINE_RANDFUNC(fmpz, randbits)
    FLINTXX_DEFINE_RANDFUNC(fmpz, randtest)
    FLINTXX_DEFINE_RANDFUNC(fmpz, randtest_unsigned)
    FLINTXX_DEFINE_RANDFUNC(fmpz, randtest_not_zero)

    template<class Fmpz>
    static fmpzxx_expression randm(frandxx& state, const Fmpz& m)
    {
        fmpzxx_expression res;
        fmpz_randm(res._fmpz(), state._data(), m.evaluate()._fmpz());
        return res;
    }
    template<class Fmpz>
    static fmpzxx_expression randtest_mod(frandxx& state, const Fmpz& m)
    {
        fmpzxx_expression res;
        fmpz_randtest_mod(res._fmpz(), state._data(), m.evaluate()._fmpz());
        return res;
    }
    template<class Fmpz>
    static fmpzxx_expression randtest_mod_signed(frandxx& state, const Fmpz& m)
    {
        fmpzxx_expression res;
        fmpz_randtest_mod_signed(res._fmpz(), state._data(), m.evaluate()._fmpz());
        return res;
    }

    // TODO would these make more sense static?
    void set_ui_smod(mp_limb_t x, mp_limb_t m)
    {
        fmpz_set_ui_smod(this->_fmpz(), x, m);
    }
    void set_uiui(mp_limb_t hi, mp_limb_t lo)
    {
        fmpz_set_uiui(this->_fmpz(), hi, lo);
    }
    void neg_uiui(mp_limb_t hi, mp_limb_t lo)
    {
        fmpz_neg_uiui(this->_fmpz(), hi, lo);
    }

    // These make sense with all expressions, but cause evaluation
    double get_d_2exp(long& exp) const
    {
        return fmpz_get_d_2exp(&exp, this->evaluate()._fmpz());
    }
    bool is_zero() const
    {
        return fmpz_is_zero(this->evaluate()._fmpz());
    }
    bool is_one() const
    {
        return fmpz_is_one(this->evaluate()._fmpz());
    }
    bool is_pm1() const
    {
        return fmpz_is_pm1(this->evaluate()._fmpz());
    }
    bool is_even() const
    {
        return fmpz_is_even(this->evaluate()._fmpz());
    }
    bool is_odd() const
    {
        return fmpz_is_odd(this->evaluate()._fmpz());
    }
    bool is_square() const
    {
        return fmpz_is_square(this->evaluate()._fmpz());
    }
};

namespace detail {
struct fmpz_data;
}

typedef fmpzxx_expression<operations::immediate, detail::fmpz_data> fmpzxx;
typedef fmpzxx_expression<operations::immediate,
            flint_classes::ref_data<fmpzxx, fmpz> > fmpzxx_ref;
typedef fmpzxx_expression<operations::immediate,
            flint_classes::srcref_data<fmpzxx, fmpzxx_ref, fmpz> > fmpzxx_srcref;

namespace detail {
struct fmpz_data
{
    typedef fmpz_t& data_ref_t;
    typedef const fmpz_t& data_srcref_t;

    fmpz_t inner;

    fmpz_data() {fmpz_init(inner);}
    ~fmpz_data() {fmpz_clear(inner);}
    fmpz_data(const fmpz_data& o) {fmpz_init_set(inner, o.inner);}

    fmpz_data(const char* str)
    {
        fmpz_init(inner);
        fmpz_set_str(inner, str, 10);
    }

    template<class T>
    fmpz_data(const T& t)
    {
        init(t);
    }

    template<class T>
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type init(T t)
    {
        fmpz_init_set_ui(inner, t);
    }
    template<class T>
    typename mp::enable_if<traits::is_signed_integer<T> >::type init(T t)
    {
        fmpz_init(inner);
        fmpz_set_si(inner, t);
    }

    void init(const fmpzxx_srcref& r)
    {
        fmpz_init_set(inner, r._fmpz());
    }
};

} // detail

///////////////////////////////////////////////////////////////////
// HELPERS
///////////////////////////////////////////////////////////////////
namespace traits {
template<class T> struct is_fmpzxx : mp::or_<
     traits::is_T_expr<T, fmpzxx>,
     flint_classes::is_source<fmpzxx, T> > { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpzxx : mp::and_<all_fmpzxx<T1>, all_fmpzxx<T2, T3, T4> > { };
template<class T>
struct all_fmpzxx<T, void, void, void> : traits::is_fmpzxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpzxx
    : mp::enable_if<all_fmpzxx<T1, T2, T3, T4>, Out> { };
} // mp

///////////////////////////////////////////////////////////////////
// RULES
///////////////////////////////////////////////////////////////////
namespace rules {

#define FMPZXX_COND_S FLINTXX_COND_S(fmpzxx)
#define FMPZXX_COND_T FLINTXX_COND_T(fmpzxx)

FLINT_DEFINE_DOIT_COND2(assignment, FMPZXX_COND_T, FMPZXX_COND_S,
        fmpz_set(to._fmpz(), from._fmpz()))

FLINT_DEFINE_DOIT_COND2(assignment,
        FMPZXX_COND_T, traits::is_unsigned_integer,
        fmpz_set_ui(to._fmpz(), from))

FLINT_DEFINE_DOIT_COND2(assignment,
        FMPZXX_COND_T, traits::is_signed_integer,
        fmpz_set_si(to._fmpz(), from))

template<class T, int n>
struct assignment<T, char[n],
    typename mp::enable_if<FMPZXX_COND_S<T> >::type>
{
    static void doit(T& target, const char* source)
    {
        fmpz_set_str(target._fmpz(), const_cast<char*>(source), 10);
    }
};

FLINTXX_DEFINE_CMP(fmpzxx, fmpz_cmp(e1._fmpz(), e2._fmpz()))

template<class T, class U>
struct cmp<T, U,
    typename mp::enable_if<mp::and_<
        FMPZXX_COND_S<T>, traits::is_signed_integer<U> > >::type>
{
    static int get(const T& v, const U& t)
    {
        return fmpz_cmp_si(v._fmpz(), t);
    }
};

template<class T>
struct cmp<fmpzxx, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static int get(const fmpzxx& v, const T& t)
    {
        return fmpz_cmp_ui(v._fmpz(), t);
    }
};

FLINTXX_DEFINE_TO_STR(fmpzxx, fmpz_get_str(0,  base, from._fmpz()))

FLINT_DEFINE_GET_COND(conversion, slong, FMPZXX_COND_S,
        fmpz_get_si(from._fmpz()))
FLINT_DEFINE_GET_COND(conversion, ulong, FMPZXX_COND_S,
        fmpz_get_ui(from._fmpz()))
FLINT_DEFINE_GET_COND(conversion, double, FMPZXX_COND_S,
        fmpz_get_d(from._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_add(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_CBINARY_EXPR_COND2(plus, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_add_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_mul(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_mul_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpzxx,
        FMPZXX_COND_S, traits::is_signed_integer,
        fmpz_mul_si(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_sub(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_sub_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_fdiv_q(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_fdiv_q_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpzxx,
        FMPZXX_COND_S, traits::is_signed_integer,
        fmpz_fdiv_q_si(to._fmpz(), e1._fmpz(), e2))

// TODO this interpretation of mod is not the same as for builtin types!
FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_mod(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(modulo, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_mod_ui(to._fmpz(), e1._fmpz(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpzxx, FMPZXX_COND_S,
        fmpz_neg(to._fmpz(), from._fmpz()))

namespace rdetail {
template<class Fmpz1, class Fmpz2, class T>
void fmpzxx_shift(Fmpz1& to, const Fmpz2& from, T howmuch)
{
    if(howmuch < 0)
        fmpz_fdiv_q_2exp(to._fmpz(), from._fmpz(), -howmuch);
    else
        fmpz_mul_2exp(to._fmpz(), from._fmpz(), howmuch);
}
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(shift, fmpzxx,
        FMPZXX_COND_S, traits::is_integer,
        rdetail::fmpzxx_shift(to, e1, e2))
} // rules

FLINTXX_DEFINE_TERNARY(fmpzxx,
        fmpz_addmul(to._fmpz(), e1._fmpz(), e2._fmpz()),
        fmpz_submul(to._fmpz(), e1._fmpz(), e2._fmpz()),
        FLINTXX_UNADORNED_MAKETYPES)

///////////////////////////////////////////////////////////////////////////
// FUNCTIONS
///////////////////////////////////////////////////////////////////////////

// These functions evaluate immediately, and (often) do not yield fmpzxxs

template<class T1, class T2>
inline typename mp::enable_all_fmpzxx<bool, T1, T2>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible(t1.evaluate()._fmpz(), t2.evaluate()._fmpz());
}
template<class T1, class T2>
inline typename mp::enable_if<mp::and_<
    traits::is_fmpzxx<T1>, traits::fits_into_slong<T2> >, bool>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible_si(t1.evaluate()._fmpz(), t2);
}

template<class Fmpz1, class Fmpz2>
inline typename mp::enable_all_fmpzxx<long, Fmpz1, Fmpz2>::type
clog(const Fmpz1& x, const Fmpz2& b)
{
    return fmpz_clog(x.evaluate()._fmpz(), b.evaluate()._fmpz());
}
template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, long>::type
clog(const Fmpz& x, ulong b)
{
    return fmpz_clog_ui(x.evaluate()._fmpz(), b);
}

template<class Fmpz1, class Fmpz2>
inline typename mp::enable_all_fmpzxx<long, Fmpz1, Fmpz2>::type
flog(const Fmpz1& x, const Fmpz2& b)
{
    return fmpz_flog(x.evaluate()._fmpz(), b.evaluate()._fmpz());
}
template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, long>::type
flog(const Fmpz& x, ulong b)
{
    return fmpz_flog_ui(x.evaluate()._fmpz(), b);
}

template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, double>::type
dlog(const Fmpz& x)
{
    return fmpz_dlog(x.evaluate()._fmpz());
}

// TODO These cannot yet be made lazy since we have no code for ternary ops...
// In any case, the uiui suffix should probably be dropped.
template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, fmpzxx>::type
mul2_uiui(const Fmpz& g, ulong x, ulong y)
{
    fmpzxx res;
    fmpz_mul2_uiui(res._fmpz(), g.evaluate()._fmpz(), x, y);
    return res;
}
template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, fmpzxx>::type
divexact2_uiui(const Fmpz& g, ulong x, ulong y)
{
    fmpzxx res;
    fmpz_divexact2_uiui(res._fmpz(), g.evaluate()._fmpz(), x, y);
    return res;
}
template<class Fmpz1, class Fmpz2>
inline typename mp::enable_all_fmpzxx<fmpzxx, Fmpz1, Fmpz2>::type
powm(const Fmpz1& g, ulong e, const Fmpz2& m)
{
    fmpzxx res;
    fmpz_powm_ui(res._fmpz(), g.evaluate()._fmpz(), e, m.evaluate()._fmpz());
    return res;
}
template<class Fmpz1, class Fmpz2, class Fmpz3>
inline typename mp::enable_all_fmpzxx<fmpzxx, Fmpz1, Fmpz2, Fmpz3>::type
powm(const Fmpz1& g, const Fmpz2 e, const Fmpz3& m)
{
    fmpzxx res;
    fmpz_powm(res._fmpz(), g.evaluate()._fmpz(), e.evaluate()._fmpz(),
            m.evaluate()._fmpz());
    return res;
}
template<class Fmpz1, class Fmpz2>
inline typename mp::enable_all_fmpzxx<fmpzxx, Fmpz1, Fmpz2>::type
mul_tdiv_q_2exp(const Fmpz1& g, const Fmpz2& x, ulong exp)
{
    fmpzxx res;
    fmpz_mul_tdiv_q_2exp(res._fmpz(), g.evaluate()._fmpz(),
            x.evaluate()._fmpz(), exp);
    return res;
}
template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, fmpzxx>::type
mul_tdiv_q_2exp(const Fmpz& g, long x, ulong exp)
{
    fmpzxx res;
    fmpz_mul_si_tdiv_q_2exp(res._fmpz(), g.evaluate()._fmpz(), x, exp);
    return res;
}
// TODO addmul, submul?

// These cannot be lazy because we do not support two return values
template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>, FMPZXX_COND_T<Fmpz2>,
        traits::is_fmpzxx<Fmpz3>, traits::is_fmpzxx<Fmpz4>
    > >::type
fdiv_qr(Fmpz1& f, Fmpz2& s, const Fmpz3& g, const Fmpz4& h)
{
    fmpz_fdiv_qr(f._fmpz(), s._fmpz(), g.evaluate()._fmpz(), h.evaluate()._fmpz());
}
template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>, FMPZXX_COND_T<Fmpz2>,
        traits::is_fmpzxx<Fmpz3>, traits::is_fmpzxx<Fmpz4>
    > >::type
tdiv_qr(Fmpz1& f, Fmpz2& s, const Fmpz3& g, const Fmpz4& h)
{
    fmpz_tdiv_qr(f._fmpz(), s._fmpz(), g.evaluate()._fmpz(), h.evaluate()._fmpz());
}

template<class Fmpz1, class Fmpz2, class Fmpz3>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>,
        traits::is_fmpzxx<Fmpz2>, traits::is_fmpzxx<Fmpz3>
    >, bool>::type
sqrtmod(Fmpz1& b, const Fmpz2& a, const Fmpz3& p)
{
    return fmpz_sqrtmod(b._fmpz(), a.evaluate()._fmpz(), p.evaluate()._fmpz());
}
template<class Fmpz1, class Fmpz2, class Fmpz3>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>, FMPZXX_COND_T<Fmpz2>,
        traits::is_fmpzxx<Fmpz3>
    > >::type
sqrtrem(Fmpz1& f, Fmpz2& s, const Fmpz3& g)
{
    fmpz_sqrtrem(f._fmpz(), s._fmpz(), g.evaluate()._fmpz());
}

template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>, FMPZXX_COND_T<Fmpz2>,
        traits::is_fmpzxx<Fmpz3>, traits::is_fmpzxx<Fmpz4>
    > >::type
gcdinv(Fmpz1& d, Fmpz2& a, const Fmpz3& f, const Fmpz4& g)
{
    fmpz_gcdinv(d._fmpz(), a._fmpz(), f.evaluate()._fmpz(), g.evaluate()._fmpz());
}
template<class Fmpz1, class Fmpz2, class Fmpz3, class Fmpz4, class Fmpz5>
inline typename mp::enable_if<mp::and_<
        FMPZXX_COND_T<Fmpz1>, FMPZXX_COND_T<Fmpz2>, FMPZXX_COND_T<Fmpz3>,
        traits::is_fmpzxx<Fmpz4>, traits::is_fmpzxx<Fmpz5>
    > >::type
gcdinv(Fmpz1& d, Fmpz2& a, Fmpz3& b, const Fmpz4& f, const Fmpz5& g)
{
    fmpz_xgcd(d._fmpz(), a._fmpz(), b._fmpz(),
            f.evaluate()._fmpz(), g.evaluate()._fmpz());
}

// These functions are evaluated lazily

FLINT_DEFINE_BINOP(cdiv_q)
FLINT_DEFINE_BINOP(fdiv_r)
FLINT_DEFINE_BINOP(tdiv_q)
FLINT_DEFINE_BINOP(fdiv_r_2exp)
FLINT_DEFINE_BINOP(tdiv_q_2exp)
FLINT_DEFINE_BINOP(divexact)
FLINT_DEFINE_UNOP(fac)
FLINT_DEFINE_UNOP(fib)
FLINT_DEFINE_BINOP(rfac)
FLINT_DEFINE_BINOP(bin)
FLINT_DEFINE_BINOP(gcd)
FLINT_DEFINE_BINOP(lcm)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(rfac_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_rfac_ui(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(fac_op, fmpzxx, traits::is_unsigned_integer,
        fmpz_fac_ui(to._fmpz(), from))
FLINT_DEFINE_UNARY_EXPR_COND(fib_op, fmpzxx, traits::is_unsigned_integer,
        fmpz_fib_ui(to._fmpz(), from))
FLINT_DEFINE_BINARY_EXPR_COND2(gcd_op, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_gcd(to._fmpz(), e1._fmpz(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(lcm_op, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_lcm(to._fmpz(), e1._fmpz(), e2._fmpz()))

template<class T1, class T2>
struct binary_expression<
    T1,
    typename mp::enable_if<
        mp::and_<
            traits::is_unsigned_integer<T1>,
            traits::is_unsigned_integer<T2> >,
        operations::bin_op>::type,
    T2>
{
    typedef fmpzxx return_t;
    template<class V>
    static void doit(V& to, const T1& t1, const T2& t2)
    {
        fmpz_bin_uiui(to._fmpz(), t1, t2);
    }
};

#define FMPZXX_DEFINE_DIVFUNCS(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S, \
        fmpz_##name(to._fmpz(), e1._fmpz(), e2._fmpz())) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpzxx, FMPZXX_COND_S, \
        traits::is_signed_integer, \
        fmpz_##name##_si(to._fmpz(), e1._fmpz(), e2)) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, fmpzxx, FMPZXX_COND_S, \
        traits::is_unsigned_integer, \
        fmpz_##name##_ui(to._fmpz(), e1._fmpz(), e2))

FMPZXX_DEFINE_DIVFUNCS(cdiv_q)
FMPZXX_DEFINE_DIVFUNCS(tdiv_q)
FMPZXX_DEFINE_DIVFUNCS(divexact)
FLINT_DEFINE_BINARY_EXPR_COND2(fdiv_r_op, fmpzxx, FMPZXX_COND_S, FMPZXX_COND_S,
        fmpz_fdiv_r(to._fmpz(), e1._fmpz(), e2._fmpz()))

FLINT_DEFINE_BINARY_EXPR_COND2(tdiv_q_2exp_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_tdiv_q_2exp(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(fdiv_r_2exp_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_fdiv_r_2exp(to._fmpz(), e1._fmpz(), e2))

// standard math functions (c/f stdmath.h)
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpzxx,
        FMPZXX_COND_S, traits::is_unsigned_integer,
        fmpz_pow_ui(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(root_op, fmpzxx,
        FMPZXX_COND_S, traits::fits_into_slong,
        fmpz_root(to._fmpz(), e1._fmpz(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(sqrt_op, fmpzxx, FMPZXX_COND_S,
        fmpz_sqrt(to._fmpz(), from._fmpz()))
FLINT_DEFINE_UNARY_EXPR_COND(abs_op, fmpzxx, FMPZXX_COND_S,
        fmpz_abs(to._fmpz(), from._fmpz()))
} // rules

} // flint

#endif
