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

#ifndef CXX_PROTOTYPE_H
#define CXX_PROTOTYPE_H

#include <cstdlib>

// TODO
// * non-explicit ctos
// * contexts

#include "cxx/evaluation_tools.h"
#include "cxx/expression.h"
#include "cxx/expression_traits.h"
#include "cxx/stdmath.h"

#include "fmpz.h"

namespace flint {
// "concrete" expression classes

template<class Operation, class Data>
class mpz_expression
    : public expression<derived_wrapper<mpz_expression>, Operation, Data>
{
public:
    mpz_expression() {}
    template<class T>
    explicit mpz_expression(const T& t)
        : expression<derived_wrapper< ::flint::mpz_expression>,
              Operation, Data>(t) {}

    template<class T>
    mpz_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit mpz_expression(const Data& d)
        : expression<derived_wrapper< ::flint::mpz_expression>,
              Operation, Data>(d) {}

    template<class D, class O, class Da>
    friend class expression;
};
typedef mpz_expression<operations::immediate, fmpz_t> mpz;

///////////////////////////////////////////////////////////////////
// HELPERS
///////////////////////////////////////////////////////////////////
namespace traits {
template<class T> struct is_mpz : is_T_expr<T, mpz> { };
} // traits
namespace mp {
template<class Out, class T1, class T2 = void>
struct enable_all_mpz
    : mp::enable_if<mp::and_<traits::is_mpz<T1>, traits::is_mpz<T2> >, Out> { };
template<class Out, class T>
struct enable_all_mpz<Out, T, void>
    : mp::enable_if<traits::is_mpz<T>, Out> { };
} // mp

///////////////////////////////////////////////////////////////////
// RULES
///////////////////////////////////////////////////////////////////
namespace rules {
template<>
struct empty_initialization<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_init(v._data());
    }
};

template<>
struct destruction<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_clear(v._data());
    }
};

FLINT_DEFINE_DOIT(initialization, mpz, mpz,
        fmpz_init_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_unsigned_integer<T>,
        fmpz_init_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(initialization, mpz, traits::is_signed_integer<T>,
        fmpz_init(to._data());fmpz_set_si(to._data(), from))

// TODO should expression automatically deduce this from the assignment
// implementation (or the other way round)?
template<int n>
struct initialization<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_init(target._data());
        // TODO what about different bases
        fmpz_set_str(target._data(), source, 10);
    }
};

FLINT_DEFINE_DOIT(assignment, mpz, mpz, fmpz_set(to._data(), from._data()))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_unsigned_integer<T>,
        fmpz_set_ui(to._data(), from))

FLINT_DEFINE_DOIT_COND(assignment, mpz, traits::is_signed_integer<T>,
        fmpz_set_si(to._data(), from))

template<int n>
struct assignment<mpz, char[n]>
{
    static void doit(mpz& target, const char* source)
    {
        fmpz_set_str(target._data(), const_cast<char*>(source), 10);
    }
};

template<>
struct cmp<mpz, mpz>
{
    static int get(const mpz& l, const mpz& r)
    {
        return fmpz_cmp(l._data(), r._data());
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_signed_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_si(v._data(), t);
    }
};

template<class T>
struct cmp<mpz, T,
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type>
{
    static int get(const mpz& v, const T& t)
    {
        return fmpz_cmp_ui(v._data(), t);
    }
};

template<>
struct to_string<mpz>
{
    static std::string get(const mpz& v, int base)
    {
        char* str = fmpz_get_str(0, base, v._data());
        std::string res(str);
        std::free(str);
        return res;
    }
};

FLINT_DEFINE_GET(conversion, slong, mpz, fmpz_get_si(from._data()))
FLINT_DEFINE_GET(conversion, ulong, mpz, fmpz_get_ui(from._data()))
FLINT_DEFINE_GET(conversion, double, mpz, fmpz_get_d(from._data()))

FLINT_DEFINE_BINARY_EXPR(plus, mpz,
        fmpz_add(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(plus, mpz, traits::is_unsigned_integer<T>,
        fmpz_add_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(times, mpz,
        fmpz_mul(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_unsigned_integer<T>,
        fmpz_mul_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND(times, mpz, traits::is_signed_integer<T>,
        fmpz_mul_si(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(minus, mpz,
        fmpz_sub(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(minus, mpz, traits::is_unsigned_integer<T>,
        fmpz_sub_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR(divided_by, mpz,
        fmpz_fdiv_q(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_unsigned_integer<T>,
        fmpz_fdiv_q_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_BINARY_EXPR_COND(divided_by, mpz, traits::is_signed_integer<T>,
        fmpz_fdiv_q_si(to._data(), e1._data(), e2))

// TODO this interpretation of mod is not the same as for builtin types!
FLINT_DEFINE_BINARY_EXPR(modulo, mpz,
        fmpz_mod(to._data(), e1._data(), e2._data()))

FLINT_DEFINE_BINARY_EXPR_COND(modulo, mpz, traits::is_unsigned_integer<T>,
        fmpz_mod_ui(to._data(), e1._data(), e2))

FLINT_DEFINE_UNARY_EXPR(negate, mpz, fmpz_neg(to._data(), from._data()))


// Optimized evaluation rules using ternary arithmetic (addmul, submul)
// a +- b*c
template<class Op, class Left, class Right1, class Right2>
struct evaluation<Op,
    tuple<Left, tuple<
        mpz_expression<operations::times,
            tuple<Right1, tuple<Right2, empty_tuple> > >,
        // NB: there is no particular reason to have the enable_if here,
        //     many other similar places would do
        typename mp::enable_if<mp::or_<
                mp::equal_types<Op, operations::plus>,
                mp::equal_types<Op, operations::minus> >,
            empty_tuple>::type> >,
    true, 1,
    typename tools::ternary_helper<mpz, Left, Right1, Right2>::enable::type>
{
    // Helpful for testing.
    static const unsigned TERNARY_OP_MARKER = 0;

    typedef mpz return_t;
    typedef tools::ternary_helper<mpz, Left, Right1, Right2> th;
    typedef typename th::temporaries_t temporaries_t;
    typedef tuple<Left, tuple<
        mpz_expression<operations::times,
            tuple<Right1, tuple<Right2, empty_tuple> > >,
        empty_tuple> > data_t;
    static const bool is_add = mp::equal_types<Op, operations::plus>::val;

    static void doit(const data_t& input, temporaries_t temps, return_t* res)
    {
        const mpz* left = 0;
        const mpz* right = 0;
        th::doit(input.first(), input.second()._data().first(),
                input.second()._data().second(), temps, res, right, left);
        if(is_add)
            fmpz_addmul(res->_data(), left->_data(), right->_data());
        else
            fmpz_submul(res->_data(), left->_data(), right->_data());
    }
};

// b*c + a
template<class Right, class Left1, class Left2>
struct evaluation<operations::plus,
    tuple<mpz_expression<operations::times,
            tuple<Left1, tuple<Left2, empty_tuple> > >,
        tuple<Right, empty_tuple> >,
    true, 1,
    typename tools::ternary_helper<mpz, 
        Right, Left1, Left2, operations::times>::enable::type>
{
    // Helpful for testing.
    static const unsigned TERNARY_OP_MARKER = 0;

    typedef mpz return_t;
    typedef tools::ternary_helper<mpz, Right, Left1, Left2> th;
    typedef typename th::temporaries_t temporaries_t;
    typedef tuple<mpz_expression<operations::times,
            tuple<Left1, tuple<Left2, empty_tuple> > >,
        tuple<Right, empty_tuple> > data_t;

    static void doit(const data_t& input, temporaries_t temps, return_t* res)
    {
        const mpz* left = 0;
        const mpz* right = 0;
        th::doit(input.second(), input.first()._data().first(),
                input.first()._data().second(), temps, res, right, left);
        fmpz_addmul(res->_data(), left->_data(), right->_data());
    }
};
} // rules

// Assignment-arithmetic using ternaries.
namespace detail {
template<class Right1, class Right2>
struct ternary_assign_helper
{
    typedef tools::evaluate_2<Right1, Right2> ev2_t;
    typedef typename ev2_t::temporaries_t temporaries_t;
    typedef mp::back_tuple<temporaries_t> back_t;

    typename back_t::type backing;
    ev2_t ev2;

    static temporaries_t backtemps(typename back_t::type& backing)
    {
        temporaries_t temps;
        back_t::init(temps, backing);
        return temps;
    }

    ternary_assign_helper(typename ev2_t::arg1_t r1, typename ev2_t::arg2_t r2)
        : ev2(backtemps(backing), r1, r2) {}
    const mpz& getleft() {return ev2.get1();}
    const mpz& getright() {return ev2.get2();}
};

template<class Right1, class Right2>
struct enable_ternary_assign
    : mp::enable_all_mpz<mpz&,
        typename traits::basetype<Right1>::type,
        typename traits::basetype<Right2>::type> { };
} // detail

// a += b*c
template<class Right1, class Right2>
inline typename detail::enable_ternary_assign<Right1, Right2>::type
operator+=(mpz& left, const mpz_expression<operations::times,
        tuple<Right1, tuple<Right2, empty_tuple> > >& other)
{
    detail::ternary_assign_helper<Right1, Right2> tah(
            other._data().first(), other._data().second());
    fmpz_addmul(left._data(), tah.getleft()._data(),
            tah.getright()._data());
    return left;
}

// a -= b*c
template<class Right1, class Right2>
inline typename detail::enable_ternary_assign<Right1, Right2>::type
operator-=(mpz& left, const mpz_expression<operations::times,
        tuple<Right1, tuple<Right2, empty_tuple> > >& other)
{
    detail::ternary_assign_helper<Right1, Right2> tah(
            other._data().first(), other._data().second());
    fmpz_submul(left._data(), tah.getleft()._data(),
            tah.getright()._data());
    return left;
}


///////////////////////////////////////////////////////////////////////////
// FUNCTIONS
///////////////////////////////////////////////////////////////////////////

// These functions evaluate immediately and do not yield mpzs

template<class T1, class T2>
inline typename mp::enable_all_mpz<bool, T1, T2>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible(t1.evaluate()._data(), t2.evaluate()._data());
}
template<class T1, class T2>
inline typename mp::enable_if<mp::and_<
    traits::is_mpz<T1>, traits::fits_into_slong<T2> >, bool>::type
divisible(const T1& t1, const T2& t2)
{
    return fmpz_divisible_si(t1.evaluate()._data(), t2);
}

// These functions are evaluated lazily

FLINT_DEFINE_UNOP(fac)
FLINT_DEFINE_BINOP(rfac)
FLINT_DEFINE_BINOP(bin)
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND(rfac_op, mpz, traits::is_unsigned_integer<T>,
        fmpz_rfac_ui(to._data(), e1._data(), e2))
FLINT_DEFINE_UNARY_EXPR_COND(fac_op, mpz, traits::is_unsigned_integer<T>,
        fmpz_fac_ui(to._data(), from))

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
    typedef mpz return_t;
    static void doit(mpz& to, const T1& t1, const T2& t2)
    {
        fmpz_bin_uiui(to._data(), t1, t2);
    }
};

// standard math functions (c/f stdmath.h)
FLINT_DEFINE_BINARY_EXPR_COND(pow_op, mpz, traits::is_unsigned_integer<T>,
        fmpz_pow_ui(to._data(), e1._data(), e2))
FLINT_DEFINE_BINARY_EXPR_COND(root_op, mpz, traits::fits_into_slong<T>,
        fmpz_root(to._data(), e1._data(), e2))
FLINT_DEFINE_UNARY_EXPR(sqrt_op, mpz, fmpz_sqrt(to._data(), from._data()))
} // rules

// TODO many more functions

} // flint

#endif
