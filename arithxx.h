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

#ifndef ARITHXX_H
#define ARITHXX_H

#include "arith.h"

#include "fmpq_polyxx.h"
#include "fmpqxx.h"
#include "fmpz_matxx.h"
#include "fmpz_vecxx.h"
#include "fmpzxx.h"
#include "nmod_vecxx.h"

// TODO namespace arith?
// TODO arith_hrr_expsum_factored
// TODO codegen / vector improvements

namespace flint {
namespace detail {
template<class T>
inline typename T::wrapped_traits::data_srcref_t extract_data(const T& t,
        typename mp::enable_if<flint_classes::is_flint_class<T> >::type* = 0)
{
    return t._data().inner;
}
template<class T>
inline typename T::wrapped_traits::data_ref_t extract_data(T& t,
        typename mp::enable_if<flint_classes::is_flint_class<T> >::type* = 0)
{
    return t._data().inner;
}
template<class T>
const T& extract_data(const T& t,
        typename mp::disable_if<flint_classes::is_flint_class<T> >::type* = 0)
{
    return t;
}
template<class T>
typename T::arrayref_t extract_data(T& t,
        typename mp::disable_if<flint_classes::is_flint_class<T> >::type* = 0,
        typename mp::enable_if<mp::or_<rules::FMPZ_VECXX_COND_S<T>,
            rules::FMPQ_VECXX_COND_S<T> > >::type* = 0)
{
    return t._array();
}
}
#define ARITHXX_DEFINE_UNOP(name, Return, Cond) \
FLINT_DEFINE_UNOP(name) \
namespace rules { \
FLINT_DEFINE_UNARY_EXPR_COND(name##_op, Return, Cond, \
        arith_##name(detail::extract_data(to), detail::extract_data(from))) \
}
#define ARITHXX_DEFINE_BINOP(name, Return, Cond1, Cond2) \
FLINT_DEFINE_BINOP(name) \
namespace rules { \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, Return, Cond1, Cond2, \
        arith_##name(detail::extract_data(to), detail::extract_data(e1), \
            detail::extract_data(e2))) \
}

namespace at {
template<class T> struct slong : traits::fits_into_slong<T> { };
template<class T> struct ulong : traits::is_unsigned_integer<T> { };
} // at

ARITHXX_DEFINE_UNOP(primorial, fmpzxx, at::slong)
ARITHXX_DEFINE_UNOP(harmonic_number, fmpqxx, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_1u, fmpzxx, at::slong, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_1, fmpzxx, at::slong, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_2, fmpzxx, at::slong, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_1u_vec, fmpz_vecxx, at::slong, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_1_vec, fmpz_vecxx, at::slong, at::slong)
ARITHXX_DEFINE_BINOP(stirling_number_2_vec, fmpz_vecxx, at::slong, at::slong)
FLINT_DEFINE_BINOP(stirling_number_1u_vec_next)
FLINT_DEFINE_BINOP(stirling_number_1_vec_next)
FLINT_DEFINE_BINOP(stirling_number_2_vec_next)
FLINT_DEFINE_BINOP(stirling_matrix_1u)
FLINT_DEFINE_BINOP(stirling_matrix_1)
FLINT_DEFINE_BINOP(stirling_matrix_2)

namespace vectors {
template<>
struct outsize<operations::stirling_number_1u_vec_op>
{
    template<class Expr>
    static unsigned get(const Expr& e) {return e._data().second();}
};
template<> struct outsize<operations::stirling_number_1_vec_op>
    : outsize<operations::stirling_number_1u_vec_op> { };
template<> struct outsize<operations::stirling_number_2_vec_op>
    : outsize<operations::stirling_number_1u_vec_op> { };

template<>
struct outsize<operations::stirling_number_1u_vec_next_op>
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        slong r = e._data().first().size();
        if(r == e._data().second())
            return r + 1;
        return r;
    }
};
template<> struct outsize<operations::stirling_number_1_vec_next_op>
    : outsize<operations::stirling_number_1u_vec_next_op> { };
template<> struct outsize<operations::stirling_number_2_vec_next_op>
    : outsize<operations::stirling_number_1u_vec_next_op> { };
} // vectors

namespace matrices {
template<>
struct outsize<operations::stirling_matrix_1u_op>
{
    template<class Expr>
    static slong rows(const Expr& e) {return e._data().first();}
    template<class Expr>
    static slong cols(const Expr& e) {return e._data().second();}
};
template<> struct outsize<operations::stirling_matrix_1_op>
    : outsize<operations::stirling_matrix_1u_op> { };
template<> struct outsize<operations::stirling_matrix_2_op>
    : outsize<operations::stirling_matrix_1u_op> { };
} // matrices

namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(stirling_number_1u_vec_next_op, fmpz_vecxx,
        FMPZ_VECXX_COND_S, at::slong,
        arith_stirling_number_1u_vec_next(to._array(), e1._array(), e2,
            e1.size() + (e1.size() == e2)))
FLINT_DEFINE_BINARY_EXPR_COND2(stirling_number_1_vec_next_op, fmpz_vecxx,
        FMPZ_VECXX_COND_S, at::slong,
        arith_stirling_number_1_vec_next(to._array(), e1._array(), e2,
            e1.size() + (e1.size() == e2)))
FLINT_DEFINE_BINARY_EXPR_COND2(stirling_number_2_vec_next_op, fmpz_vecxx,
        FMPZ_VECXX_COND_S, at::slong,
        arith_stirling_number_2_vec_next(to._array(), e1._array(), e2,
            e1.size() + (e1.size() == e2)))

FLINT_DEFINE_BINARY_EXPR_COND2(stirling_matrix_1u_op, fmpz_matxx,
        at::slong, at::slong, arith_stirling_matrix_1u(to._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(stirling_matrix_1_op, fmpz_matxx,
        at::slong, at::slong, arith_stirling_matrix_1(to._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(stirling_matrix_2_op, fmpz_matxx,
        at::slong, at::slong, arith_stirling_matrix_2(to._mat()))
} // rules

ARITHXX_DEFINE_UNOP(bell_number, fmpzxx, at::ulong)
ARITHXX_DEFINE_UNOP(bell_number_bsplit, fmpzxx, at::ulong)
ARITHXX_DEFINE_UNOP(bell_number_multi_mod, fmpzxx, at::ulong)
ARITHXX_DEFINE_UNOP(bell_number_vec, fmpz_vecxx, at::slong)
ARITHXX_DEFINE_UNOP(bell_number_vec_recursive, fmpz_vecxx, at::slong)
ARITHXX_DEFINE_UNOP(bell_number_vec_multi_mod, fmpz_vecxx, at::slong)
FLINT_DEFINE_BINOP(bell_number_nmod)
FLINT_DEFINE_BINOP(bell_number_nmod_vec)
FLINT_DEFINE_BINOP(bell_number_nmod_vec_recursive)
FLINT_DEFINE_BINOP(bell_number_nmod_vec_series)

namespace vectors {
template<>
struct outsize<operations::bell_number_vec_op>
{
    template<class Expr> static unsigned get(const Expr& e)
        {return e._data().first();}
};
template<> struct outsize<operations::bell_number_vec_recursive_op>
    : outsize<operations::bell_number_vec_op> { };
template<> struct outsize<operations::bell_number_vec_multi_mod_op>
    : outsize<operations::bell_number_vec_op> { };
template<> struct outsize<operations::bell_number_nmod_vec_op>
    : outsize<operations::bell_number_vec_op> { };
template<> struct outsize<operations::bell_number_nmod_vec_recursive_op>
    : outsize<operations::bell_number_vec_op> { };
template<> struct outsize<operations::bell_number_nmod_vec_series_op>
    : outsize<operations::bell_number_vec_op> { };
} // vectors

namespace at {
template<class T> struct nmod : mp::or_<mp::equal_types<T, nmodxx_ctx_srcref>,
    mp::equal_types<T, nmodxx_ctx> > { };
} // at

#define NMODXX_DEFINE_FIND_CTX_BY_OP(exprname, eval) \
namespace traits { \
template<class Data> \
struct has_nmodxx_ctx< exprname, Data> > : mp::true_ { }; \
} \
namespace detail { \
template<class Data> \
struct get_nmodxx_ctx<exprname, Data> > \
{ \
    template<class T> \
    static nmodxx_ctx_srcref get(const T& e) {return eval;} \
}; \
}
#define FLINTXX_COMMA() ,
NMODXX_DEFINE_FIND_CTX_BY_OP(nmodxx_expression<operations::bell_number_nmod_op,
        e._data().second())
NMODXX_DEFINE_FIND_CTX_BY_OP(vector_expression<detail::nmod_vector_traits
        FLINTXX_COMMA() operations::bell_number_nmod_vec_op, e._data().second())
NMODXX_DEFINE_FIND_CTX_BY_OP(vector_expression<detail::nmod_vector_traits
        FLINTXX_COMMA() operations::bell_number_nmod_vec_series_op,
        e._data().second())
NMODXX_DEFINE_FIND_CTX_BY_OP(vector_expression<detail::nmod_vector_traits
        FLINTXX_COMMA() operations::bell_number_nmod_vec_recursive_op,
        e._data().second())
namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(bell_number_nmod_op, nmodxx, at::ulong, at::nmod,
        to.set_nored(arith_bell_number_nmod(e1, e2._nmod())))
FLINT_DEFINE_BINARY_EXPR_COND2(bell_number_nmod_vec_op, nmod_vecxx,
        at::ulong, at::nmod,
        arith_bell_number_nmod_vec(to._array(), e1, e2._nmod()))
FLINT_DEFINE_BINARY_EXPR_COND2(bell_number_nmod_vec_recursive_op, nmod_vecxx,
        at::ulong, at::nmod,
        arith_bell_number_nmod_vec_recursive(to._array(), e1, e2._nmod()))
FLINT_DEFINE_BINARY_EXPR_COND2(bell_number_nmod_vec_series_op, nmod_vecxx,
        at::ulong, at::nmod,
        arith_bell_number_nmod_vec_series(to._array(), e1, e2._nmod()))
} // rules

ARITHXX_DEFINE_UNOP(bernoulli_number, fmpqxx, at::ulong)
ARITHXX_DEFINE_UNOP(bernoulli_polynomial, fmpq_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(bernoulli_number_vec, fmpq_vecxx, at::ulong)
ARITHXX_DEFINE_UNOP(bernoulli_number_denom, fmpzxx, at::ulong)

namespace vectors {
template<>
struct outsize<operations::bernoulli_number_vec_op>
{
    template<class Expr> static unsigned get(const Expr& e)
        {return e._data().first();}
};
} // vectors

ARITHXX_DEFINE_UNOP(euler_number, fmpzxx, at::ulong)
ARITHXX_DEFINE_UNOP(euler_number_vec, fmpz_vecxx, at::ulong)
ARITHXX_DEFINE_UNOP(euler_polynomial, fmpq_polyxx, at::ulong)

namespace vectors {
template<> struct outsize<operations::euler_number_vec_op>
    : outsize<operations::bernoulli_number_vec_op> { };
}

ARITHXX_DEFINE_UNOP(legendre_polynomial, fmpq_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(chebyshev_t_polynomial, fmpz_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(chebyshev_u_polynomial, fmpz_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(euler_phi, fmpzxx, FMPZXX_COND_S)
ARITHXX_DEFINE_BINOP(divisor_sigma, fmpzxx, FMPZXX_COND_S, at::ulong)
ARITHXX_DEFINE_UNOP(divisors, fmpz_polyxx, FMPZXX_COND_S)
ARITHXX_DEFINE_UNOP(ramanujan_tau, fmpzxx, FMPZXX_COND_S)
ARITHXX_DEFINE_UNOP(ramanujan_tau_series, fmpz_polyxx, at::slong)
ARITHXX_DEFINE_UNOP(cyclotomic_polynomial, fmpz_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(cos_minpoly, fmpz_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(swinnerton_dyer_polynomial, fmpz_polyxx, at::ulong)
ARITHXX_DEFINE_UNOP(landau_function_vec, fmpz_vecxx, at::slong)

namespace vectors {
template<> struct outsize<operations::landau_function_vec_op>
    : outsize<operations::bernoulli_number_vec_op> { };
}

ARITHXX_DEFINE_BINOP(dedekind_sum_naive, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S)
ARITHXX_DEFINE_BINOP(dedekind_sum_coprime_large, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S)
ARITHXX_DEFINE_BINOP(dedekind_sum_coprime, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S)
ARITHXX_DEFINE_BINOP(dedekind_sum, fmpqxx,
        FMPZXX_COND_S, FMPZXX_COND_S)

ARITHXX_DEFINE_UNOP(number_of_partitions_vec, fmpz_vecxx, at::slong)
ARITHXX_DEFINE_UNOP(number_of_partitions, fmpzxx, at::ulong)
FLINT_DEFINE_BINOP(number_of_partitions_nmod_vec)

namespace vectors {
template<> struct outsize<operations::number_of_partitions_vec_op>
    : outsize<operations::bernoulli_number_vec_op> { };
template<> struct outsize<operations::number_of_partitions_nmod_vec_op>
    : outsize<operations::bell_number_vec_op> { };
}

NMODXX_DEFINE_FIND_CTX_BY_OP(vector_expression<detail::nmod_vector_traits
        FLINTXX_COMMA() operations::number_of_partitions_nmod_vec_op,
        e._data().second())

namespace rules {
FLINT_DEFINE_BINARY_EXPR_COND2(number_of_partitions_nmod_vec_op, nmod_vecxx,
        at::ulong, at::nmod,
        arith_number_of_partitions_nmod_vec(to._array(), e1, e2._nmod()))
}

ARITHXX_DEFINE_BINOP(sum_of_squares, fmpzxx, at::ulong, FMPZXX_COND_S)
ARITHXX_DEFINE_BINOP(sum_of_squares_vec, fmpz_vecxx, at::ulong, at::slong)

namespace vectors {
template<> struct outsize<operations::sum_of_squares_vec_op>
    : outsize<operations::stirling_number_1u_vec_op> { };
}

inline double bell_number_size(ulong n) {return arith_bell_number_size(n);}
inline double bernoulli_number_size(ulong n)
    {return arith_bernoulli_number_size(n);}
inline double euler_number_size(ulong n) {return arith_euler_number_size(n);}
inline double dedekind_sum_coprime_d(double h, double k)
    {return arith_dedekind_sum_coprime_d(h, k);}

template<class Fmpz>
inline typename mp::enable_if<traits::is_fmpzxx<Fmpz>, int>::type moebius_mu(
        const Fmpz& n)
{
    return arith_moebius_mu(n.evaluate()._fmpz());
}
} // flint

#endif
