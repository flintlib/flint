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

#ifndef FMPQ_MATXX_H
#define FMPQ_MATXX_H FMPQ_MATXX_H

#include "fmpq_mat.h"

#include "fmpqxx.h"
#include "fmpz_matxx.h"
#include "fmpz_vecxx.h"

#include "flintxx/ltuple.h"
#include "flintxx/matrix.h"

// TODO wrap entry_num, entry_den?
// TODO numden_rowwise_2
// TODO rref members

namespace flint {
FLINT_DEFINE_BINOP(hilbert_matrix)
FLINT_DEFINE_BINOP(mul_direct)
FLINT_DEFINE_BINOP(mul_cleared)
FLINT_DEFINE_BINOP(solve_fraction_free)

FLINT_DEFINE_UNOP(numden_colwise)
FLINT_DEFINE_UNOP(numden_entrywise)
FLINT_DEFINE_UNOP(numden_matwise)
FLINT_DEFINE_UNOP(numden_rowwise)
FLINT_DEFINE_UNOP(num_rowwise)
FLINT_DEFINE_UNOP(num_colwise)

FLINT_DEFINE_UNOP(rref_classical)
FLINT_DEFINE_UNOP(rref_fraction_free)

namespace detail {
template<class Mat>
struct fmpq_matxx_traits : matrices::generic_traits<Mat> { };

typedef make_ltuple<mp::make_tuple<fmpz_matxx, fmpz_matxx>::type>::type
    fmpq_matxx_numden_entrywise_rt;
typedef make_ltuple<mp::make_tuple<fmpz_matxx, fmpzxx>::type>::type
    fmpq_matxx_numden_matwise_rt;
typedef make_ltuple<mp::make_tuple<fmpz_matxx, fmpz_vecxx>::type>::type
    fmpq_matxx_numden_rowwise_rt;
typedef fmpq_matxx_numden_rowwise_rt fmpq_matxx_numden_colwise_rt;
} // detail

template<class Operation, class Data>
class fmpq_matxx_expression
    : public expression<derived_wrapper<fmpq_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpq_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpq_matxx_traits<fmpq_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(fmpq_matxx_expression)
    FLINTXX_DEFINE_CTORS(fmpq_matxx_expression)
    FLINTXX_DEFINE_C_REF(fmpq_matxx_expression, fmpq_mat_struct, _mat)

    template<class Expr>
    static evaluated_t create_temporary_rowscols(
            const Expr&, slong rows, slong cols)
    {
        return evaluated_t(rows, cols);
    }
    FLINTXX_DEFINE_MATRIX_METHODS(traits_t)

    template<class Fmpz_mat, class Fmpz>
    static fmpq_matxx_expression reconstruct(const Fmpz_mat& mat, const Fmpz& mod,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpq_matxx_expression res(mat.rows(), mat.cols());
        res.set_reconstruct(mat, mod);
        return res;
    }
    template<class Fmpz_mat, class Fmpz>
    static fmpq_matxx_expression frac(const Fmpz_mat& num, const Fmpz& den,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpq_matxx_expression res(num.rows(), num.cols());
        res.set_frac(num, den);
        return res;
    }
    template<class Fmpz_mat>
    static fmpq_matxx_expression integer_matrix(const Fmpz_mat& mat,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0)
    {
        fmpq_matxx_expression res(mat.rows(), mat.cols());
        res = mat;
        return res;
    }

    static fmpq_matxx_expression randbits(slong rows, slong cols,
            frandxx& state, mp_bitcnt_t bits)
    {
        fmpq_matxx_expression res(rows, cols);
        res.set_randbits(state, bits);
        return res;
    }
    static fmpq_matxx_expression randtest(slong rows, slong cols,
            frandxx& state, mp_bitcnt_t bits)
    {
        fmpq_matxx_expression res(rows, cols);
        res.set_randtest(state, bits);
        return res;
    }

    static fmpq_matxx_expression zero(slong rows, slong cols)
        {return fmpq_matxx_expression(rows, cols);}
    static fmpq_matxx_expression one(slong rows, slong cols)
    {
        fmpq_matxx_expression res(rows, cols);
        res.set_one();
        return res;
    }

    // these only make sense with targets
    void set_randbits(frandxx& state, mp_bitcnt_t bits)
        {fmpq_mat_randbits(_mat(), state._data(), bits);}
    void set_randtest(frandxx& state, mp_bitcnt_t bits)
        {fmpq_mat_randtest(_mat(), state._data(), bits);}
    void set_hilbert_matrix()
        {fmpq_mat_hilbert_matrix(_mat());}
    void set_zero()
        {fmpq_mat_zero(_mat());}
    void set_one()
        {fmpq_mat_one(_mat());}

    template<class Fmpz_mat, class Fmpz>
    void set_frac(const Fmpz_mat& num, const Fmpz& den,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpq_mat_set_fmpz_mat_div_fmpz(num.evaluate()._mat(),
                den.evaluate()._fmpz());
    }
    template<class Fmpz_mat, class Fmpz>
    void set_reconstruct(const Fmpz_mat& mat, const Fmpz& mod,
            typename mp::enable_if<traits::is_fmpz_matxx<Fmpz_mat> >::type* = 0,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        execution_check(fmpq_mat_set_fmpz_mat_mod_fmpz(
                    _mat(), mat.evaluate()._mat(), mod.evaluate()._fmpz()),
                "reconstruct", "fmpq_matxx");
    }

    bool pivot(slong r, slong c, permxx* perm = 0)
        {return fmpq_mat_pivot(maybe_perm_data(perm), _mat(), r, c);}

    // these cause evaluation
    slong rank() const {return fmpq_mat_rank(this->evaluate()._mat());}
    bool is_zero() const {return fmpq_mat_is_zero(this->evaluate()._mat());}
    bool is_empty() const {return fmpq_mat_is_empty(this->evaluate()._mat());}
    bool is_square() const {return fmpq_mat_is_square(this->evaluate()._mat());}
    bool is_integral() const
        {return fmpq_mat_is_integral(this->evaluate()._mat());}

    // forwarded lazy ops
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpqxx, det)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpqxx, trace)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(detail::fmpq_matxx_numden_entrywise_rt,
            numden_entrywise)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(detail::fmpq_matxx_numden_matwise_rt,
            numden_matwise)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(detail::fmpq_matxx_numden_rowwise_rt,
            numden_rowwise)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(detail::fmpq_matxx_numden_colwise_rt,
            numden_colwise)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(fmpz_matxx, num_rowwise,
            num_rowwise)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE_(fmpz_matxx, num_colwise,
            num_colwise)

    FLINTXX_DEFINE_MEMBER_BINOP(mul_cleared)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_direct)
    FLINTXX_DEFINE_MEMBER_BINOP(solve_dixon)
    FLINTXX_DEFINE_MEMBER_BINOP(solve_fraction_free)
};

namespace detail {
struct fmpq_mat_data;
} // detail

typedef fmpq_matxx_expression<operations::immediate, detail::fmpq_mat_data> fmpq_matxx;
typedef fmpq_matxx_expression<operations::immediate,
            flint_classes::ref_data<fmpq_matxx, fmpq_mat_struct> > fmpq_matxx_ref;
typedef fmpq_matxx_expression<operations::immediate, flint_classes::srcref_data<
    fmpq_matxx, fmpq_matxx_ref, fmpq_mat_struct> > fmpq_matxx_srcref;

template<>
struct matrix_traits<fmpq_matxx>
{
    template<class M> static slong rows(const M& m)
    {
        return fmpq_mat_nrows(m._mat());
    }
    template<class M> static slong cols(const M& m)
    {
        return fmpq_mat_ncols(m._mat());
    }

    template<class M> static fmpqxx_srcref at(const M& m, slong i, slong j)
    {
        return fmpqxx_srcref::make(fmpq_mat_entry(m._mat(), i, j));
    }
    template<class M> static fmpqxx_ref at(M& m, slong i, slong j)
    {
        return fmpqxx_ref::make(fmpq_mat_entry(m._mat(), i, j));
    }
};

namespace detail {
template<>
struct fmpq_matxx_traits<fmpq_matxx_srcref>
    : matrices::generic_traits_srcref<fmpqxx_srcref> { };
template<>
struct fmpq_matxx_traits<fmpq_matxx_ref>
    : matrices::generic_traits_ref<fmpqxx_ref> { };
template<> struct fmpq_matxx_traits<fmpq_matxx>
    : matrices::generic_traits_nonref<fmpqxx_ref, fmpqxx_srcref> { };

struct fmpq_mat_data
{
    typedef fmpq_mat_t& data_ref_t;
    typedef const fmpq_mat_t& data_srcref_t;

    fmpq_mat_t inner;

    fmpq_mat_data(slong m, slong n)
    {
        fmpq_mat_init(inner, m, n);
    }

    fmpq_mat_data(const fmpq_mat_data& o)
    {
        fmpq_mat_init(inner, fmpq_mat_nrows(o.inner), fmpq_mat_ncols(o.inner));
        fmpq_mat_set(inner, o.inner);
    }

    fmpq_mat_data(fmpq_matxx_srcref o)
    {
        fmpq_mat_init(inner, o.rows(), o.cols());
        fmpq_mat_set(inner, o._data().inner);
    }

    ~fmpq_mat_data() {fmpq_mat_clear(inner);}
};
} // detail

#define FMPQ_MATXX_COND_S FLINTXX_COND_S(fmpq_matxx)
#define FMPQ_MATXX_COND_T FLINTXX_COND_T(fmpq_matxx)

namespace traits {
template<class T> struct is_fmpq_matxx
    : flint_classes::is_Base<fmpq_matxx, T> { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpq_matxx : mp::and_<all_fmpq_matxx<T1>, all_fmpq_matxx<T2, T3, T4> > { };
template<class T>
struct all_fmpq_matxx<T, void, void, void> : traits::is_fmpq_matxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpq_matxx
    : mp::enable_if<all_fmpq_matxx<T1, T2, T3, T4>, Out> { };
} // mp

namespace matrices {
template<>
struct outsize<operations::mul_direct_op>
    : outsize<operations::times> { };
template<>
struct outsize<operations::mul_cleared_op>
    : outsize<operations::times> { };

template<> struct outsize<operations::solve_fraction_free_op>
    : outsize<operations::solve_op> { };

template<>
struct outsize<operations::hilbert_matrix_op>
{
    template<class Expr>
    static slong rows(const Expr& e) {return e._data().first();}
    template<class Expr>
    static slong cols(const Expr& e) {return e._data().second();}
};

template<>
struct outsize<operations::numden_entrywise_op>
{
    template<class Expr>
    static slong rows(const Expr& e) {return e._data().first().rows();}
    template<class Expr>
    static slong cols(const Expr& e) {return e._data().first().cols();}
};
template<> struct outsize<operations::numden_colwise_op>
    : outsize<operations::numden_entrywise_op> { };
template<> struct outsize<operations::numden_rowwise_op>
    : outsize<operations::numden_entrywise_op> { };
template<> struct outsize<operations::numden_matwise_op>
    : outsize<operations::numden_entrywise_op> { };
template<> struct outsize<operations::num_rowwise_op>
    : outsize<operations::numden_entrywise_op> { };
template<> struct outsize<operations::num_colwise_op>
    : outsize<operations::numden_entrywise_op> { };
} // matrices

namespace vectors {
template<>
struct outsize<operations::numden_rowwise_op>
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        return e._data().first().rows();
    }
};
template<>
struct outsize<operations::numden_colwise_op>
{
    template<class Expr>
    static unsigned get(const Expr& e)
    {
        return e._data().first().cols();
    }
};
} // vectors

// temporary instantiation stuff
FLINTXX_DEFINE_TEMPORARY_RULES(fmpq_matxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_MATXX_COND_T, FMPQ_MATXX_COND_S,
        fmpq_mat_set(to._mat(), from._mat()))
FLINT_DEFINE_DOIT_COND2(assignment, FMPQ_MATXX_COND_T, FMPZ_MATXX_COND_S,
        fmpq_mat_set_fmpz_mat(to._mat(), from._mat()))

FLINTXX_DEFINE_SWAP(fmpq_matxx, fmpq_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(fmpq_matxx, fmpq_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_PRINT_COND(FMPQ_MATXX_COND_S, (fmpq_mat_print(from._mat()), 1))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpq_matxx,
        FMPZ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_mul_r_fmpz_mat(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpq_mat_mul_fmpz_mat(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPZXX_COND_S,
        fmpq_mat_scalar_mul_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(divided_by, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPZXX_COND_S,
        fmpq_mat_scalar_div_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(mul_direct_op, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_mul_direct(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(mul_cleared_op, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_mul_cleared(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        fmpq_mat_sub(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpq_matxx, FMPQ_MATXX_COND_S,
        fmpq_mat_neg(to._mat(), from._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, fmpq_matxx, FMPQ_MATXX_COND_S,
        fmpq_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, fmpqxx, FMPQ_MATXX_COND_S,
        fmpq_mat_trace(to._fmpq(), from._mat()))

FLINT_DEFINE_THREEARY_EXPR_COND3(mat_at_op, fmpqxx,
        FMPQ_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        fmpq_set(to._fmpq(), fmpq_mat_entry(e1._mat(), e2, e3)))

FLINT_DEFINE_BINARY_EXPR_COND2(hilbert_matrix_op, fmpq_matxx,
        traits::fits_into_slong, traits::fits_into_slong,
        to.set_hilbert_matrix())

FLINT_DEFINE_UNARY_EXPR_COND(det_op, fmpqxx, FMPQ_MATXX_COND_S,
        fmpq_mat_det(to._fmpq(), from._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(solve_fraction_free_op, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        execution_check(fmpq_mat_solve_fraction_free(
                to._mat(), e1._mat(), e2._mat()),
            "solve", "fmpq_mat"))
FLINT_DEFINE_BINARY_EXPR_COND2(solve_dixon_op, fmpq_matxx,
        FMPQ_MATXX_COND_S, FMPQ_MATXX_COND_S,
        execution_check(fmpq_mat_solve_dixon(to._mat(), e1._mat(), e2._mat()),
            "solve", "fmpq_mat"))
FLINT_DEFINE_UNARY_EXPR_COND(inv_op, fmpq_matxx, FMPQ_MATXX_COND_S,
        execution_check(fmpq_mat_inv(to._mat(), from._mat()),
            "inv", "fmpq_mat"))

FLINT_DEFINE_UNARY_EXPR_COND(num_rowwise_op, fmpz_matxx,
        FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_rowwise(to._mat(), 0, from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(num_colwise_op, fmpz_matxx,
        FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_colwise(to._mat(), 0, from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(numden_entrywise_op,
        detail::fmpq_matxx_numden_entrywise_rt, FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_entrywise(to.template get<0>()._mat(),
            to.template get<1>()._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(numden_matwise_op,
        detail::fmpq_matxx_numden_matwise_rt, FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_matwise(to.template get<0>()._mat(),
            to.template get<1>()._fmpz(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(numden_rowwise_op,
        detail::fmpq_matxx_numden_rowwise_rt, FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_rowwise(to.template get<0>()._mat(),
            to.template get<1>()._array(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(numden_colwise_op,
        detail::fmpq_matxx_numden_colwise_rt, FMPQ_MATXX_COND_S,
        fmpq_mat_get_fmpz_mat_colwise(to.template get<0>()._mat(),
            to.template get<1>()._array(), from._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<slong, fmpq_matxx>::type>::type
    fmpq_matxx_rref_rt;
}

FLINT_DEFINE_UNARY_EXPR_COND(rref_op, rdetail::fmpq_matxx_rref_rt,
        FMPQ_MATXX_COND_S,
        to.template get<0>() =
             fmpq_mat_rref(to.template get<1>()._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(rref_classical_op, rdetail::fmpq_matxx_rref_rt,
        FMPQ_MATXX_COND_S,
        to.template get<0>() =
             fmpq_mat_rref_classical(to.template get<1>()._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(rref_fraction_free_op, rdetail::fmpq_matxx_rref_rt,
        FMPQ_MATXX_COND_S,
        to.template get<0>() =
             fmpq_mat_rref_fraction_free(to.template get<1>()._mat(), from._mat()))
} // rules
} // flint

#endif

