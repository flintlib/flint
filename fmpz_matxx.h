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

#ifndef FMPZ_MATXX_H
#define FMPZ_MATXX_H FMPZ_MATXX_H

#include <algorithm> // max

#include "fmpz_mat.h"

#include "fmpzxx.h"
#include "fmpz_polyxx.h"
#include "flintxx/ltuple.h"

// TODO input and output
// TODO modular reduction
// TODO addmul
// TODO row reduction
// TODO modular gaussian elimination
// TODO echelon form
// TODO randtest static versions

namespace flint {
FLINT_DEFINE_BINOP(det_modular)
FLINT_DEFINE_BINOP(det_modular_accelerated)
FLINT_DEFINE_BINOP(mul_multi_mod)
FLINT_DEFINE_BINOP(mat_solve)
FLINT_DEFINE_BINOP(mat_solve_fflu)
FLINT_DEFINE_BINOP(mat_solve_dixon)
FLINT_DEFINE_BINOP(mat_solve_cramer)
FLINT_DEFINE_BINOP(mat_solve_bound)
FLINT_DEFINE_THREEARY(det_modular_given_divisor)
FLINT_DEFINE_THREEARY(fmpz_matxx_at)
FLINT_DEFINE_UNOP(det_bareiss)
FLINT_DEFINE_UNOP(det_bound)
FLINT_DEFINE_UNOP(det_cofactor)
FLINT_DEFINE_UNOP(det_divisor)

namespace detail {
template<class Operation>
struct fmpz_matxx_outsize;

struct fmpz_matxx_traits_generic
{
    template<class M>
    static slong rows(const M& m)
    {
        return fmpz_matxx_outsize<typename M::operation_t>::rows(m);
    }
    template<class M>
    static slong cols(const M& m)
    {
        return fmpz_matxx_outsize<typename M::operation_t>::cols(m);
    }
};

template<class Mat>
struct fmpz_matxx_traits
    : fmpz_matxx_traits_generic
{
    template<class T, class U>
    struct at
    {
        typedef FLINT_THREEARY_ENABLE_RETTYPE(fmpz_matxx_at, Mat, T, U)
            entry_ref_t;
        typedef entry_ref_t entry_srcref_t;

        static entry_srcref_t get(const Mat& m, T i, U j)
        {
            return fmpz_matxx_at(m, i, j);
        }
    };
};
} // detail

template<class Operation, class Data>
class fmpz_matxx_expression
    : public expression<derived_wrapper<fmpz_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_matxx_traits<fmpz_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(fmpz_matxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_matxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_matxx_expression, fmpz_mat_struct, _mat)

public:
    template<class T, class U>
    typename traits_t::template at<T, U>::entry_ref_t at(T i, U j)
        {return traits_t::template at<T, U>::get(*this, i, j);}
    template<class T, class U>
    typename traits_t::template at<T, U>::entry_ref_t at(T i, U j) const
        {return traits_t::template at<T, U>::get(*this, i, j);}

    slong rows() const {return traits_t::rows(*this);}
    slong cols() const {return traits_t::cols(*this);}

    evaluated_t create_temporary() const
    {
        return evaluated_t(rows(), cols());
    }

    // these only make sense with targets
    void set_randbits(frandxx& state, mp_bitcnt_t bits)
        {fmpz_mat_randbits(_mat(), state._data(), bits);}
    void set_randtest(frandxx& state, mp_bitcnt_t bits)
        {fmpz_mat_randtest(_mat(), state._data(), bits);}
    void set_randintrel(frandxx& state, mp_bitcnt_t bits)
        {fmpz_mat_randintrel(_mat(), state._data(), bits);}
    void set_randsimdioph(frandxx& state, mp_bitcnt_t bits, mp_bitcnt_t bits2)
        {fmpz_mat_randsimdioph(_mat(), state._data(), bits, bits2);}
    void set_randntrulike(frandxx& state, mp_bitcnt_t bits, ulong q)
        {fmpz_mat_randntrulike(_mat(), state._data(), bits, q);}
    void set_randntrulike2(frandxx& state, mp_bitcnt_t bits, ulong q)
        {fmpz_mat_randntrulike2(_mat(), state._data(), bits, q);}
    void set_randajtai(frandxx& state, mp_bitcnt_t bits, double alpha)
        {fmpz_mat_randajtai(_mat(), state._data(), bits, alpha);}
    void set_randrank(frandxx& state, slong rank, mp_bitcnt_t bits)
        {fmpz_mat_randrank(_mat(), state._data(), rank, bits);}

    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_randdet(frandxx& state, const Fmpz& d)
        {fmpz_mat_randdet(_mat(), state._data(), d.evaluate()._fmpz());}

    template<class Vec>
    int set_randpermdiag(frandxx& state, const Vec& v)
    {
        return fmpz_mat_randpermdiag(_mat(), state._data(), v._array(), v.size());
    }

    void apply_randops(frandxx& state, slong count)
        {fmpz_mat_randops(_mat(), state._data(), count);}

    // these cause evaluation
    slong rank() const {return fmpz_mat_rank(this->evaluate()._mat());}
    bool is_zero() const {return fmpz_mat_is_zero(this->evaluate()._mat());}
    bool is_empty() const {return fmpz_mat_is_empty(this->evaluate()._mat());}
    bool is_square() const {return fmpz_mat_is_square(this->evaluate()._mat());}

    // forwarded lazy ops
    FLINTXX_DEFINE_MEMBER_BINOP(det_modular)
    FLINTXX_DEFINE_MEMBER_BINOP(det_modular_accelerated)
    FLINTXX_DEFINE_MEMBER_BINOP(divexact)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_multi_mod)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)

    FLINTXX_DEFINE_MEMBER_BINOP_(solve, mat_solve)
    FLINTXX_DEFINE_MEMBER_BINOP_(solve_bound, mat_solve_bound)
    FLINTXX_DEFINE_MEMBER_BINOP_(solve_cramer, mat_solve_cramer)
    FLINTXX_DEFINE_MEMBER_BINOP_(solve_dixon, mat_solve_dixon)
    FLINTXX_DEFINE_MEMBER_BINOP_(solve_fflu, mat_solve_fflu)

    FLINTXX_DEFINE_MEMBER_3OP(det_modular_given_divisor)

    FLINTXX_DEFINE_MEMBER_UNOP(sqr)
    FLINTXX_DEFINE_MEMBER_UNOP(inv)
    FLINTXX_DEFINE_MEMBER_UNOP(nullspace)
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)

    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpz_polyxx, charpoly)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, det)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, det_bareiss)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, det_bound)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, det_cofactor)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, det_divisor)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(fmpzxx, trace)
};

namespace detail {
struct fmpz_mat_data;
} // detail

typedef fmpz_matxx_expression<operations::immediate, detail::fmpz_mat_data> fmpz_matxx;
typedef fmpz_matxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_matxx, fmpz_mat_struct> > fmpz_matxx_ref;
typedef fmpz_matxx_expression<operations::immediate, flint_classes::srcref_data<
    fmpz_matxx, fmpz_matxx_ref, fmpz_mat_struct> > fmpz_matxx_srcref;

namespace detail {
template<>
struct fmpz_matxx_traits<fmpz_matxx_srcref>
    : fmpz_matxx_traits_generic
{
    template<class T, class U>
    struct at
    {
        typedef fmpzxx_srcref entry_ref_t;
        typedef fmpzxx_srcref entry_srcref_t;

        template<class M>
        static fmpzxx_srcref get(const M& m, T i, U j)
        {
            return fmpzxx_srcref::make(fmpz_mat_entry(m._mat(), i, j));
        }
    };
};
template<>
struct fmpz_matxx_traits<fmpz_matxx_ref>
    : fmpz_matxx_traits<fmpz_matxx_srcref>
{
    template<class T, class U>
    struct at
        : fmpz_matxx_traits<fmpz_matxx_srcref>::template at<T, U>
    {
        typedef fmpzxx_ref entry_ref_t;

        template<class M>
        static fmpzxx_ref get(M& m, T i, U j)
        {
            return fmpzxx_ref::make(fmpz_mat_entry(m._mat(), i, j));
        }
    };
};
template<> struct fmpz_matxx_traits<fmpz_matxx>
    : fmpz_matxx_traits<fmpz_matxx_ref> { };

struct fmpz_mat_data
{
    typedef fmpz_mat_t& data_ref_t;
    typedef const fmpz_mat_t& data_srcref_t;

    fmpz_mat_t inner;

    fmpz_mat_data(slong m, slong n)
    {
        fmpz_mat_init(inner, m, n);
    }

    fmpz_mat_data(const fmpz_mat_data& o)
    {
        fmpz_mat_init_set(inner, o.inner);
    }

    ~fmpz_mat_data() {fmpz_mat_clear(inner);}
};
} // detail

#define FMPZ_MATXX_COND_S FLINTXX_COND_S(fmpz_matxx)
#define FMPZ_MATXX_COND_T FLINTXX_COND_T(fmpz_matxx)

namespace traits {
template<> struct use_temporary_merging<fmpz_matxx> : mp::false_ { };

template<class T> struct is_fmpz_matxx : mp::or_<
     traits::is_T_expr<T, fmpz_matxx>,
     flint_classes::is_source<fmpz_matxx, T> > { };
} // traits
namespace mp {
template<class T1, class T2 = void, class T3 = void, class T4 = void>
struct all_fmpz_matxx : mp::and_<all_fmpz_matxx<T1>, all_fmpz_matxx<T2, T3, T4> > { };
template<class T>
struct all_fmpz_matxx<T, void, void, void> : traits::is_fmpz_matxx<T> { };

template<class Out, class T1, class T2 = void, class T3 = void, class T4 = void>
struct enable_all_fmpz_matxx
    : mp::enable_if<all_fmpz_matxx<T1, T2, T3, T4>, Out> { };
} // mp

namespace detail {
// TODO this does not actually have anything to do with *fmpz*_matxx
template<class Operation>
struct fmpz_matxx_outsize_generic
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.cols();
    }

    template<class Data1, class Data2>
    static slong rows(const fmpz_matxx_expression<
            Operation, tuple<Data1, tuple<Data2, empty_tuple> > >& m,
            typename mp::enable_if<traits::is_fmpz_matxx<
                typename traits::basetype<Data2>::type> >::type* = 0)
    {
        return m._data().tail.head.rows();
    }
    template<class Data1, class Data2>
    static slong cols(const fmpz_matxx_expression<
            Operation, tuple<Data1, tuple<Data2, empty_tuple> > >& m,
            typename mp::enable_if<traits::is_fmpz_matxx<
                typename traits::basetype<Data2>::type> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }
};
template<class Operation>
struct fmpz_matxx_outsize : fmpz_matxx_outsize_generic<Operation> { };
template<>
struct fmpz_matxx_outsize<operations::immediate>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return fmpz_mat_nrows(m._mat());
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return fmpz_mat_ncols(m._mat());
    }
};

template<class Mat>
struct both_mat : mp::all_fmpz_matxx<
    typename traits::basetype<typename Mat::data_t::head_t>::type,
    typename traits::basetype<typename Mat::data_t::tail_t::head_t>::type> { };
template<>
struct fmpz_matxx_outsize<operations::times>
{
    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::enable_if<both_mat<Mat> >::type* = 0)
    {
        return m._data().head.rows();
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::enable_if<both_mat<Mat> >::type* = 0)
    {
        return m._data().tail.head.cols();
    }

    template<class Mat>
    static slong rows(const Mat& m,
            typename mp::disable_if<both_mat<Mat> >::type* = 0)
    {
        return fmpz_matxx_outsize_generic<operations::times>::rows(m);
    }
    template<class Mat>
    static slong cols(const Mat& m,
            typename mp::disable_if<both_mat<Mat> >::type* = 0)
    {
        return fmpz_matxx_outsize_generic<operations::times>::cols(m);
    }
};
template<>
struct fmpz_matxx_outsize<operations::mul_classical_op>
    : fmpz_matxx_outsize<operations::times> { };
template<>
struct fmpz_matxx_outsize<operations::mul_multi_mod_op>
    : fmpz_matxx_outsize<operations::times> { };

template<>
struct fmpz_matxx_outsize<operations::transpose_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.rows();
    }
};

template<>
struct fmpz_matxx_outsize<operations::nullspace_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().head.cols();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().head.cols();
    }
};

template<unsigned n>
struct fmpz_matxx_outsize<operations::ltuple_get_op<n> >
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return fmpz_matxx_outsize<
            typename Mat::data_t::head_t::operation_t>::rows(m._data().head);
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return fmpz_matxx_outsize<
            typename Mat::data_t::head_t::operation_t>::cols(m._data().head);
    }
};

// This is not actually a matrix expression, but called by the above ...
template<>
struct fmpz_matxx_outsize<operations::mat_solve_op>
{
    template<class Mat>
    static slong rows(const Mat& m)
    {
        return m._data().second().rows();
    }
    template<class Mat>
    static slong cols(const Mat& m)
    {
        return m._data().second().cols();
    }
};

template<> struct fmpz_matxx_outsize<operations::mat_solve_cramer_op>
    : fmpz_matxx_outsize<operations::mat_solve_op> { };
template<> struct fmpz_matxx_outsize<operations::mat_solve_dixon_op>
    : fmpz_matxx_outsize<operations::mat_solve_op> { };
} // detail

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MATXX_COND_T, FMPZ_MATXX_COND_S,
        fmpz_mat_set(to._mat(), from._mat()))

FLINTXX_DEFINE_SWAP(fmpz_matxx, fmpz_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(fmpz_matxx, fmpz_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZXX_COND_S,
        fmpz_mat_scalar_mul_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, traits::is_unsigned_integer,
        fmpz_mat_scalar_mul_ui(to._mat(), e1._mat(), e2))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, fmpz_matxx,
        FMPZ_MATXX_COND_S, traits::is_signed_integer,
        fmpz_mat_scalar_mul_si(to._mat(), e1._mat(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(divexact_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZXX_COND_S,
        fmpz_mat_scalar_divexact_fmpz(to._mat(), e1._mat(), e2._fmpz()))
FLINT_DEFINE_BINARY_EXPR_COND2(divexact_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, traits::is_unsigned_integer,
        fmpz_mat_scalar_divexact_ui(to._mat(), e1._mat(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(divexact_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, traits::is_signed_integer,
        fmpz_mat_scalar_divexact_si(to._mat(), e1._mat(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_sub(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_matxx, FMPZ_MATXX_COND_S,
        fmpz_mat_neg(to._mat(), from._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, fmpz_matxx, FMPZ_MATXX_COND_S,
        fmpz_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, fmpzxx, FMPZ_MATXX_COND_S,
        fmpz_mat_trace(to._fmpz(), from._mat()))

#define FMPZ_MATXX_DEFINE_DET(name) \
FLINT_DEFINE_UNARY_EXPR_COND(name##_op, fmpzxx, FMPZ_MATXX_COND_S, \
        fmpz_mat_##name(to._fmpz(), from._mat()))
FMPZ_MATXX_DEFINE_DET(det)
FMPZ_MATXX_DEFINE_DET(det_cofactor)
FMPZ_MATXX_DEFINE_DET(det_bareiss)
FMPZ_MATXX_DEFINE_DET(det_divisor)
FMPZ_MATXX_DEFINE_DET(det_bound)

namespace rdetail {
template<class T>
struct is_bool : mp::equal_types<T, bool> { };
} // rdetail
FLINT_DEFINE_BINARY_EXPR_COND2(det_modular_op, fmpzxx,
        FMPZ_MATXX_COND_S, rdetail::is_bool,
        fmpz_mat_det_modular(to._fmpz(), e1._mat(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(det_modular_accelerated_op, fmpzxx,
        FMPZ_MATXX_COND_S, rdetail::is_bool,
        fmpz_mat_det_modular_accelerated(to._fmpz(), e1._mat(), e2))
FLINT_DEFINE_THREEARY_EXPR_COND3(det_modular_given_divisor_op, fmpzxx,
        FMPZ_MATXX_COND_S, FMPZXX_COND_S, rdetail::is_bool,
        fmpz_mat_det_modular_given_divisor(to._fmpz(), e1._mat(), e2._fmpz(), e3))

FLINT_DEFINE_THREEARY_EXPR_COND3(fmpz_matxx_at_op, fmpzxx,
        FMPZ_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        fmpz_set(to._fmpz(), fmpz_mat_entry(e1._mat(), e2, e3)))

FLINT_DEFINE_BINARY_EXPR_COND2(mul_classical_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(mul_multi_mod_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_mul(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(sqr_op, fmpz_matxx, FMPZ_MATXX_COND_S,
        fmpz_mat_sqr(to._mat(), from._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, fmpz_matxx,
        FMPZ_MATXX_COND_S, traits::is_unsigned_integer,
        fmpz_mat_pow(to._mat(), e1._mat(), e2))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<bool, fmpz_matxx, fmpzxx>::type >::type
    fmpz_mat_inv_rt;

typedef make_ltuple<mp::make_tuple<fmpzxx, fmpzxx>::type >::type
    fmpz_mat_solve_bound_rt;

typedef make_ltuple<mp::make_tuple<slong, fmpz_matxx>::type >::type
    fmpz_mat_nullspace_rt;
} // rdetail

FLINT_DEFINE_UNARY_EXPR_COND(inv_op, rdetail::fmpz_mat_inv_rt,
        FMPZ_MATXX_COND_S,
        to.template get<0>() = fmpz_mat_inv(to.template get<1>()._mat(),
            to.template get<2>()._fmpz(), from._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(charpoly_op, fmpz_polyxx, FMPZ_MATXX_COND_S,
        fmpz_mat_charpoly(to._poly(), from._mat()))

#define FMPZ_MATXX_DEFINE_SOLVE(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, rdetail::fmpz_mat_inv_rt, \
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S, \
        to.template get<0>() = fmpz_##name(to.template get<1>()._mat(), \
            to.template get<2>()._fmpz(), e1._mat(), e2._mat()))
FMPZ_MATXX_DEFINE_SOLVE(mat_solve)
FMPZ_MATXX_DEFINE_SOLVE(mat_solve_dixon)
FMPZ_MATXX_DEFINE_SOLVE(mat_solve_cramer)
FMPZ_MATXX_DEFINE_SOLVE(mat_solve_fflu)

FLINT_DEFINE_BINARY_EXPR_COND2(mat_solve_bound_op,
        rdetail::fmpz_mat_solve_bound_rt,
        FMPZ_MATXX_COND_S, FMPZ_MATXX_COND_S,
        fmpz_mat_solve_bound(to.template get<0>()._fmpz(),
            to.template get<1>()._fmpz(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(nullspace_op, rdetail::fmpz_mat_nullspace_rt,
        FMPZ_MATXX_COND_S, to.template get<0>() = fmpz_mat_nullspace(
            to.template get<1>()._mat(), from._mat()))

// temporary instantiation stuff

template<class Expr>
struct use_default_temporary_instantiation<Expr, fmpz_matxx> : mp::false_ { };
template<class Expr>
struct instantiate_temporaries<Expr, fmpz_matxx,
    typename mp::disable_if<traits::is_fmpz_matxx<Expr> >::type>
{
    // The only case where this should ever happen is if Expr is an ltuple
    // TODO static assert this
    static fmpz_matxx get(const Expr& e)
    {
        typedef typename Expr::operation_t op_t;
        return fmpz_matxx(
                detail::fmpz_matxx_outsize<op_t>::rows(e),
                detail::fmpz_matxx_outsize<op_t>::cols(e));
    }
};
template<class Expr>
struct instantiate_temporaries<Expr, fmpz_matxx,
    typename mp::enable_if<traits::is_fmpz_matxx<Expr> >::type>
{
    static fmpz_matxx get(const Expr& e)
    {
        return e.create_temporary();
    }
};
} // rules

// immediate functions
template<class Mat>
inline typename mp::enable_if<traits::is_fmpz_matxx<Mat>, slong>::type
rank(const Mat& m)
{
    return m.rank();
}
} // flint

#endif
