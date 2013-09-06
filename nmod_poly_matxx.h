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

#ifndef NMOD_POLY_MATXX_H
#define NMOD_POLY_MATXX_H

#include "nmod_poly_mat.h"

#include "nmod_matxx.h"
#include "nmod_polyxx.h"
#include "permxx.h"

#include "flintxx/matrix.h"
#include "flintxx/stdmath.h"

// NOTE: it is *not* valid to use empty nmod_poly_matxx matrices!
// TODO nullspace member

namespace flint {
FLINT_DEFINE_UNOP(sqr_interpolate)
FLINT_DEFINE_BINOP(mul_interpolate)

namespace detail {
template<class Mat>
struct nmod_poly_matxx_traits : matrices::generic_traits<Mat> { };
} // detail

template<class Operation, class Data>
class nmod_poly_matxx_expression
    : public expression<derived_wrapper<nmod_poly_matxx_expression>, Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::nmod_poly_matxx_expression>,
              Operation, Data> base_t;
    typedef detail::nmod_poly_matxx_traits<nmod_poly_matxx_expression> traits_t;

    FLINTXX_DEFINE_BASICS(nmod_poly_matxx_expression)
    FLINTXX_DEFINE_CTORS(nmod_poly_matxx_expression)
    FLINTXX_DEFINE_C_REF(nmod_poly_matxx_expression, nmod_poly_mat_struct, _mat)

    // These only make sense with immediates
    nmodxx_ctx_srcref _ctx() const
    {
        return nmodxx_ctx_srcref::make(nmod_poly_mat_entry(_mat(), 0, 0)->mod);
    }

    // These work on any expression without evaluation
    nmodxx_ctx_srcref estimate_ctx() const
    {
        return tools::find_nmodxx_ctx(*this);
    }
    mp_limb_t modulus() const {return estimate_ctx().n();}

    template<class Expr>
    static evaluated_t create_temporary_rowscols(
            const Expr& e, slong rows, slong cols)
    {
        return evaluated_t(rows, cols, tools::find_nmodxx_ctx(e).n());
    }
    FLINTXX_DEFINE_MATRIX_METHODS(traits_t)

    FLINTXX_DEFINE_FORWARD_STATIC(from_ground)

    static nmod_poly_matxx_expression randtest(slong rows, slong cols,
            mp_limb_t M, frandxx& state, slong len)
    {
        nmod_poly_matxx_expression res(rows, cols, M);
        res.set_randtest(state, len);
        return res;
    }
    static nmod_poly_matxx_expression randtest_sparse(slong rows, slong cols,
            mp_limb_t M, frandxx& state, slong len, float density)
    {
        nmod_poly_matxx_expression res(rows, cols, M);
        res.set_randtest_sparse(state, len, density);
        return res;
    }

    static nmod_poly_matxx_expression zero(slong rows, slong cols, mp_limb_t n)
        {return nmod_poly_matxx_expression(rows, cols, n);}
    static nmod_poly_matxx_expression one(slong rows, slong cols, mp_limb_t n)
    {
        nmod_poly_matxx_expression res(rows, cols, n);
        res.set_one();
        return res;
    }

    // these only make sense with targets
    void set_randtest(frandxx& state, slong len)
        {nmod_poly_mat_randtest(_mat(), state._data(), len);}
    void set_randtest_sparse(frandxx& state, slong len, float density)
        {nmod_poly_mat_randtest_sparse(_mat(), state._data(), len, density);}
    void set_zero() {nmod_poly_mat_zero(_mat());}
    void set_one() {nmod_poly_mat_one(_mat());}

    // these cause evaluation
    bool is_zero() const
        {return nmod_poly_mat_is_zero(this->evaluate()._mat());}
    bool is_one() const
        {return nmod_poly_mat_is_one(this->evaluate()._mat());}
    bool is_square() const
        {return nmod_poly_mat_is_square(this->evaluate()._mat());}
    bool is_empty() const
        {return nmod_poly_mat_is_empty(this->evaluate()._mat());}
    slong max_length() const
        {return nmod_poly_mat_max_length(this->evaluate()._mat());}
    slong rank() const {return nmod_poly_mat_rank(this->evaluate()._mat());}
    slong find_pivot_any(slong start, slong end, slong c) const
    {
        return nmod_poly_mat_find_pivot_any(
                this->evaluate()._mat(), start, end, c);
    }
    slong find_pivot_partial(slong start, slong end, slong c) const
    {
        return nmod_poly_mat_find_pivot_partial(
                this->evaluate()._mat(), start, end, c);
    }

    // lazy members
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmod_polyxx, det)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmod_polyxx, det_fflu)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmod_polyxx, det_interpolate)
    FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(nmod_polyxx, trace)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_classical)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_interpolate)
    FLINTXX_DEFINE_MEMBER_UNOP(sqr_KS)
    FLINTXX_DEFINE_MEMBER_UNOP(transpose)
    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(???, nullspace) // TODO
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)
    FLINTXX_DEFINE_MEMBER_BINOP(solve)
    FLINTXX_DEFINE_MEMBER_BINOP(solve_fflu)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_classical)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_interpolate)
    FLINTXX_DEFINE_MEMBER_BINOP(mul_KS)
    FLINTXX_DEFINE_MEMBER_BINOP(pow)

    //FLINTXX_DEFINE_MEMBER_UNOP_RTYPE(???, rref) // TODO

    FLINTXX_DEFINE_MEMBER_FFLU
};

namespace detail {
struct nmod_poly_mat_data;
} // detail

typedef nmod_poly_matxx_expression<
    operations::immediate, detail::nmod_poly_mat_data> nmod_poly_matxx;
typedef nmod_poly_matxx_expression<operations::immediate,
            flint_classes::ref_data<
                nmod_poly_matxx, nmod_poly_mat_struct> > nmod_poly_matxx_ref;
typedef nmod_poly_matxx_expression<operations::immediate,
        flint_classes::srcref_data<
            nmod_poly_matxx, nmod_poly_matxx_ref,
            nmod_poly_mat_struct> > nmod_poly_matxx_srcref;

template<>
struct matrix_traits<nmod_poly_matxx>
{
    template<class M> static slong rows(const M& m)
    {
        return nmod_poly_mat_nrows(m._mat());
    }
    template<class M> static slong cols(const M& m)
    {
        return nmod_poly_mat_ncols(m._mat());
    }

    template<class M> static nmod_polyxx_srcref at(const M& m, slong i, slong j)
    {
        return nmod_polyxx_srcref::make(nmod_poly_mat_entry(m._mat(), i, j));
    }
    template<class M> static nmod_polyxx_ref at(M& m, slong i, slong j)
    {
        return nmod_polyxx_ref::make(nmod_poly_mat_entry(m._mat(), i, j));
    }
};

namespace traits {
template<> struct has_nmodxx_ctx<nmod_poly_matxx> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_poly_matxx_ref> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_poly_matxx_srcref> : mp::true_ { };
} // traits

namespace detail {
template<>
struct nmod_poly_matxx_traits<nmod_poly_matxx_srcref>
    : matrices::generic_traits_srcref<nmod_polyxx_srcref> { };
template<>
struct nmod_poly_matxx_traits<nmod_poly_matxx_ref>
    : matrices::generic_traits_ref<nmod_polyxx_ref> { };
template<> struct nmod_poly_matxx_traits<nmod_poly_matxx>
    : matrices::generic_traits_nonref<nmod_polyxx_ref, nmod_polyxx_srcref> { };

struct nmod_poly_mat_data
{
    typedef nmod_poly_mat_t& data_ref_t;
    typedef const nmod_poly_mat_t& data_srcref_t;

    nmod_poly_mat_t inner;

    nmod_poly_mat_data(slong m, slong n, mp_limb_t modulus)
    {
        nmod_poly_mat_init(inner, m, n, modulus);
    }

    nmod_poly_mat_data(const nmod_poly_mat_data& o)
    {
        nmod_poly_mat_init_set(inner, o.inner);
    }

    nmod_poly_mat_data(nmod_poly_matxx_srcref o)
    {
        nmod_poly_mat_init_set(inner, o._data().inner);
    }

    ~nmod_poly_mat_data() {nmod_poly_mat_clear(inner);}

    template<class Nmod_mat>
    static nmod_poly_mat_data _from_ground(const Nmod_mat& m)
    {
        nmod_poly_mat_data res(m.rows(), m.cols(), m.modulus());
        for(slong i = 0;i < m.rows();++i)
            for(slong j = 0;j < m.cols();++j)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(res.inner, i, j), 0,
                        nmod_mat_entry(m._mat(), i, j));
        return res;
    }
    template<class Nmod_mat>
    static nmod_poly_mat_data from_ground(const Nmod_mat& m,
            typename mp::enable_if<traits::is_nmod_matxx<Nmod_mat> >::type* = 0)
    {
        return _from_ground(m.evaluate());
    }
};
} // detail

// temporary instantiation stuff
FLINTXX_DEFINE_TEMPORARY_RULES(nmod_poly_matxx)

#define NMOD_POLY_MATXX_COND_S FLINTXX_COND_S(nmod_poly_matxx)
#define NMOD_POLY_MATXX_COND_T FLINTXX_COND_T(nmod_poly_matxx)

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, NMOD_POLY_MATXX_COND_T, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_set(to._mat(), from._mat()))

FLINTXX_DEFINE_SWAP(nmod_poly_matxx, nmod_poly_mat_swap(e1._mat(), e2._mat()))

FLINTXX_DEFINE_EQUALS(nmod_poly_matxx, nmod_poly_mat_equal(e1._mat(), e2._mat()))

FLINT_DEFINE_PRINT_PRETTY_COND_2(NMOD_POLY_MATXX_COND_S, const char*,
        (nmod_poly_mat_print(from._mat(), extra), 1))

FLINT_DEFINE_THREEARY_EXPR_COND3(mat_at_op, nmod_polyxx,
        NMOD_POLY_MATXX_COND_S, traits::fits_into_slong, traits::fits_into_slong,
        nmod_poly_set(to._poly(), nmod_poly_mat_entry(e1._mat(), e2, e3)))

FLINT_DEFINE_BINARY_EXPR_COND2(times, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_mul(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, NMOD_POLYXX_COND_S,
        nmod_poly_mat_scalar_mul_nmod_poly(to._mat(), e1._mat(), e2._poly()))
FLINT_DEFINE_CBINARY_EXPR_COND2(times, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, NMODXX_COND_S,
        nmod_poly_mat_scalar_mul_nmod(to._mat(), e1._mat(), e2._limb()))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_add(to._mat(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_sub(to._mat(), e1._mat(), e2._mat()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, nmod_poly_matxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_neg(to._mat(), from._mat()))

namespace rdetail {
inline void nmod_poly_mat_transpose(nmod_poly_mat_t to,
    const nmod_poly_mat_t from)
{
    if(from == to) // guaranteed to be square
    {
        for(slong i = 0;i < nmod_poly_mat_nrows(to) - 1;++i)
          for(slong j = i + 1;j < nmod_poly_mat_ncols(to);++j)
              nmod_poly_swap(nmod_poly_mat_entry(to, i, j),
                      nmod_poly_mat_entry(to, j, i));
    }
    else
    {
        for(slong i = 0;i < nmod_poly_mat_nrows(to);++i)
          for(slong j = 0;j < nmod_poly_mat_ncols(to);++j)
            nmod_poly_set(nmod_poly_mat_entry(to, i, j),
                nmod_poly_mat_entry(from, j, i));
    }
}
}
// TODO update this when nmod_poly_mat has transpose
FLINT_DEFINE_UNARY_EXPR_COND(transpose_op, nmod_poly_matxx, NMOD_POLY_MATXX_COND_S,
        rdetail::nmod_poly_mat_transpose(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(trace_op, nmod_polyxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_trace(to._poly(), from._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(evaluate_op, nmod_matxx,
        NMOD_POLY_MATXX_COND_S, NMODXX_COND_S,
        nmod_poly_mat_evaluate_nmod(to._mat(), e1._mat(), e2._limb()))

#define NMOD_POLY_MATXX_DEFINE_MUL(name) \
FLINT_DEFINE_BINARY_EXPR_COND2(name##_op, nmod_poly_matxx, \
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S, \
        nmod_poly_mat_##name(to._mat(), e1._mat(), e2._mat()))
NMOD_POLY_MATXX_DEFINE_MUL(mul_classical)
NMOD_POLY_MATXX_DEFINE_MUL(mul_KS)
NMOD_POLY_MATXX_DEFINE_MUL(mul_interpolate)

FLINT_DEFINE_UNARY_EXPR_COND(sqr_op, nmod_poly_matxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_sqr(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_KS_op, nmod_poly_matxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_sqr_KS(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_classical_op, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_sqr_classical(to._mat(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(sqr_interpolate_op, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_sqr_interpolate(to._mat(), from._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(pow_op, nmod_poly_matxx,
        NMOD_POLY_MATXX_COND_S, traits::is_unsigned_integer,
        nmod_poly_mat_pow(to._mat(), e1._mat(), e2))

FLINT_DEFINE_UNARY_EXPR_COND(det_op, nmod_polyxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_det(to._poly(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(det_fflu_op, nmod_polyxx, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_det_fflu(to._poly(), from._mat()))
FLINT_DEFINE_UNARY_EXPR_COND(det_interpolate_op, nmod_polyxx,
        NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_det_interpolate(to._poly(), from._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<bool, nmod_poly_matxx, nmod_polyxx>::type >::type
    nmod_poly_mat_inv_rt;
} // rdetail

FLINT_DEFINE_UNARY_EXPR_COND(inv_op, rdetail::nmod_poly_mat_inv_rt,
        NMOD_POLY_MATXX_COND_S,
        to.template get<0>() = nmod_poly_mat_inv(to.template get<1>()._mat(),
            to.template get<2>()._poly(), from._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<slong, nmod_poly_matxx>::type >::type
    nmod_poly_mat_nullspace_rt;
} // rdetail
FLINT_DEFINE_UNARY_EXPR_COND(nullspace_op, rdetail::nmod_poly_mat_nullspace_rt,
        NMOD_POLY_MATXX_COND_S, to.template get<0>() = nmod_poly_mat_nullspace(
            to.template get<1>()._mat(), from._mat()))

FLINT_DEFINE_BINARY_EXPR_COND2(solve_op, rdetail::nmod_poly_mat_inv_rt,
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        to.template get<0>() = nmod_poly_mat_solve(to.template get<1>()._mat(),
            to.template get<2>()._poly(), e1._mat(), e2._mat()))
FLINT_DEFINE_BINARY_EXPR_COND2(solve_fflu_op, rdetail::nmod_poly_mat_inv_rt,
        NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        to.template get<0>() = nmod_poly_mat_solve_fflu(
            to.template get<1>()._mat(),
            to.template get<2>()._poly(), e1._mat(), e2._mat()))

namespace rdetail {
typedef make_ltuple<mp::make_tuple<slong, nmod_poly_matxx, nmod_polyxx>::type>::type
    nmod_poly_matxx_fflu_rt;
} // rdetail

FLINT_DEFINE_THREEARY_EXPR_COND3(fflu_op, rdetail::nmod_poly_matxx_fflu_rt,
        NMOD_POLY_MATXX_COND_S, traits::is_maybe_perm, tools::is_bool,
        to.template get<0>() = nmod_poly_mat_fflu(to.template get<1>()._mat(),
            to.template get<2>()._poly(), maybe_perm_data(e2), e1._mat(), e3))

FLINT_DEFINE_UNARY_EXPR_COND(rref_op, rdetail::nmod_poly_matxx_fflu_rt,
        NMOD_POLY_MATXX_COND_S,
        to.template get<0>() = nmod_poly_mat_rref(to.template get<1>()._mat(),
            to.template get<2>()._poly(), from._mat()))

FLINT_DEFINE_THREEARY_EXPR_COND3(solve_fflu_precomp_op, nmod_poly_matxx,
        traits::is_permxx, NMOD_POLY_MATXX_COND_S, NMOD_POLY_MATXX_COND_S,
        nmod_poly_mat_solve_fflu_precomp(to._mat(), e1._data(),
            e2._mat(), e3._mat()))
} // rules
} // flint

#endif
