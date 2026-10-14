/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arb_poly.h"
#include "fmpz_poly.h"

// f must be preallocated with size at least n+1, where n=a.degree().
// For i=0..n, f[i] will be set to some integer such that 2^f[i] >= a[i].
// return 1 if succeed, 0 if fail (e.g. exponent overflows slong).
// somewhat inefficient, could be made more efficient by not handling COEFF_IS_MPZ() case,
// but this is probably not the bottleneck
static int _arb_poly_exponent_sequence(slong* f, const arb_poly_t a) {
	slong n = arb_poly_degree(a);
	for (slong i = 0; i <= n; ++i) {
		arb_srcptr c = &a->coeffs[i];

		slong mid_exp = arf_abs_bound_lt_2exp_si(arb_midref(c));
		if (mid_exp == ARF_PREC_EXACT)
			return 0;

		slong rad_exp;
		if (mag_is_special(arb_radref(c))) {
			rad_exp = -ARF_PREC_EXACT;
		} else if (mag_is_inf(arb_radref(c))) {
			return 0;
		} else {
			if (fmpz_fits_si(MAG_EXPREF(arb_radref(c))))
				rad_exp = fmpz_get_si(MAG_EXPREF(arb_radref(c)));
			else
				rad_exp = fmpz_sgn(MAG_EXPREF(arb_radref(c))) < 0 ? -ARF_PREC_EXACT : ARF_PREC_EXACT;

			if (rad_exp < -ARF_PREC_EXACT)
				rad_exp = -ARF_PREC_EXACT;
			if (rad_exp > ARF_PREC_EXACT)
				rad_exp = ARF_PREC_EXACT;
			if (rad_exp == ARF_PREC_EXACT)
				return 0;
		}

		f[i] = FLINT_MAX(mid_exp, rad_exp) + 1;
	}

	return 1;
}

// hat_f_indices gets assigned the indices f[i] on the upper convex hull.
// Subsequent entries may get assigned random values.
// f[0..n] (n inclusive) are considered.
// hat_f_indices must have preallocated size at least n+1.
// The returned value is the number of assigned entries.
// caveat: we don't try to check for integer overflow in the multiplication.
static size_t _slong_vec_convex_hull(slong* hat_f_indices, const slong* f, size_t n) {
	size_t hull_size = 1;
	hat_f_indices[0] = 0;
	for (size_t i = 1; i <= n; ++i) {
		while (hull_size >= 2) {
			size_t j1 = hat_f_indices[hull_size-1];
			size_t j2 = hat_f_indices[hull_size-2];
			if ((slong)(i - j2) * (f[j1] - f[j2]) <= (slong)(j1 - j2) * (f[i] - f[j2])) {
				hull_size--;
			} else {
				break;
			}
		}
		hat_f_indices[hull_size++] = i;
	}
	return hull_size;
}

static void _slong_vec_fill_linear_interpolation(slong* f, slong n, int should_round_up) {
	// f[0] and f[n] are input, mutate f[1..n-1]
	if (n <= 1) return;
	if (f[0] == f[n]) {
		for (slong i = 1; i < n; ++i)
			f[i] = f[0];
	} else {
		slong q = (f[n] - f[0]) / n, r = (f[n] - f[0]) % n;
		slong v = f[0], remain = 0;
		if (f[0] < f[n]) {
			for (slong i = 1; i < n; ++i) {
				v += q;
				remain += r;
				if (remain >= n) { remain -= n; ++v; }
				f[i] = v + (remain && should_round_up);
			}
		} else {
			for (slong i = 1; i < n; ++i) {
				v += q;
				remain += r;
				if (remain <= -n) { remain += n; --v; }
				f[i] = v - (remain && !should_round_up);
			}
		}
	}
}

static void _slong_vec_fill_linear_interpolation_round_down(slong* f, slong n) {
	_slong_vec_fill_linear_interpolation(f, n, 0);
}

static void _slong_vec_fill_linear_interpolation_round_up(slong* f, slong n) {
	_slong_vec_fill_linear_interpolation(f, n, 1);
}

static void build_ceil_hat(slong *ceil_hat, const slong *f, const slong *hat_f_indices, size_t hat_f_indices_size)
{
	FLINT_ASSERT(hat_f_indices[0] == 0);
	ceil_hat[0] = f[0];
	slong a = 0;
	for (size_t t = 1; t < hat_f_indices_size; ++t) {
		slong b = hat_f_indices[t];
		ceil_hat[b] = f[b];
		_slong_vec_fill_linear_interpolation_round_up(ceil_hat + a, b - a);
		a = b;
	}
}

// given f[0..n inclusive], g[0..m inclusive],
// hat_f_indices[0..hat_f_indices_size-1],
// hat_g_indices[0..hat_g_indices_size-1],
// compute the convolution of hat_f and hat_g.
static void _slong_vec_convex_maxplus_convolve(
		slong* j_down, // [0..n inclusive]
		// [i] = j^↓_i i.e. find a witness path of pairs (i, j) such that hat_f[i] + hat_g[j] == hat_h[i+j]
		// for all such pairs, then j^↓_i is the minimum j value for each i value
		slong* floor_h, // [0..n+m inclusive]. [k] = floor(hat_h_k)
		const slong* f, size_t n, const slong* hat_f_indices, size_t hat_f_indices_size,
		const slong* g, size_t m, const slong* hat_g_indices, size_t hat_g_indices_size
		) {
	size_t a = 0, b = 0;
	FLINT_ASSERT(hat_f_indices[0] == 0); 
	FLINT_ASSERT(hat_g_indices[0] == 0);
	floor_h[0] = f[0] + g[0];
	j_down[0] = 0;
	while (a < hat_f_indices_size - 1 || b < hat_g_indices_size - 1) {
		if (a < hat_f_indices_size - 1 && (b >= hat_g_indices_size - 1 || 
					// note: multiplication overflow not checked!
					(f[hat_f_indices[a+1]] - f[hat_f_indices[a]]) * (hat_g_indices[b+1] - hat_g_indices[b]) >
					(g[hat_g_indices[b+1]] - g[hat_g_indices[b]]) * (hat_f_indices[a+1] - hat_f_indices[a]))) {
			floor_h[hat_f_indices[a+1] + hat_g_indices[b]] = f[hat_f_indices[a+1]] + g[hat_g_indices[b]];
			slong j = hat_g_indices[b];
			for (slong i = hat_f_indices[a] + 1; i <= hat_f_indices[a+1]; ++i)
				j_down[i] = j;
			_slong_vec_fill_linear_interpolation_round_down(floor_h + hat_f_indices[a] + hat_g_indices[b],
					hat_f_indices[a+1] - hat_f_indices[a]);
			++a;
		} else {
			floor_h[hat_f_indices[a] + hat_g_indices[b+1]] = f[hat_f_indices[a]] + g[hat_g_indices[b+1]];
			_slong_vec_fill_linear_interpolation_round_down(floor_h + hat_f_indices[a] + hat_g_indices[b],
					hat_g_indices[b+1] - hat_g_indices[b]);
			++b;
		}
	}
}

static slong lower_d(const slong* floor_h, const slong* ceil_hat_f, const slong* ceil_hat_g,
		size_t i, size_t j) { return floor_h[i+j] - ceil_hat_f[i] - ceil_hat_g[j]; }

static int is_below_witness_path(const slong* j_down, size_t i, size_t j) { return j < j_down[i]; }
static int is_above_witness_path(const slong* j_down, size_t i, size_t j) { return j > j_down[i]; }

static int _arb_poly_should_add_product_early_return(
		size_t ia, size_t ib, size_t ja, size_t jb,
		const slong* j_down, const slong* floor_h, const slong* ceil_hat_f, const slong* ceil_hat_g, slong prec
		) {
	if (is_below_witness_path(j_down, ia, jb) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ia, jb) >= prec) return 1;
	if (is_above_witness_path(j_down, ib, ja) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ib, ja) >= prec) return 1;
	return 0;
}

static int _arb_poly_are_all_coefficients_nonnegative(const arb_poly_t a)
{
	for (slong i = 0; i < arb_poly_length(a); ++i)
		if (!arb_is_nonnegative(&a->coeffs[i])) return 0;
	return 1;
}

// underestimate by at most 1
// output vectors alloc guaranteed to be at least length of input
// lo[i] is guaranteed to be <= up[i]
static void _arb_poly_lower_and_upper_underestimate(fmpz_poly_t lo, fmpz_poly_t up, const arb_poly_t a)
{
	slong len = arb_poly_length(a);
	fmpz_poly_fit_length(lo, len);
	fmpz_poly_fit_length(up, len);

	arf_t L, U;
	arf_init(L);
	arf_init(U);
	for (slong i = 0; i < len; ++i) {
		if (fmpz_cmp_si(MAG_EXPREF(arb_radref(&a->coeffs[i])), 0) < 0) {
			// arb_radref(&a->coeffs[i]) < 2^-1 = 0.5
			arf_get_fmpz(&lo->coeffs[i], arb_midref(&a->coeffs[i]), ARF_RND_NEAR);
			fmpz_set(&up->coeffs[i], &lo->coeffs[i]);
		} else {
			// in this branch, lo[i] <= up[i] is guaranteed
			slong tmp = arf_abs_bound_lt_2exp_si(arb_midref(&a->coeffs[i]));
			arb_get_interval_arf(L, U, &a->coeffs[i], FLINT_MAX(tmp, 0) + 2);
			arf_get_fmpz(&lo->coeffs[i], L, ARF_RND_UP); // because underestimate
			arf_get_fmpz(&up->coeffs[i], U, ARF_RND_DOWN);
			if (fmpz_cmp(&lo->coeffs[i], &up->coeffs[i]) > 0) // slightly inefficient
				fmpz_set(&lo->coeffs[i], &up->coeffs[i]);
		}
		FLINT_ASSERT(fmpz_cmp(&lo->coeffs[i], &up->coeffs[i]) <= 0);
	}
	arf_clear(L);
	arf_clear(U);
	_fmpz_poly_set_length(lo, len);
	_fmpz_poly_set_length(up, len);
	_fmpz_poly_normalise(lo);
	_fmpz_poly_normalise(up);
}

// underestimate by at most 0.5
// output vectors alloc guaranteed to be at least length of input
static void arb_poly_get_mid_rad_fmpz_poly(fmpz_poly_t mid, fmpz_poly_t rad, const arb_poly_t a)
{
	slong len = arb_poly_length(a);
	fmpz_poly_fit_length(mid, len);
	fmpz_poly_fit_length(rad, len);

	for (slong i = 0; i < len; ++i) {
		arb_srcptr c = &a->coeffs[i];
		arf_get_fmpz(&mid->coeffs[i], arb_midref(c), ARF_RND_NEAR);
		mag_get_fmpz(&rad->coeffs[i], arb_radref(c));
	}

	_fmpz_poly_set_length(mid, len);
	_fmpz_poly_set_length(rad, len);
	_fmpz_poly_normalise(mid);
	_fmpz_poly_normalise(rad);
}

// compare: _arb_poly_addmullow_rad
static void _arb_poly_convolve_using_fmpz_multiply(arb_poly_t c, const arb_poly_t a, const arb_poly_t b, slong prec)
{
	slong len = arb_poly_degree(a) + arb_poly_degree(b) + 1;
	arb_poly_fit_length(c, len);

	if (_arb_poly_are_all_coefficients_nonnegative(a) && _arb_poly_are_all_coefficients_nonnegative(b)) {
		// use only 2 multiplications (both are high-precision. Not sure if it's actually faster)
		fmpz_poly_t a_lower, a_upper, b_lower, b_upper, c_lower, c_upper;

		fmpz_poly_init(a_lower);
		fmpz_poly_init(a_upper);
		fmpz_poly_init(b_lower);
		fmpz_poly_init(b_upper);
		_arb_poly_lower_and_upper_underestimate(a_lower, a_upper, a);
		_arb_poly_lower_and_upper_underestimate(b_lower, b_upper, b);

		for (slong i = 0; i < fmpz_poly_length(a_lower); ++i) {
			FLINT_ASSERT(fmpz_cmp_si(&a_lower->coeffs[i], 0) >= 0);
			FLINT_ASSERT(fmpz_cmp(&a_lower->coeffs[i], &a_upper->coeffs[i]) <= 0);
		}
		for (slong i = 0; i < fmpz_poly_length(b_lower); ++i) {
			FLINT_ASSERT(fmpz_cmp_si(&b_lower->coeffs[i], 0) >= 0);
			FLINT_ASSERT(fmpz_cmp(&b_lower->coeffs[i], &b_upper->coeffs[i]) <= 0);
		}

		fmpz_poly_init(c_lower);
		fmpz_poly_mul(c_lower, a_lower, b_lower);
		fmpz_poly_clear(a_lower);
		fmpz_poly_clear(b_lower);

		fmpz_poly_init(c_upper);
		fmpz_poly_mul(c_upper, a_upper, b_upper);
		fmpz_poly_clear(a_upper);
		fmpz_poly_clear(b_upper);

		arf_t lower_arf, upper_arf;
		arf_init(lower_arf);
		arf_init(upper_arf);
		slong len1 = fmpz_poly_length(c_lower);
		FLINT_ASSERT(len1 <= fmpz_poly_length(c_upper));
		for (slong i = 0; i < len1; ++i) {
			arf_set_fmpz(lower_arf, &c_lower->coeffs[i]);
			arf_set_fmpz(upper_arf, &c_upper->coeffs[i]);
			arb_set_interval_arf(&c->coeffs[i], lower_arf, upper_arf, prec);
		}
		arf_zero(lower_arf);
		for (slong i = len1; i < fmpz_poly_length(c_upper); ++i) {
			arf_set_fmpz(upper_arf, &c_upper->coeffs[i]);
			arb_set_interval_arf(&c->coeffs[i], lower_arf, upper_arf, prec);
		}
		fmpz_poly_clear(c_lower);
		fmpz_poly_clear(c_upper);
		arf_clear(lower_arf);
		arf_clear(upper_arf);
	} else {
		// use 1 high-precision multiplication and 2 low-precision multiplications
		/* (xm + xr)*(ym + yr) = (xm*ym) + (xr*ym + xm*yr + xr*yr)
						   = (xm*ym) + (xm*yr + xr*(ym + yr))  */
		fmpz_poly_t a_mid, a_rad, b_mid, b_rad, c_mid, c_rad;

		fmpz_poly_init(a_mid);
		fmpz_poly_init(a_rad);
		fmpz_poly_init(b_mid);
		fmpz_poly_init(b_rad);
		fmpz_poly_init(c_mid);
		fmpz_poly_init(c_rad);

		arb_poly_get_mid_rad_fmpz_poly(a_mid, a_rad, a);
		arb_poly_get_mid_rad_fmpz_poly(b_mid, b_rad, b);

		fmpz_poly_mul(c_mid, a_mid, b_mid);

		for (slong i = 0; i < fmpz_poly_length(a_mid); ++i)
			fmpz_abs(&a_mid->coeffs[i], &a_mid->coeffs[i]);
		for (slong i = 0; i < fmpz_poly_length(b_mid); ++i)
			fmpz_abs(&b_mid->coeffs[i], &b_mid->coeffs[i]);

		fmpz_poly_mul(c_rad, a_mid, b_rad);
		fmpz_poly_add(a_mid, b_mid, b_rad);
		fmpz_poly_clear(b_rad);
		fmpz_poly_mul(b_mid, a_mid, a_rad);
		fmpz_poly_clear(a_mid);
		fmpz_poly_clear(a_rad);
		fmpz_poly_add(c_rad, c_rad, b_mid);
		fmpz_poly_clear(b_mid);

		slong len1 = fmpz_poly_length(c_mid);
		FLINT_ASSERT(len1 <= arb_poly_length(c));
		for (slong i = 0; i < len1; ++i)
			arf_set_fmpz(arb_midref(&c->coeffs[i]), &c_mid->coeffs[i]);

		len1 = fmpz_poly_length(c_rad);
		FLINT_ASSERT(len1 <= arb_poly_length(c));
		for (slong i = 0; i < len1; ++i)
			mag_set_fmpz(arb_radref(&c->coeffs[i]), &c_rad->coeffs[i]);

		fmpz_poly_clear(c_mid);
		fmpz_poly_clear(c_rad);
	}

	_arb_poly_set_length(c, len);
	_arb_poly_normalise(c);
}

// given as input double a, find a real number numer/denom such that
// 2^(numer/denom*k) can be computed quickly.
// Internal implementation detail: Since multiplication by a power of 2 is fast,
// we want denom to be small, precompute 2^(numer/denom*k) for 0 <= k < denom
typedef struct nice_powers_struct {
	slong numer, denom /* always positive */;
	arb_struct* powers; // size = denom
} nice_powers_struct;

static double f_prime(size_t i, size_t ja, size_t jb,
		const slong *ceil_hat_f, const slong *ceil_hat_g, const slong *floor_h, slong prec)
{
	slong mind = LONG_MAX;
	for (size_t j = ja; j <= jb; ++j) {
		slong d = lower_d(floor_h, ceil_hat_f, ceil_hat_g, i, j);
		if (d < mind) mind = d;
	}
	return (double)(ceil_hat_f[i] - prec + mind - 2);
}

static double g_prime(size_t j, size_t ia, size_t ib,
		const slong *ceil_hat_f, const slong *ceil_hat_g, const slong *floor_h, slong prec)
{
	slong mind = LONG_MAX;
	for (size_t i = ia; i <= ib; ++i) {
		slong d = lower_d(floor_h, ceil_hat_f, ceil_hat_g, i, j);
		if (d < mind) mind = d;
	}
	return (double)(ceil_hat_g[j] - prec + mind - 2);
}

// max_error may not be exactly satisfied, unlike ball arithmetic
static void nice_powers_init(nice_powers_struct* data, double a, double max_error, slong prec) {
	FLINT_ASSERT(max_error > 0);
	FLINT_ASSERT(ceil(0.5 / max_error) < ((ulong) -1) / 2);
	data->denom = ceil(0.5 / max_error);
	data->numer = round(a * data->denom);
	// TODO use best_rational_fast should be better
	data->powers = _arb_vec_init(data->denom);

	arb_t two_pow_numer_over_denom, two;
	arb_init(two_pow_numer_over_denom);
	arb_set_si(two_pow_numer_over_denom, data->numer);
	arb_div_si(two_pow_numer_over_denom, two_pow_numer_over_denom, data->denom, prec);
	arb_init(two);
	arb_set_ui(two, 2);
	arb_pow(two_pow_numer_over_denom, two, two_pow_numer_over_denom, prec);
	arb_clear(two);

	_arb_vec_set_powers(data->powers, two_pow_numer_over_denom, data->denom, prec);
	arb_clear(two_pow_numer_over_denom);
}

static void nice_powers_clear(nice_powers_struct* data) {
	_arb_vec_clear(data->powers, data->denom);
}

// y = x * 2^(numer/denom*k + extra). Allow aliasing between x and y
static void nice_powers_mul(arb_t y, nice_powers_struct *data, arb_srcptr x, slong k, slong prec, slong extra)
{
	slong q = k / data->denom;
	slong r = k % data->denom;
	if (r < 0) { q -= 1; r += data->denom; }
	arb_mul_2exp_si(y, x, q * data->numer + extra); // because of aliasing, x may be invalid now
	if (r) {
		arb_mul(y, y, data->powers + r, prec);
	}
}

static void _arb_poly_mullow_add_product_fft(
		arb_poly_t c,
		size_t ia, size_t ib, size_t ja, size_t jb,
		const arb_poly_t a, const arb_poly_t b,
		const slong* ceil_hat_f, const slong* ceil_hat_g,
		const slong* floor_h, slong prec) {
	double fpia /* f'(i_A) */ = f_prime(ia, ja, jb, ceil_hat_f, ceil_hat_g, floor_h, prec);
	double fpib = f_prime(ib, ja, jb, ceil_hat_f, ceil_hat_g, floor_h, prec);
	double gpja = g_prime(ja, ia, ib, ceil_hat_f, ceil_hat_g, floor_h, prec);
	double gpjb = g_prime(jb, ia, ib, ceil_hat_f, ceil_hat_g, floor_h, prec);
	double lc_a /* lowercase a */ = (fpib - fpia + gpjb - gpja) / (double)(ib - ia + jb - ja);
	nice_powers_struct data;
	nice_powers_init(&data, lc_a, 1. / (double)(ib - ia + jb - ja), prec);
	lc_a = (double)data.numer / (double)data.denom * (data.numer > 0 ? 1 + 1e-15 : 1 - 1e-15);
	slong lc_b = floor(FLINT_MIN(fpia - lc_a * ia, fpib - lc_a * ib) * (1 - 1e-15)),
		  lc_c = floor(FLINT_MIN(gpja - lc_a * ja, gpjb - lc_a * jb) * (1 - 1e-15));
	// might be an underestimation, but must not be an overestimation

	arb_poly_t a_prime, b_prime;
	arb_poly_init2(a_prime, ib - ia + 1);
	arb_poly_init2(b_prime, jb - ja + 1);

	for (size_t i = ia; i <= ib; ++i)
		// a_prime[i-ia] = a[i] / 2^(lc_a*i+lc_b)
		nice_powers_mul(&a_prime->coeffs[i - ia], &data, &a->coeffs[i], -(slong)i, prec, -lc_b);
	for (size_t j = ja; j <= jb; ++j)
		// b_prime[j-ja] = b[j] / 2^(lc_a*j+lc_c)
		nice_powers_mul(&b_prime->coeffs[j - ja], &data, &b->coeffs[j], -(slong)j, prec, -lc_c);

	_arb_poly_set_length(a_prime, ib - ia + 1);
	_arb_poly_set_length(b_prime, jb - ja + 1);
	_arb_poly_normalise(a_prime);
	_arb_poly_normalise(b_prime);

	arb_poly_t c_prime;
	arb_poly_init2(c_prime, ib - ia + jb - ja + 1);
	_arb_poly_convolve_using_fmpz_multiply(c_prime, a_prime, b_prime, prec);
	arb_poly_clear(a_prime);
	arb_poly_clear(b_prime);

	arb_t tmp;
	arb_init(tmp);
	for (size_t k = ia + ja; k <= ib + jb; ++k)
	{
		// c[k] += c_prime[k - ia - ja] * 2^(a*k + lc_b + lc_c)
		nice_powers_mul(tmp, &data, &c_prime->coeffs[k - ia - ja], (slong)k, prec, lc_b + lc_c);
		arb_add(&c->coeffs[k], &c->coeffs[k], tmp, prec);
	}
	arb_clear(tmp);
	nice_powers_clear(&data);
	arb_poly_clear(c_prime);
}

static void _arb_poly_mullow_add_product(
		arb_poly_t c,
		const arb_poly_t a, const arb_poly_t b,
		size_t ia, size_t ib, size_t ja, size_t jb,
		const slong* j_down, const slong* floor_h, const slong* ceil_hat_f, const slong* ceil_hat_g, slong prec) {
	FLINT_ASSERT(ia <= ib);
	FLINT_ASSERT(ja <= jb);
	if (_arb_poly_should_add_product_early_return(ia, ib, ja, jb, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec)) return;
	if (is_below_witness_path(j_down, ia, ja) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ia, ja) >= prec) {
		size_t jm;
		// find max jm >= ja such that is_below_witness_path(ia, jm) && lower_d(ia, jm) >= prec
		jm = jb;
		while (!(is_below_witness_path(j_down, ia, jm) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ia, jm) >= prec)) {
			jm--;
		}
		FLINT_ASSERT(_arb_poly_should_add_product_early_return(ia, ib, ja, jm, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec));
		_arb_poly_mullow_add_product(c, a, b, ia, ib, jm + 1, jb, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
		return;
	}
	if (is_above_witness_path(j_down, ib, jb) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ib, jb) >= prec) {
		size_t jm;
		// find min jm <= jb such that is_above_witness_path(ib, jm) && lower_d(ib, jb) >= prec
		jm = ja;
		while (!(is_above_witness_path(j_down, ib, jm) && lower_d(floor_h, ceil_hat_f, ceil_hat_g, ib, jm) >= prec)) {
			jm++;
		}
		_arb_poly_mullow_add_product(c, a, b, ia, ib, ja, jm - 1, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
		FLINT_ASSERT(_arb_poly_should_add_product_early_return(ia, ib, jm, jb, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec));
		return;
	}
	const size_t NAIVE_THRESHOLD = 3;
	if (ib - ia <= NAIVE_THRESHOLD || jb - ja <= NAIVE_THRESHOLD) {
		arb_t tmp;
		arb_init(tmp);

		for (size_t i = ia; i <= ib; ++i)
			for (size_t j = ja; j <= jb; ++j) {
				// c[i+j] = c[i+j] + a[i] * b[j]
				arb_mul(tmp, &a->coeffs[i], &b->coeffs[j], prec);
				arb_add(&c->coeffs[i + j], &c->coeffs[i + j], tmp, prec);
			}

		arb_clear(tmp);
		return;
	}
	if (FLINT_MAX(lower_d(floor_h, ceil_hat_f, ceil_hat_g, ia, jb), lower_d(floor_h, ceil_hat_f, ceil_hat_g, ib, ja)) <= 3*prec) {
		_arb_poly_mullow_add_product_fft(c, ia, ib, ja, jb, a, b, ceil_hat_f, ceil_hat_g, floor_h, prec);
		return;
	}
	size_t im = (ia+ib)/2, jm = (ja+jb)/2;
	_arb_poly_mullow_add_product(c, a, b, ia, im, ja, jm, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
	_arb_poly_mullow_add_product(c, a, b, im+1, ib, ja, jm, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
	_arb_poly_mullow_add_product(c, a, b, ia, im, jm+1, jb, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
	_arb_poly_mullow_add_product(c, a, b, im+1, ib, jm+1, jb, j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);
}

// main function, compute c = a * b. Despite the name this doesn't yet support specifying output length.
void arb_poly_mullow_convex(arb_poly_t c, const arb_poly_t a, const arb_poly_t b, slong prec) {
	if (arb_poly_length(a) == 0 || arb_poly_length(b) == 0) {
		arb_poly_zero(c);
		return;
	}

	slong *f = FLINT_ARRAY_ALLOC(arb_poly_length(a), slong);
	_arb_poly_exponent_sequence(f, a);

	slong *hat_f_indices = FLINT_ARRAY_ALLOC(arb_poly_length(a), slong);
	size_t hat_f_indices_size = _slong_vec_convex_hull(hat_f_indices, f, arb_poly_degree(a));
	slong *ceil_hat_f = FLINT_ARRAY_ALLOC(arb_poly_length(a), slong);
	build_ceil_hat(ceil_hat_f, f, hat_f_indices, hat_f_indices_size);

	slong *g = FLINT_ARRAY_ALLOC(arb_poly_length(b), slong);
	_arb_poly_exponent_sequence(g, b);

	slong *hat_g_indices = FLINT_ARRAY_ALLOC(arb_poly_length(b), slong);
	size_t hat_g_indices_size = _slong_vec_convex_hull(hat_g_indices, g, arb_poly_degree(b));
	slong *ceil_hat_g = FLINT_ARRAY_ALLOC(arb_poly_length(b), slong);
	build_ceil_hat(ceil_hat_g, g, hat_g_indices, hat_g_indices_size);

	slong *j_down = FLINT_ARRAY_ALLOC(arb_poly_length(a), slong);
	slong *floor_h = FLINT_ARRAY_ALLOC(arb_poly_degree(a) + arb_poly_degree(b) + 1, slong);

	_slong_vec_convex_maxplus_convolve(j_down, floor_h,
			f, arb_poly_degree(a), hat_f_indices, hat_f_indices_size,
			g, arb_poly_degree(b), hat_g_indices, hat_g_indices_size);

	slong len = arb_poly_degree(a) + arb_poly_degree(b) + 1;
	arb_poly_fit_length(c, len);
	{
		mag_t t, tmp;
		mag_init(t);
		mag_init(tmp);
		mag_set_ui(t, (ulong) FLINT_MIN(arb_poly_length(a), arb_poly_length(b)));
		for (size_t k = 0; k < (size_t)len; ++k) {
			arb_zero(&c->coeffs[k]);
			mag_set(tmp, t);
			mag_mul_2exp_si(tmp, tmp, floor_h[k] + 1 - prec);
			arb_add_error_mag(&c->coeffs[k], tmp);
		}
		mag_clear(t);
		mag_clear(tmp);
	}

	_arb_poly_mullow_add_product(c, a, b, 0, arb_poly_degree(a), 0, arb_poly_degree(b), j_down, floor_h, ceil_hat_f, ceil_hat_g, prec);

	_arb_poly_set_length(c, len);
	_arb_poly_normalise(c);

	flint_free(f);
	flint_free(hat_f_indices);
	flint_free(g);
	flint_free(hat_g_indices);
	flint_free(j_down);
	flint_free(floor_h);
	flint_free(ceil_hat_f);
	flint_free(ceil_hat_g);
}
