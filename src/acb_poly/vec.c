#include "acb_poly.h"

acb_poly_struct *
_acb_poly_vec_init(slong n)
{
    acb_poly_struct *vec = (acb_poly_struct *) flint_malloc(sizeof(acb_poly_struct) * n);

    for (slong i = 0; i < n; i++)
        acb_poly_init(vec + i);

    return vec;
}

void
_acb_poly_vec_clear(acb_poly_struct *vec, slong n)
{
    for (slong i = 0; i < n; i++)
        acb_poly_clear(vec + i);
    flint_free(vec);
}

void
_acb_poly_vec_set(acb_poly_struct *dest, const acb_poly_struct *src, slong n)
{
    for (slong i = 0; i < n; i++)
        acb_poly_set(dest + i, src + i);
}

void
_acb_poly_vec_set_block(acb_poly_struct *dest, const acb_poly_struct *src,
                        slong n, slong start, slong len)
{
    for (slong k = 0; k < n; k++)
    {
        if (start >= (src + k)->length)
            continue;
        slong len1 = FLINT_MIN(len, (src + k)->length - start);
        acb_poly_fit_length(dest + k, len1);
        _acb_vec_set((dest + k)->coeffs, (src + k)->coeffs + start, len1);
        _acb_poly_set_length(dest + k, len1);
        _acb_poly_normalise(dest + k);
    }
}

void
_acb_poly_vec_fit_length(acb_poly_struct *vec, slong n, slong len)
{
    for (slong i = 0; i < n; i++)
        acb_poly_fit_length(vec + i, len);
}

void
_acb_poly_vec_set_length(acb_poly_struct *vec, slong n, slong len)
{
    for (slong i = 0; i < n; i++)
        _acb_poly_set_length(vec + i, len);
}

void
_acb_poly_vec_normalise(acb_poly_struct *vec, slong n)
{
    for (slong i = 0; i < n; i++)
        _acb_poly_normalise(vec + i);
}

int
_acb_poly_vec_overlaps(acb_poly_struct *vec1, acb_poly_struct *vec2, slong n)
{
    for (slong i = 0; i < n; i++)
        if (!acb_poly_overlaps(vec1 + i, vec2 + i))
            return 0;
    return 1;
}
