#include <math.h>
#include "test_helpers.h"
#include "arb.h"

#define PRINT_PRECISION 20

typedef enum {
    ALL_POSITIVE       = 0,
    CONTAINS_ZERO      = 1,
    ALL_NEGATIVE       = 2,
    ANY_FINITE         = 3,         /* used only to represent either of the previous three */
    POSITIVE_INFINITY  = 4,
    NEGATIVE_INFINITY  = 5,
    WHOLE_LINE         = 6,
    ANY_INFINITE       = 7,         /* used only to represent either of the previous three */
    NOT_A_NUMBER       = 8,
    NUM_VALUE_TYPES    = 9          /* used only to count the number of enum values */
} value_type;

typedef struct {
    const char *value;
    value_type type;
} string_with_type;

const string_with_type test_values_arb[] = {
    {"5.3 +/- 1.2",           ALL_POSITIVE},
    {"1.2 +/- 5.3",           CONTAINS_ZERO},
    {"-1.2 +/- 5.3",          CONTAINS_ZERO},
    {"-5.3 +/- 1.2",          ALL_NEGATIVE},
    {"inf +/- 17",            POSITIVE_INFINITY},
    {"-inf +/- 17",           NEGATIVE_INFINITY},
    {"17 +/- inf",            WHOLE_LINE},
    {"-inf +/- inf",          WHOLE_LINE},
    {"nan +/- 17",            NOT_A_NUMBER},
    {"nan +/- inf",           NOT_A_NUMBER},
    {NULL,                    NOT_A_NUMBER},
};

const string_with_type test_values_arf[] = {
    {"5.3",                   ALL_POSITIVE},
    {"0",                     CONTAINS_ZERO},
    {"-5.3",                  ALL_NEGATIVE},
    {"inf",                   POSITIVE_INFINITY},
    {"-inf",                  NEGATIVE_INFINITY},
    {"nan",                   NOT_A_NUMBER},
    {NULL,                    NOT_A_NUMBER},
};

/* Multiplying value_type i with value_type j should yield value_type multiplication[i][j]. */
const value_type multiplication[NUM_VALUE_TYPES][NUM_VALUE_TYPES] = {
    {ALL_POSITIVE,     CONTAINS_ZERO,    ALL_NEGATIVE,     ANY_FINITE,       POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {CONTAINS_ZERO,    CONTAINS_ZERO,    CONTAINS_ZERO,    ANY_FINITE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {ALL_NEGATIVE,     CONTAINS_ZERO,    ALL_POSITIVE,     ANY_FINITE,       NEGATIVE_INFINITY,POSITIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {ANY_FINITE,       ANY_FINITE,       ANY_FINITE,       ANY_FINITE,       ANY_INFINITE,     ANY_INFINITE,     WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {POSITIVE_INFINITY,WHOLE_LINE,       NEGATIVE_INFINITY,ANY_INFINITE,     POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {NEGATIVE_INFINITY,WHOLE_LINE,       POSITIVE_INFINITY,ANY_INFINITE,     NEGATIVE_INFINITY,POSITIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       NOT_A_NUMBER},
    {ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER},
};

/* Adding value_type i to value_type j should yield value_type addition[i][j]. */
const value_type addition[NUM_VALUE_TYPES][NUM_VALUE_TYPES] = {
    {ALL_POSITIVE,     ANY_FINITE,       ANY_FINITE,       ANY_FINITE,       POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {ANY_FINITE,       CONTAINS_ZERO,    ANY_FINITE,       ANY_FINITE,       POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {ANY_FINITE,       ANY_FINITE,       ALL_NEGATIVE,     ANY_FINITE,       POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {ANY_FINITE,       ANY_FINITE,       ANY_FINITE,       ANY_FINITE,       POSITIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {POSITIVE_INFINITY,POSITIVE_INFINITY,POSITIVE_INFINITY,POSITIVE_INFINITY,POSITIVE_INFINITY,WHOLE_LINE,       WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {NEGATIVE_INFINITY,NEGATIVE_INFINITY,NEGATIVE_INFINITY,NEGATIVE_INFINITY,WHOLE_LINE,       NEGATIVE_INFINITY,WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       WHOLE_LINE,       NOT_A_NUMBER},
    {ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     ANY_INFINITE,     WHOLE_LINE,       ANY_INFINITE,     NOT_A_NUMBER},
    {NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER,     NOT_A_NUMBER},
};

int arb_satisfies_value_type(const arb_t x, const value_type t)
{
    switch(t) {
    case ALL_POSITIVE:
        return arb_is_finite(x) && arb_is_positive(x);
    case CONTAINS_ZERO:
        return arb_is_finite(x) && arb_contains_zero(x);
    case ALL_NEGATIVE:
        return arb_is_finite(x) && arb_is_negative(x);
    case ANY_FINITE:
        return arb_is_finite(x);
    case POSITIVE_INFINITY:
        return !arb_is_finite(x) && arb_is_positive(x);
    case NEGATIVE_INFINITY:
        return !arb_is_finite(x) && arb_is_negative(x);
    case WHOLE_LINE:
        return !arb_is_finite(x) && arb_contains_zero(x);
    case ANY_INFINITE:
        return !arb_is_finite(x) && !arf_is_nan(arb_midref(x));
    case NOT_A_NUMBER:
        return arf_is_nan(arb_midref(x));
    default:
        printf("unexpected condition\n");
        flint_abort();
    }
}

void print_arb_and_type(arb_t x, const char *s, const value_type t)
{
    flint_printf("%s: ", s);
    arb_printd(x, PRINT_PRECISION);
    flint_printf(" of type: %d\n", t);
}
void print_arf_and_type(arf_t x, const char *s, const value_type t)
{
    flint_printf("%s: ", s);
    arf_printd(x, PRINT_PRECISION);
    flint_printf(" of type: %d\n", t);
}

TEST_FUNCTION_START(arb_pos_times_posinf, state)
{
    slong i, j, k, prec;
    arb_t t, u, v, w;
    arf_t x;

    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(w);
    arf_init(x);

    for(prec = 20; prec <= 100; prec += 80)
        for(i = 0; test_values_arb[i].value != NULL; ++i)
        {
            for(j = 0; test_values_arb[j].value != NULL; ++j)
            {
                value_type vt = test_values_arb[i].type, vu = test_values_arb[j].type, expected = multiplication[vt][vu];
                arb_set_str(t, test_values_arb[i].value, prec);
                arb_set_str(u, test_values_arb[j].value, prec);

                arb_mul(v, t, u, prec);

                if(!arb_satisfies_value_type(v, expected))
                {
                    flint_printf("FAIL (arb_mul): %s, %s\n", test_values_arb[i].value, test_values_arb[j].value);
                    print_arb_and_type(t, "t", vt);
                    print_arb_and_type(u, "u", vu);
                    flint_printf("v: ");
                    arb_printd(v, PRINT_PRECISION);
                    flint_printf("\nexpected type: %d\n", expected);
                    flint_abort();
                }

                for(k = 0; test_values_arb[k].value != NULL; ++k)
                {
                    value_type vw = test_values_arb[k].type, expected = addition[vw][multiplication[vt][vu]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_addmul(w, t, u, prec);

                    if(!arb_satisfies_value_type(w, expected))
                    {
                        flint_printf("FAIL (arb_addmul): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arb[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arb_and_type(u, "u", vu);
                        flint_printf("w (after): ");
                        arb_printd(w, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }

                    expected = addition[vw][multiplication[ALL_NEGATIVE][multiplication[vt][vu]]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_submul(w, t, u, prec);

                    if(!arb_satisfies_value_type(w, expected))
                    {
                        flint_printf("FAIL (arb_submul): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arb[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arb_and_type(u, "u", vu);
                        flint_printf("w (after): ");
                        arb_printd(w, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }

                    expected = addition[vw][multiplication[vt][vu]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_fma(v, t, u, w, prec);
                    if(!arb_satisfies_value_type(v, expected))
                    {
                        flint_printf("FAIL (arb_fma): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arb[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arb_and_type(u, "u", vu);
                        flint_printf("v (after): ");
                        arb_printd(v, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }
                }
            }

            for(j = 0; test_values_arf[j].value != NULL; ++j)
            {
                value_type vt = test_values_arb[i].type, vx = test_values_arf[j].type;
                arb_set_str(t, test_values_arb[i].value, prec);
                arb_set_str(w, test_values_arf[j].value, prec);
                arf_set(x, arb_midref(w));

                arb_mul_arf(v, t, x, prec);

                if(!arb_satisfies_value_type(v, multiplication[vt][vx]))
                {
                    flint_printf("FAIL (arb_mul_arf): %s, %s\n", test_values_arb[i].value, test_values_arf[j].value);
                    print_arb_and_type(t, "t", vt);
                    print_arf_and_type(x, "x", vx);
                    flint_printf("v: ");
                    arb_printd(v, PRINT_PRECISION);
                    flint_printf("\nexpected type: %d\n", multiplication[vt][vx]);
                    flint_abort();
                }

                for(k = 0; test_values_arb[k].value != NULL; ++k)
                {
                    value_type vw = test_values_arb[k].type, expected = addition[vw][multiplication[vt][vx]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_addmul_arf(w, t, x, prec);

                    if(!arb_satisfies_value_type(w, expected))
                    {
                        flint_printf("FAIL (arb_addmul_arf): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arf[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arf_and_type(x, "x", vx);
                        flint_printf("w (after): ");
                        arb_printd(w, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }

                    expected = addition[vw][multiplication[ALL_NEGATIVE][multiplication[vt][vx]]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_submul_arf(w, t, x, prec);

                    if(!arb_satisfies_value_type(w, expected))
                    {
                        flint_printf("FAIL (arb_submul_arf): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arf[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arf_and_type(x, "x", vx);
                        flint_printf("w (after): ");
                        arb_printd(w, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }

                    expected = addition[vw][multiplication[vt][vx]];
                    arb_set_str(w, test_values_arb[k].value, prec);

                    arb_fma_arf(v, t, x, w, prec);
                    if(!arb_satisfies_value_type(v, expected))
                    {
                        flint_printf("FAIL (arb_fma): %s, %s, %s\n", test_values_arb[k].value, test_values_arb[i].value, test_values_arf[j].value);
                        print_arb_and_type(t, "t", vt);
                        print_arf_and_type(x, "x", vx);
                        flint_printf("v (after): ");
                        arb_printd(v, PRINT_PRECISION);
                        flint_printf("\nexpected type: %d\n", expected);
                        flint_abort();
                    }
                }
            }
        }

    arf_clear(x);
    arb_clear(w);
    arb_clear(v);
    arb_clear(u);
    arb_clear(t);

    TEST_FUNCTION_END(state);
}
