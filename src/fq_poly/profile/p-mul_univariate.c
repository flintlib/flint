#include "flint.h"
#include "fq_poly.h"
#include "profiler.h"

#define nalgs 2
#define ncases 10
#define cpumin 2

int
main(int argc, char** argv)
{
    double s[nalgs];

    int c, n, lenf, leng, ext, reps = 0;
    fmpz_t p, temp;
    fq_poly_t f, g, h;
    fq_ctx_t ctx;
    
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_set_str(p, argv[1], 10);

    fmpz_init(temp);
       
    fmpz_set_str(temp, argv[2], 10);
    ext = fmpz_get_si(temp);

    lenf = atol(argv[3]);
    leng = atol(argv[4]);

    fq_ctx_init(ctx, p, ext, "a");

    fq_poly_init(f, ctx);
    fq_poly_init(g, ctx);
    fq_poly_init(h, ctx);

    for (c = 0; c < nalgs; c++)
    {
        s[c] = 0.0;
    }
       
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int l, loops = 1;

        /*
           Construct random elements of fq
        */
        {
            fq_poly_randtest_monic(f, state, lenf, ctx);
            fq_poly_randtest_monic(g, state, leng, ctx);
        }
        
    loop:
        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            fq_poly_mul_classical(h, f, g, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (l = 0; l < loops; l++)
        {
            fq_poly_mul_univariate(h, f, g, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        for (c = 0; c < nalgs; c++)
            if (t[c] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 10;
                goto loop;
            }
        
        for (c = 0; c < nalgs; c++)
            s[c] += t[c];
        reps += loops;
    }
        
    for (c = 0; c < nalgs; c++)
    {
        flint_printf("%20f ", s[c] / (double) reps);
        fflush(stdout);
    }
    printf("\n");
        
    fq_poly_clear(h, ctx);
    fq_poly_clear(f, ctx);
    fq_poly_clear(g, ctx);
    fq_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(temp);

    FLINT_TEST_CLEANUP(state);
    
    return 0;
}
