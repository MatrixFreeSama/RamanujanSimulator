#include "ramanujan_c_common.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

static long double minus_log10_abs(const mpf_t z) {
    if (mpf_sgn(z) == 0) return INFINITY;
    mpf_t a;
    mpf_init2(a, mpf_get_prec(z));
    mpf_abs(a, z);
    mp_exp_t e2 = 0;
    double m = mpf_get_d_2exp(&e2, a);
    mpf_clear(a);
    return -(log10l(fabsl((long double)m)) +
             (long double)e2 * log10l(2.0L));
}

static int cmp_double(const void *aa, const void *bb) {
    const double a = *(const double *)aa;
    const double b = *(const double *)bb;
    return (a > b) - (a < b);
}

static double timed_step_median(const RjPoly *poly, int p,
                                const RjState *in, RjState *out,
                                int newton_iters, int repeats) {
    double *samples = (double *)calloc((size_t)repeats, sizeof(double));
    if (!samples) return -1.0;
    for (int r = 0; r < repeats; ++r) {
        double t0 = now_seconds();
        if (!rj_transform_step(poly, p, newton_iters, in, out)) {
            free(samples);
            return -1.0;
        }
        double t1 = now_seconds();
        samples[r] = t1 - t0;
    }
    qsort(samples, (size_t)repeats, sizeof(double), cmp_double);
    double med = samples[repeats / 2];
    free(samples);
    return med;
}

static int run_chain(int p, unsigned long dps, int layers,
                     int newton_iters, int repeats, const char *label) {
    RjPoly poly = {0};
    if (!rj_get_z_poly(p, &poly)) {
        fprintf(stderr, "no frozen z-polynomial for p=%d\n", p);
        return 0;
    }

    mp_bitcnt_t bits = rj_bits_for_decimal(dps, 160);
    RjState a, b;
    rj_state_init(&a, bits);
    rj_state_init(&b, bits);
    rj_initial_state(&a);
    RjState *cur = &a, *nxt = &b;

    printf("# %s p=%d dps=%lu bits=%lu newton=%d repeats=%d terms=%zu\n",
           label, p, dps, (unsigned long)bits, newton_iters, repeats, poly.nterms);
    printf("label,p,layer,depth_before,depth_after,gain_digits,multiplier,median_seconds,gain_digits_per_second,log_multiplier_per_second\n");

    for (int k = 1; k <= layers; ++k) {
        long double before = minus_log10_abs(cur->z);
        double sec = timed_step_median(&poly, p, cur, nxt, newton_iters, repeats);
        if (sec <= 0.0) {
            fprintf(stderr, "transform timing failure p=%d layer=%d\n", p, k);
            rj_state_clear(&a); rj_state_clear(&b); rj_poly_clear(&poly);
            return 0;
        }
        long double after = minus_log10_abs(nxt->z);
        long double gain = after - before;
        long double mult = after / before;
        long double gps = gain / (long double)sec;
        long double eta = logl(mult) / (long double)sec;
        printf("%s,%d,%d,%.18Lg,%.18Lg,%.18Lg,%.12Lg,%.9g,%.12Lg,%.12Lg\n",
               label, p, k, before, after, gain, mult, sec, gps, eta);
        RjState *sw = cur; cur = nxt; nxt = sw;
    }

    rj_state_clear(&a);
    rj_state_clear(&b);
    rj_poly_clear(&poly);
    return 1;
}

int main(void) {
    /*
     * Baseline: repeated p=2 transforms at one fixed working precision.
     * This isolates whether gain-per-wall-time grows as depth doubles.
     */
    if (!run_chain(2, 6000UL, 8, 14, 5, "p2_chain")) return 1;

    /*
     * Dimension sweep: same seed, same precision, same Newton budget.
     * p=2,3,5 are the frozen exact modular dimensions currently present in
     * the repository. Three layers keep all runs well inside the chosen
     * 6000-decimal-digit working precision.
     */
    if (!run_chain(2, 6000UL, 3, 14, 5, "dimension_sweep")) return 1;
    if (!run_chain(3, 6000UL, 3, 14, 5, "dimension_sweep")) return 1;
    if (!run_chain(5, 6000UL, 3, 14, 5, "dimension_sweep")) return 1;
    return 0;
}
