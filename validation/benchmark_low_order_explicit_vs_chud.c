#include "ramanujan_c_common.h"

#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CHUD_A 13591409UL
#define CHUD_B 545140134UL
#define CHUD_C3_OVER_24 10939058860032000ULL

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

static int files_equal(const char *a, const char *b) {
    FILE *fa = fopen(a, "rb");
    FILE *fb = fopen(b, "rb");
    if (!fa || !fb) {
        if (fa) fclose(fa);
        if (fb) fclose(fb);
        return 0;
    }
    int ok = 1;
    unsigned char ba[8192], bb[8192];
    for (;;) {
        size_t na = fread(ba, 1, sizeof(ba), fa);
        size_t nb = fread(bb, 1, sizeof(bb), fb);
        if (na != nb || memcmp(ba, bb, na) != 0) {
            ok = 0;
            break;
        }
        if (na == 0) break;
    }
    fclose(fa);
    fclose(fb);
    return ok;
}

static void bs_chud(unsigned long a, unsigned long b,
                    mpz_t P, mpz_t Q, mpz_t T) {
    if (b - a == 1) {
        if (a == 0) {
            mpz_set_ui(P, 1);
            mpz_set_ui(Q, 1);
            mpz_set_ui(T, CHUD_A);
            return;
        }

        mpz_set_ui(P, 6UL * a - 5UL);
        mpz_mul_ui(P, P, 2UL * a - 1UL);
        mpz_mul_ui(P, P, 6UL * a - 1UL);

        mpz_set_ui(Q, a);
        mpz_mul_ui(Q, Q, a);
        mpz_mul_ui(Q, Q, a);
        mpz_mul_ui(Q, Q, (unsigned long)CHUD_C3_OVER_24);

        mpz_set_ui(T, CHUD_B);
        mpz_mul_ui(T, T, a);
        mpz_add_ui(T, T, CHUD_A);
        mpz_mul(T, T, P);
        if (a & 1UL) mpz_neg(T, T);
        return;
    }

    unsigned long m = (a + b) / 2UL;
    mpz_t P1, Q1, T1, P2, Q2, T2, tmp;
    mpz_init(P1); mpz_init(Q1); mpz_init(T1);
    mpz_init(P2); mpz_init(Q2); mpz_init(T2); mpz_init(tmp);

    bs_chud(a, m, P1, Q1, T1);
    bs_chud(m, b, P2, Q2, T2);

    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_mul(tmp, P1, T2);
    mpz_add(T, T, tmp);

    mpz_clear(P1); mpz_clear(Q1); mpz_clear(T1);
    mpz_clear(P2); mpz_clear(Q2); mpz_clear(T2); mpz_clear(tmp);
}

static int chudnovsky_explicit(uint64_t digits, uint64_t guard,
                               const char *path, double *seconds,
                               unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    const long double digits_per_term = 14.1816474627254776555L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard) /
                            digits_per_term) + 2UL;

    mpz_t P, Q, T;
    mpz_init(P); mpz_init(Q); mpz_init(T);
    bs_chud(0, terms, P, Q, T);

    mpf_t fq, ft, root, pi;
    mpf_init2(fq, bits); mpf_init2(ft, bits);
    mpf_init2(root, bits); mpf_init2(pi, bits);
    mpf_set_z(fq, Q);
    mpf_set_z(ft, T);
    mpf_set_ui(root, 10005UL);
    mpf_sqrt(root, root);
    mpf_mul_ui(root, root, 426880UL);
    mpf_mul(pi, fq, root);
    mpf_div(pi, pi, ft);

    FILE *fp = fopen(path, "wb");
    if (!fp) return 0;
    if (!rj_write_decimal_mpf(fp, pi, (size_t)digits)) {
        fclose(fp);
        return 0;
    }
    fputc('\n', fp);
    fclose(fp);

    double t1 = now_seconds();
    if (seconds) *seconds = t1 - t0;
    if (terms_out) *terms_out = terms;

    mpf_clear(fq); mpf_clear(ft); mpf_clear(root); mpf_clear(pi);
    mpz_clear(P); mpz_clear(Q); mpz_clear(T);
    return 1;
}

static int modular_explicit(int p, uint64_t digits, uint64_t guard,
                            const char *path, double *seconds,
                            unsigned long *layers_out,
                            long double *depth_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    RjPoly poly = {0};
    if (!rj_get_z_poly(p, &poly)) return 0;

    RjState a, b;
    rj_state_init(&a, bits);
    rj_state_init(&b, bits);
    rj_initial_state(&a);
    RjState *cur = &a, *nxt = &b;

    const long double target = (long double)digits + (long double)guard;
    unsigned long layers = 0;
    while (minus_log10_abs(cur->z) < target) {
        if (layers >= 128UL) return 0;
        if (!rj_transform_step(&poly, p, 14, cur, nxt)) return 0;
        RjState *sw = cur; cur = nxt; nxt = sw;
        ++layers;
    }

    long double depth = minus_log10_abs(cur->z);
    mpf_t F, thetaF, inv_pi, pi;
    mpf_init2(F, bits); mpf_init2(thetaF, bits);
    mpf_init2(inv_pi, bits); mpf_init2(pi, bits);
    if (!rj_hyper_F_theta(cur->z, (unsigned long)(digits + guard),
                          F, thetaF, 10000000UL)) return 0;
    mpf_mul(inv_pi, cur->u, F);
    mpf_mul(pi, cur->v, thetaF);
    mpf_add(inv_pi, inv_pi, pi);
    mpf_ui_div(pi, 1UL, inv_pi);

    FILE *fp = fopen(path, "wb");
    if (!fp) return 0;
    if (!rj_write_decimal_mpf(fp, pi, (size_t)digits)) {
        fclose(fp);
        return 0;
    }
    fputc('\n', fp);
    fclose(fp);

    double t1 = now_seconds();
    if (seconds) *seconds = t1 - t0;
    if (layers_out) *layers_out = layers;
    if (depth_out) *depth_out = depth;

    mpf_clear(F); mpf_clear(thetaF); mpf_clear(inv_pi); mpf_clear(pi);
    rj_state_clear(&a); rj_state_clear(&b); rj_poly_clear(&poly);
    return 1;
}

static int run_case(uint64_t digits) {
    const uint64_t guard = 160;
    char chud_path[256], p2_path[256], p3_path[256];
    snprintf(chud_path, sizeof(chud_path), "build/bench_chud_%llu.txt",
             (unsigned long long)digits);
    snprintf(p2_path, sizeof(p2_path), "build/bench_p2_%llu.txt",
             (unsigned long long)digits);
    snprintf(p3_path, sizeof(p3_path), "build/bench_p3_%llu.txt",
             (unsigned long long)digits);

    double tc = 0.0, t2 = 0.0, t3 = 0.0;
    unsigned long terms = 0, l2 = 0, l3 = 0;
    long double d2 = 0.0L, d3 = 0.0L;

    if (!chudnovsky_explicit(digits, guard, chud_path, &tc, &terms)) {
        fprintf(stderr, "Chudnovsky failed at %llu digits\n",
                (unsigned long long)digits);
        return 0;
    }
    if (!modular_explicit(2, digits, guard, p2_path, &t2, &l2, &d2)) {
        fprintf(stderr, "p=2 failed at %llu digits\n",
                (unsigned long long)digits);
        return 0;
    }
    if (!modular_explicit(3, digits, guard, p3_path, &t3, &l3, &d3)) {
        fprintf(stderr, "p=3 failed at %llu digits\n",
                (unsigned long long)digits);
        return 0;
    }

    int ok2 = files_equal(chud_path, p2_path);
    int ok3 = files_equal(chud_path, p3_path);
    printf("%llu,%lu,%.9f,%lu,%.12Lg,%.9f,%.6f,%d,%lu,%.12Lg,%.9f,%.6f,%d\n",
           (unsigned long long)digits,
           terms, tc,
           l2, d2, t2, t2 / tc, ok2,
           l3, d3, t3, t3 / tc, ok3);
    fflush(stdout);
    return ok2 && ok3;
}

int main(int argc, char **argv) {
    printf("digits,chud_terms,chud_seconds,p2_layers,p2_depth,p2_seconds,p2_over_chud,p2_equal,p3_layers,p3_depth,p3_seconds,p3_over_chud,p3_equal\n");
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            uint64_t digits = strtoull(argv[i], NULL, 10);
            if (digits == 0 || !run_case(digits)) return 1;
        }
        return 0;
    }

    const uint64_t cases[] = {1000ULL, 10000ULL, 30000ULL, 100000ULL};
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        if (!run_case(cases[i])) return 1;
    }
    return 0;
}
