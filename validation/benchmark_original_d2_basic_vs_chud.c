#include "ramanujan_c_common.h"

#include <gmp.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define CHUD_A 13591409UL
#define CHUD_B 545140134UL
#define CHUD_C3_OVER_24 10939058860032000ULL

/*
 * Original low-dimensional plug-in family member:
 *   dimension d = h(Delta) = 2, Delta = -427.
 *
 * H(J) = J^2 + H1 J + H2
 * L(U) = U^2 + L1 U + L2
 *
 * The evaluator uses only these two quadratic algebraic roots and the common
 * hypergeometric series.  It does NOT use modular-polynomial state transport.
 */
static const char *D2_H1 = "15611455512523783919812608000";
static const char *D2_H2 = "155041756222618916546936832000000";
static const char *D2_L1 = "6049980860956530737897555251200";
static const char *D2_L2 = "9989238519497195119784714748139929600";

static double now_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

static int files_equal(const char *a, const char *b) {
    FILE *fa = fopen(a, "rb");
    FILE *fb = fopen(b, "rb");
    if (!fa || !fb) {
        if (fa) fclose(fa);
        if (fb) fclose(fb);
        return 0;
    }

    unsigned char ba[8192], bb[8192];
    int ok = 1;
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

/* Stable tiny root of 1 + a1*x + a2*x^2 = 0. */
static int tiny_quadratic_root(const char *a1s, const char *a2s,
                               mpf_t x, mp_bitcnt_t bits) {
    mpz_t a1z, a2z, discz, tmpz;
    mpz_inits(a1z, a2z, discz, tmpz, NULL);
    if (mpz_set_str(a1z, a1s, 10) != 0 || mpz_set_str(a2z, a2s, 10) != 0) {
        mpz_clears(a1z, a2z, discz, tmpz, NULL);
        return 0;
    }

    mpz_mul(discz, a1z, a1z);
    mpz_mul_ui(tmpz, a2z, 4UL);
    mpz_sub(discz, discz, tmpz);
    if (mpz_sgn(discz) <= 0) {
        mpz_clears(a1z, a2z, discz, tmpz, NULL);
        return 0;
    }

    mpf_t a1, disc, root, den;
    mpf_init2(a1, bits);
    mpf_init2(disc, bits);
    mpf_init2(root, bits);
    mpf_init2(den, bits);
    mpf_set_z(a1, a1z);
    mpf_set_z(disc, discz);
    mpf_sqrt(root, disc);
    mpf_add(den, a1, root);
    mpf_set_si(x, -2);
    mpf_div(x, x, den);

    mpf_clears(a1, disc, root, den, NULL);
    mpz_clears(a1z, a2z, discz, tmpz, NULL);
    return 1;
}

static void hyper_term_advance(mpf_t term, const mpf_t z, unsigned long n) {
    unsigned long k = n - 1UL;
#if ULONG_MAX > 0xffffffffUL
    unsigned long long num =
        (unsigned long long)(6UL * k + 1UL) *
        (unsigned long long)(2UL * k + 1UL) *
        (unsigned long long)(6UL * k + 5UL);
    unsigned long long den =
        72ULL * (unsigned long long)n * (unsigned long long)n * (unsigned long long)n;
    mpf_mul_ui(term, term, (unsigned long)num);
    mpf_div_ui(term, term, (unsigned long)den);
#else
    mpf_mul_ui(term, term, 6UL * k + 1UL);
    mpf_mul_ui(term, term, 2UL * k + 1UL);
    mpf_mul_ui(term, term, 6UL * k + 5UL);
    mpf_div_ui(term, term, 72UL);
    mpf_div_ui(term, term, n);
    mpf_div_ui(term, term, n);
    mpf_div_ui(term, term, n);
#endif
    mpf_mul(term, term, z);
}

static int write_pi(const char *path, const mpf_t pi, uint64_t digits) {
    FILE *fp = fopen(path, "wb");
    if (!fp) return 0;
    int ok = rj_write_decimal_mpf(fp, pi, (size_t)digits);
    if (ok) fputc('\n', fp);
    fclose(fp);
    return ok;
}

static int d2_basic_pi(uint64_t digits, uint64_t guard, const char *path,
                       double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    mpf_t y, v, z, alpha, beta, K, term, F, T, tmp, inv_pi, pi;
    mpf_init2(y, bits);      mpf_init2(v, bits);
    mpf_init2(z, bits);      mpf_init2(alpha, bits);
    mpf_init2(beta, bits);   mpf_init2(K, bits);
    mpf_init2(term, bits);   mpf_init2(F, bits);
    mpf_init2(T, bits);      mpf_init2(tmp, bits);
    mpf_init2(inv_pi, bits); mpf_init2(pi, bits);

    if (!tiny_quadratic_root(D2_H1, D2_H2, y, bits) ||
        !tiny_quadratic_root(D2_L1, D2_L2, v, bits)) {
        return 0;
    }

    mpf_mul_ui(z, y, 1728UL);
    mpf_ui_sub(tmp, 1UL, z);
    mpf_mul(tmp, tmp, v);
    mpf_mul_ui(tmp, tmp, 427UL);
    mpf_div(alpha, y, tmp);
    mpf_ui_sub(beta, 1UL, alpha);

    mpf_set_ui(K, 427UL);
    mpf_sqrt(K, K);
    mpf_div_ui(K, K, 6UL);
    mpf_ui_sub(tmp, 1UL, z);
    mpf_sqrt(tmp, tmp);
    mpf_mul(K, K, tmp);

    const long double digits_per_term =
        24.95589965765426673091470594098813130578L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 20.0L) /
                            digits_per_term) + 2UL;

    mpf_set_ui(term, 1UL);
    mpf_set_ui(F, 1UL);
    mpf_set_ui(T, 0UL);
    for (unsigned long n = 1; n < terms; ++n) {
        hyper_term_advance(term, z, n);
        mpf_add(F, F, term);
        mpf_mul_ui(tmp, term, n);
        mpf_add(T, T, tmp);
    }

    mpf_mul(inv_pi, beta, F);
    mpf_mul_ui(tmp, T, 6UL);
    mpf_add(inv_pi, inv_pi, tmp);
    mpf_mul(inv_pi, inv_pi, K);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;

    mpf_clears(y, v, z, alpha, beta, K, term, F, T, tmp, inv_pi, pi, NULL);
    return ok;
}

/* Class-number-1 D=163 direct evaluator.  Algebraically this is Chudnovsky,
 * but it intentionally uses the same direct hypergeometric recurrence as D2
 * so that convergence-rate benefit can be separated from binary splitting.
 */
static int chud_direct_pi(uint64_t digits, uint64_t guard, const char *path,
                          double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    mpf_t z, alpha, beta, K, term, F, T, tmp, den, inv_pi, pi;
    mpf_init2(z, bits);      mpf_init2(alpha, bits);
    mpf_init2(beta, bits);   mpf_init2(K, bits);
    mpf_init2(term, bits);   mpf_init2(F, bits);
    mpf_init2(T, bits);      mpf_init2(tmp, bits);
    mpf_init2(den, bits);    mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    mpf_set_si(z, -1);
    mpf_set_ui(den, 53360UL);
    mpf_pow_ui(den, den, 3UL);
    mpf_div(z, z, den);
    mpf_set_ui(alpha, 77265280UL);
    mpf_div_ui(alpha, alpha, 90856689UL);
    mpf_ui_sub(beta, 1UL, alpha);

    mpf_set_ui(K, 163UL);
    mpf_sqrt(K, K);
    mpf_div_ui(K, K, 6UL);
    mpf_ui_sub(tmp, 1UL, z);
    mpf_sqrt(tmp, tmp);
    mpf_mul(K, K, tmp);

    const long double digits_per_term = 14.1816474627254776555255216782L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 20.0L) /
                            digits_per_term) + 2UL;

    mpf_set_ui(term, 1UL);
    mpf_set_ui(F, 1UL);
    mpf_set_ui(T, 0UL);
    for (unsigned long n = 1; n < terms; ++n) {
        hyper_term_advance(term, z, n);
        mpf_add(F, F, term);
        mpf_mul_ui(tmp, term, n);
        mpf_add(T, T, tmp);
    }

    mpf_mul(inv_pi, beta, F);
    mpf_mul_ui(tmp, T, 6UL);
    mpf_add(inv_pi, inv_pi, tmp);
    mpf_mul(inv_pi, inv_pi, K);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;

    mpf_clears(z, alpha, beta, K, term, F, T, tmp, den, inv_pi, pi, NULL);
    return ok;
}

static void bs_chud(unsigned long a, unsigned long b,
                    mpz_t P, mpz_t Q, mpz_t T) {
    if (b - a == 1UL) {
        if (a == 0UL) {
            mpz_set_ui(P, 1UL);
            mpz_set_ui(Q, 1UL);
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
    mpz_inits(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
    bs_chud(a, m, P1, Q1, T1);
    bs_chud(m, b, P2, Q2, T2);
    mpz_mul(P, P1, P2);
    mpz_mul(Q, Q1, Q2);
    mpz_mul(T, T1, Q2);
    mpz_mul(tmp, P1, T2);
    mpz_add(T, T, tmp);
    mpz_clears(P1, Q1, T1, P2, Q2, T2, tmp, NULL);
}

static int chud_bs_pi(uint64_t digits, uint64_t guard, const char *path,
                      double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    const long double digits_per_term = 14.1816474627254776555L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard) /
                            digits_per_term) + 2UL;

    mpz_t P, Q, T;
    mpz_inits(P, Q, T, NULL);
    bs_chud(0UL, terms, P, Q, T);

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

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;

    mpf_clears(fq, ft, root, pi, NULL);
    mpz_clears(P, Q, T, NULL);
    return ok;
}

static int run_case(uint64_t digits) {
    const uint64_t guard = 160ULL;
    char bs_path[256], direct_path[256], d2_path[256];
    snprintf(bs_path, sizeof(bs_path), "build/d2_chud_bs_%llu.txt",
             (unsigned long long)digits);
    snprintf(direct_path, sizeof(direct_path), "build/d2_chud_direct_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/d2_basic_%llu.txt",
             (unsigned long long)digits);

    double tbs = 0.0, tdirect = 0.0, td2 = 0.0;
    unsigned long nbs = 0UL, ndirect = 0UL, nd2 = 0UL;
    if (!chud_bs_pi(digits, guard, bs_path, &tbs, &nbs) ||
        !chud_direct_pi(digits, guard, direct_path, &tdirect, &ndirect) ||
        !d2_basic_pi(digits, guard, d2_path, &td2, &nd2)) {
        return 0;
    }

    int eq_bs = files_equal(bs_path, d2_path);
    int eq_direct = files_equal(direct_path, d2_path);
    printf("%llu,%lu,%.9f,%lu,%.9f,%lu,%.9f,%.6f,%.6f,%d,%d\n",
           (unsigned long long)digits,
           nbs, tbs,
           ndirect, tdirect,
           nd2, td2,
           td2 / tbs,
           td2 / tdirect,
           eq_bs,
           eq_direct);
    fflush(stdout);
    return eq_bs && eq_direct;
}

int main(int argc, char **argv) {
    printf("digits,chud_bs_terms,chud_bs_seconds,chud_direct_terms,chud_direct_seconds,d2_terms,d2_seconds,d2_over_chud_bs,d2_over_chud_direct,equal_bs,equal_direct\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {1000ULL, 10000ULL, 30000ULL, 100000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i) {
            if (!run_case(defaults[i])) return 1;
        }
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !run_case(digits)) return 1;
    }
    return 0;
}
