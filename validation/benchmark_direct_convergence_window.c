#define main original_d2_direct_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/* Convergence-Tapered Direct Recurrence (CTDR).
 *
 * This stays on the original scalar/direct evaluator.  No binary splitting
 * and no quadratic-field product tree are used.  The only acceleration is
 * precision tapering: term n is about 10^(-q n), so its relative precision
 * can fall with n while preserving a fixed absolute error target.  Precision
 * is changed only at coarse block boundaries to avoid allocator churn.
 */
static mp_bitcnt_t ct_active_bits(mp_bitcnt_t max_bits,
                                  long double digits_per_term,
                                  unsigned long n) {
    long double dropped = digits_per_term * (long double)n * 3.32192809488736234787L;
    long double keep = (long double)max_bits - dropped + 1024.0L;
    if (keep < 2048.0L) keep = 2048.0L;
    if (keep > (long double)max_bits) keep = (long double)max_bits;
    return (mp_bitcnt_t)ceill(keep);
}

static void ct_series(const mpf_t z,
                      long double digits_per_term,
                      unsigned long terms,
                      mp_bitcnt_t bits,
                      unsigned long block_terms,
                      mpf_t F,
                      mpf_t T) {
    mpf_t term, zwork, nt;
    mpf_init2(term, bits);
    mpf_init2(zwork, bits);
    mpf_init2(nt, bits);
    mpf_set_ui(term, 1UL);
    mpf_set_ui(F, 1UL);
    mpf_set_ui(T, 0UL);

    for (unsigned long a = 1UL; a < terms; a += block_terms) {
        unsigned long b = a + block_terms;
        if (b > terms) b = terms;
        mp_bitcnt_t abits = ct_active_bits(bits, digits_per_term, a);

        /* Safe precision changes at block boundaries. */
        mpf_set_prec(term, abits);
        mpf_set_prec(zwork, abits);
        mpf_set_prec(nt, abits);
        mpf_set(zwork, z);

        for (unsigned long n = a; n < b; ++n) {
            hyper_term_advance(term, zwork, n);
            mpf_add(F, F, term);
            mpf_mul_ui(nt, term, n);
            mpf_add(T, T, nt);
        }
    }

    mpf_clears(term, zwork, nt, NULL);
}

static int d2_ct_pi(uint64_t digits, uint64_t guard, unsigned long block_terms,
                    const char *path, double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    mpf_t y, v, z, alpha, beta, K, F, T, tmp, inv_pi, pi;
    mpf_init2(y, bits);      mpf_init2(v, bits);
    mpf_init2(z, bits);      mpf_init2(alpha, bits);
    mpf_init2(beta, bits);   mpf_init2(K, bits);
    mpf_init2(F, bits);      mpf_init2(T, bits);
    mpf_init2(tmp, bits);    mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    if (!tiny_quadratic_root(D2_H1, D2_H2, y, bits) ||
        !tiny_quadratic_root(D2_L1, D2_L2, v, bits)) return 0;

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

    const long double q = 24.95589965765426673091470594098813130578L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 20.0L) / q) + 2UL;
    ct_series(z, q, terms, bits, block_terms, F, T);

    mpf_mul(inv_pi, beta, F);
    mpf_mul_ui(tmp, T, 6UL);
    mpf_add(inv_pi, inv_pi, tmp);
    mpf_mul(inv_pi, inv_pi, K);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;
    mpf_clears(y, v, z, alpha, beta, K, F, T, tmp, inv_pi, pi, NULL);
    return ok;
}

static int chud_ct_pi(uint64_t digits, uint64_t guard, unsigned long block_terms,
                      const char *path, double *seconds, unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    mpf_t z, alpha, beta, K, F, T, tmp, den, inv_pi, pi;
    mpf_init2(z, bits);      mpf_init2(alpha, bits);
    mpf_init2(beta, bits);   mpf_init2(K, bits);
    mpf_init2(F, bits);      mpf_init2(T, bits);
    mpf_init2(tmp, bits);    mpf_init2(den, bits);
    mpf_init2(inv_pi, bits); mpf_init2(pi, bits);

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

    const long double q = 14.1816474627254776555255216782L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 20.0L) / q) + 2UL;
    ct_series(z, q, terms, bits, block_terms, F, T);

    mpf_mul(inv_pi, beta, F);
    mpf_mul_ui(tmp, T, 6UL);
    mpf_add(inv_pi, inv_pi, tmp);
    mpf_mul(inv_pi, inv_pi, K);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;
    mpf_clears(z, alpha, beta, K, F, T, tmp, den, inv_pi, pi, NULL);
    return ok;
}

static int run_ct(uint64_t digits, unsigned long block_terms) {
    const uint64_t guard = 192ULL;
    char ref[256], db[256], dct[256], cb[256], cct[256];
    snprintf(ref, sizeof(ref), "build/ct_ref_%llu_%lu.txt", (unsigned long long)digits, block_terms);
    snprintf(db, sizeof(db), "build/ct_d2_base_%llu_%lu.txt", (unsigned long long)digits, block_terms);
    snprintf(dct, sizeof(dct), "build/ct_d2_%llu_%lu.txt", (unsigned long long)digits, block_terms);
    snprintf(cb, sizeof(cb), "build/ct_chud_base_%llu_%lu.txt", (unsigned long long)digits, block_terms);
    snprintf(cct, sizeof(cct), "build/ct_chud_%llu_%lu.txt", (unsigned long long)digits, block_terms);

    double tr = 0.0, tdb = 0.0, tdct = 0.0, tcb = 0.0, tcct = 0.0;
    unsigned long nr = 0, ndb = 0, ndct = 0, ncb = 0, ncct = 0;
    if (!chud_bs_pi(digits, guard, ref, &tr, &nr) ||
        !d2_basic_pi(digits, guard, db, &tdb, &ndb) ||
        !d2_ct_pi(digits, guard, block_terms, dct, &tdct, &ndct) ||
        !chud_direct_pi(digits, guard, cb, &tcb, &ncb) ||
        !chud_ct_pi(digits, guard, block_terms, cct, &tcct, &ncct)) return 0;

    int e1 = files_equal(ref, db);
    int e2 = files_equal(ref, dct);
    int e3 = files_equal(ref, cb);
    int e4 = files_equal(ref, cct);
    printf("%llu,%lu,%.9f,%.9f,%.9f,%.9f,%.9f,%.6f,%.6f,%d,%d,%d,%d\n",
           (unsigned long long)digits, block_terms,
           tr, tdb, tdct, tcb, tcct,
           tdct / tdb, tdct / tcct,
           e1, e2, e3, e4);
    fflush(stdout);
    return e1 && e2 && e3 && e4;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s DIGITS BLOCK_TERMS\n", argv[0]);
        return 2;
    }
    uint64_t digits = strtoull(argv[1], NULL, 10);
    unsigned long block = strtoul(argv[2], NULL, 10);
    printf("digits,block_terms,chud_bs_s,d2_direct_s,d2_ct_s,chud_direct_s,chud_ct_s,d2_ct_over_d2_direct,d2_ct_over_chud_ct,eq_d2_base,eq_d2_ct,eq_chud_base,eq_chud_ct\n");
    return (digits && block && run_ct(digits, block)) ? 0 : 1;
}
