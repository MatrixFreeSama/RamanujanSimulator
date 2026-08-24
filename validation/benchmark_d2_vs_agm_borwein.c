#define main original_d2_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/*
 * End-to-end comparison of the original class-number-2 direct evaluator with
 * two classical iterative pi algorithms. All methods use the same GMP mpf
 * precision budget and decimal serializer, and every output is byte-compared
 * with Chudnovsky binary splitting.
 */

static unsigned long ceil_log2_u64(uint64_t n) {
    unsigned long k = 0;
    uint64_t x = 1;
    while (x < n && k < 63UL) {
        x <<= 1;
        ++k;
    }
    return k;
}

static int agm_pi(uint64_t digits, uint64_t guard, const char *path,
                  double *seconds, unsigned long *iterations_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    mpf_t a, b, t, p, an, bn, delta, tmp, sum, pi;
    mpf_init2(a, bits);     mpf_init2(b, bits);
    mpf_init2(t, bits);     mpf_init2(p, bits);
    mpf_init2(an, bits);    mpf_init2(bn, bits);
    mpf_init2(delta, bits); mpf_init2(tmp, bits);
    mpf_init2(sum, bits);   mpf_init2(pi, bits);

    mpf_set_ui(a, 1UL);
    mpf_set_ui(tmp, 2UL);
    mpf_sqrt(tmp, tmp);
    mpf_ui_div(b, 1UL, tmp);
    mpf_set_ui(t, 1UL);
    mpf_div_ui(t, t, 4UL);
    mpf_set_ui(p, 1UL);

    unsigned long iterations = ceil_log2_u64(digits + guard + 1ULL) + 4UL;
    for (unsigned long i = 0; i < iterations; ++i) {
        mpf_add(an, a, b);
        mpf_div_ui(an, an, 2UL);

        mpf_mul(tmp, a, b);
        mpf_sqrt(bn, tmp);

        mpf_sub(delta, a, an);
        mpf_mul(delta, delta, delta);
        mpf_mul(tmp, p, delta);
        mpf_sub(t, t, tmp);

        mpf_mul_ui(p, p, 2UL);
        mpf_set(a, an);
        mpf_set(b, bn);
    }

    mpf_add(sum, a, b);
    mpf_mul(sum, sum, sum);
    mpf_mul_ui(tmp, t, 4UL);
    mpf_div(pi, sum, tmp);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (iterations_out) *iterations_out = iterations;

    mpf_clears(a, b, t, p, an, bn, delta, tmp, sum, pi, NULL);
    return ok;
}

static int borwein_quartic_pi(uint64_t digits, uint64_t guard, const char *path,
                              double *seconds, unsigned long *iterations_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    mpf_t sqrt2, y, yn, a, y2, y4, root, num, den, one,
          oneplus, oneplus2, oneplus4, poly, tmp, factor, pi;
    mpf_init2(sqrt2, bits);   mpf_init2(y, bits);
    mpf_init2(yn, bits);      mpf_init2(a, bits);
    mpf_init2(y2, bits);      mpf_init2(y4, bits);
    mpf_init2(root, bits);    mpf_init2(num, bits);
    mpf_init2(den, bits);     mpf_init2(one, bits);
    mpf_init2(oneplus, bits); mpf_init2(oneplus2, bits);
    mpf_init2(oneplus4, bits);mpf_init2(poly, bits);
    mpf_init2(tmp, bits);     mpf_init2(factor, bits);
    mpf_init2(pi, bits);

    mpf_set_ui(one, 1UL);
    mpf_set_ui(sqrt2, 2UL);
    mpf_sqrt(sqrt2, sqrt2);
    mpf_sub_ui(y, sqrt2, 1UL);

    mpf_mul_ui(a, sqrt2, 4UL);
    mpf_ui_sub(a, 6UL, a);

    unsigned long log2n = ceil_log2_u64(digits + guard + 1ULL);
    unsigned long iterations = (log2n + 1UL) / 2UL + 3UL;
    unsigned long factor_ui = 8UL;

    for (unsigned long i = 0; i < iterations; ++i) {
        mpf_mul(y2, y, y);
        mpf_mul(y4, y2, y2);
        mpf_ui_sub(root, 1UL, y4);
        mpf_sqrt(root, root);
        mpf_sqrt(root, root);

        mpf_ui_sub(num, 1UL, root);
        mpf_add_ui(den, root, 1UL);
        mpf_div(yn, num, den);

        mpf_add_ui(oneplus, yn, 1UL);
        mpf_mul(oneplus2, oneplus, oneplus);
        mpf_mul(oneplus4, oneplus2, oneplus2);

        mpf_mul(y2, yn, yn);
        mpf_add(poly, one, yn);
        mpf_add(poly, poly, y2);
        mpf_mul(poly, poly, yn);

        mpf_set_ui(factor, factor_ui);
        mpf_mul(poly, poly, factor);
        mpf_mul(tmp, a, oneplus4);
        mpf_sub(a, tmp, poly);

        mpf_set(y, yn);
        if (factor_ui <= ULONG_MAX / 4UL) factor_ui *= 4UL;
    }

    mpf_ui_div(pi, 1UL, a);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (iterations_out) *iterations_out = iterations;

    mpf_clears(sqrt2, y, yn, a, y2, y4, root, num, den, one,
               oneplus, oneplus2, oneplus4, poly, tmp, factor, pi, NULL);
    return ok;
}

static int run_compare_case(uint64_t digits) {
    const uint64_t guard = 192ULL;
    char ref_path[256], d2_path[256], agm_path[256], bq_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/classic_ref_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/classic_d2_%llu.txt",
             (unsigned long long)digits);
    snprintf(agm_path, sizeof(agm_path), "build/classic_agm_%llu.txt",
             (unsigned long long)digits);
    snprintf(bq_path, sizeof(bq_path), "build/classic_bq_%llu.txt",
             (unsigned long long)digits);

    double tref = 0.0, td2 = 0.0, tagm = 0.0, tbq = 0.0;
    unsigned long nref = 0UL, nd2 = 0UL, iagm = 0UL, ibq = 0UL;

    if (!chud_bs_pi(digits, guard, ref_path, &tref, &nref) ||
        !d2_basic_pi(digits, guard, d2_path, &td2, &nd2) ||
        !agm_pi(digits, guard, agm_path, &tagm, &iagm) ||
        !borwein_quartic_pi(digits, guard, bq_path, &tbq, &ibq)) {
        return 0;
    }

    int eq_d2 = files_equal(ref_path, d2_path);
    int eq_agm = files_equal(ref_path, agm_path);
    int eq_bq = files_equal(ref_path, bq_path);

    printf("%llu,%.9f,%lu,%.9f,%lu,%.9f,%lu,%.9f,%lu,%.6f,%.6f,%d,%d,%d\n",
           (unsigned long long)digits,
           tref, nref,
           td2, nd2,
           tagm, iagm,
           tbq, ibq,
           td2 / tagm,
           td2 / tbq,
           eq_d2, eq_agm, eq_bq);
    fflush(stdout);
    return eq_d2 && eq_agm && eq_bq;
}

int main(int argc, char **argv) {
    printf("digits,chud_bs_seconds,chud_terms,d2_seconds,d2_terms,agm_seconds,agm_iterations,borwein4_seconds,borwein4_iterations,d2_over_agm,d2_over_borwein4,equal_d2,equal_agm,equal_borwein4\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {1000ULL, 10000ULL, 30000ULL, 100000ULL, 200000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i) {
            if (!run_compare_case(defaults[i])) return 1;
        }
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !run_compare_case(digits)) return 1;
    }
    return 0;
}
