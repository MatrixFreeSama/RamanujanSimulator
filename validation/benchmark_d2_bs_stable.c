#define main d2_bs_single_run_main
#include "benchmark_d2_bs_vs_classics.c"
#undef main

#include <math.h>

static int cmp_double(const void *pa, const void *pb) {
    const double a = *(const double *)pa;
    const double b = *(const double *)pb;
    return (a > b) - (a < b);
}

static double median_copy(const double *x, int n) {
    double tmp[32];
    if (n > (int)(sizeof(tmp) / sizeof(tmp[0]))) return NAN;
    for (int i = 0; i < n; ++i) tmp[i] = x[i];
    qsort(tmp, (size_t)n, sizeof(tmp[0]), cmp_double);
    if (n & 1) return tmp[n / 2];
    return 0.5 * (tmp[n / 2 - 1] + tmp[n / 2]);
}

static double mad(const double *x, int n, double med) {
    double d[32];
    if (n > (int)(sizeof(d) / sizeof(d[0]))) return NAN;
    for (int i = 0; i < n; ++i) d[i] = fabs(x[i] - med);
    return median_copy(d, n);
}

static int stable_case(uint64_t digits, int reps) {
    const uint64_t guard = 192ULL;
    char ref_path[256], d2_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/stable_ref_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/stable_d2_%llu.txt",
             (unsigned long long)digits);

    double junk = 0.0;
    unsigned long nc = 0UL, nd = 0UL;
    if (!chud_bs_pi(digits, guard, ref_path, &junk, &nc) ||
        !d2_bs_pi(digits, guard, d2_path, &junk, &nd) ||
        !files_equal(ref_path, d2_path)) return 0;

    if (!chud_bs_pi(digits, guard, "/dev/null", &junk, &nc) ||
        !d2_bs_pi(digits, guard, "/dev/null", &junk, &nd)) return 0;

    double chud[32], d2[32], ratio[32];
    if (reps > 32) return 0;
    for (int r = 0; r < reps; ++r) {
        double tc = 0.0, td = 0.0;
        unsigned long tc_n = 0UL, td_n = 0UL;
        if ((r & 1) == 0) {
            if (!chud_bs_pi(digits, guard, "/dev/null", &tc, &tc_n) ||
                !d2_bs_pi(digits, guard, "/dev/null", &td, &td_n)) return 0;
        } else {
            if (!d2_bs_pi(digits, guard, "/dev/null", &td, &td_n) ||
                !chud_bs_pi(digits, guard, "/dev/null", &tc, &tc_n)) return 0;
        }
        if (tc_n != nc || td_n != nd) return 0;
        chud[r] = tc;
        d2[r] = td;
        ratio[r] = td / tc;
    }

    const double mc = median_copy(chud, reps);
    const double md = median_copy(d2, reps);
    const double mr = median_copy(ratio, reps);
    const double madc = mad(chud, reps, mc);
    const double madd = mad(d2, reps, md);
    const double madr = mad(ratio, reps, mr);

    printf("%llu,%d,%lu,%lu,%.9f,%.9f,%.9f,%.9f,%.6f,%.6f,%.3f,%.3f\n",
           (unsigned long long)digits, reps, nc, nd,
           mc, madc, md, madd, mr, madr,
           100.0 * madc / mc, 100.0 * madd / md);
    fflush(stdout);
    return 1;
}

int main(int argc, char **argv) {
    const int reps = 9;
    printf("digits,reps,chud_terms,d2_terms,chud_median_s,chud_mad_s,d2_median_s,d2_mad_s,paired_ratio_median,paired_ratio_mad,chud_rel_mad_pct,d2_rel_mad_pct\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {100000ULL, 300000ULL, 1000000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i)
            if (!stable_case(defaults[i], reps)) return 1;
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !stable_case(digits, reps)) return 1;
    }
    return 0;
}
