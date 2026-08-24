#define main d2_bs_serial_main
#include "benchmark_d2_bs_vs_classics.c"
#undef main

#include <omp.h>

/*
 * Exact two-lane D2 execution.
 *
 * The quadratic basis has two output coordinates, a and b. A product
 *
 *   (a+b*w)(c+d*w), w^2=w+15
 *
 * is evaluated as two coordinate jobs:
 *
 *   lane A: ac + 15 bd
 *   lane B: (a+b)(c+d) - ac
 *
 * For sufficiently large operands, lane A and lane B run on two OpenMP
 * threads. The duplicated ac product is intentional: it removes an
 * inter-lane dependency and trades one extra multiplication for lower wall
 * time. Small products keep the 3-multiply serial kernel because scheduling
 * overhead dominates there.
 *
 * Arithmetic remains exact mpz arithmetic and the product-tree topology is
 * unchanged.
 */

static mp_bitcnt_t d2_parallel_threshold_bits = 65536UL;

static mp_bitcnt_t d2_ring_max_bits(const d2_ring_t *x, const d2_ring_t *y) {
    mp_bitcnt_t m = 0, v;
    v = mpz_sizeinbase(x->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(x->b, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->b, 2); if (v > m) m = v;
    return m;
}

static void d2_ring_mul_2lane(d2_ring_t *out,
                              const d2_ring_t *x,
                              const d2_ring_t *y,
                              d2_ws_t *ws) {
    if (d2_ring_max_bits(x, y) < d2_parallel_threshold_bits) {
        d2_ring_mul(out, x, y, ws);
        return;
    }

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
        {
            mpz_mul(ws->m0, x->a, y->a);
            mpz_mul(ws->m1, x->b, y->b);
            mpz_set(ws->t0, ws->m0);
            mpz_addmul_ui(ws->t0, ws->m1, 15UL);
        }
#pragma omp section
        {
            mpz_add(ws->num, x->a, x->b);
            mpz_add(ws->den, y->a, y->b);
            mpz_mul(ws->m2, ws->num, ws->den);
            mpz_mul(ws->t1, x->a, y->a);
            mpz_sub(ws->m2, ws->m2, ws->t1);
        }
    }

    mpz_set(out->a, ws->t0);
    mpz_set(out->b, ws->m2);
}

static void d2_bs_build_2lane(d2_seg_t *out,
                              unsigned long a,
                              unsigned long b,
                              const mpz_t z0,
                              const mpz_t z1,
                              const mpz_t zd,
                              const mpz_t b0,
                              const mpz_t b1,
                              const mpz_t bd,
                              d2_ws_t *ws) {
    if (b - a == 1UL) {
        d2_leaf(out, a, z0, z1, zd, b0, b1, bd, ws);
        return;
    }

    unsigned long m = a + (b - a) / 2UL;
    d2_seg_t left, right;
    d2_seg_init(&left);
    d2_seg_init(&right);
    d2_bs_build_2lane(&left, a, m, z0, z1, zd, b0, b1, bd, ws);
    d2_bs_build_2lane(&right, m, b, z0, z1, zd, b0, b1, bd, ws);

    mpz_mul(out->q, left.q, right.q);
    d2_ring_mul_2lane(&out->p, &left.p, &right.p, ws);

    /* W = W_L + P_L*W_R. */
    d2_ring_mul_2lane(&ws->rw, &left.p, &right.w, ws);

    if (mpz_sizeinbase(right.q, 2) >= d2_parallel_threshold_bits) {
#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            {
                mpz_mul(out->w.a, left.w.a, right.q);
                mpz_add(out->w.a, out->w.a, ws->rw.a);
            }
#pragma omp section
            {
                mpz_mul(out->w.b, left.w.b, right.q);
                mpz_add(out->w.b, out->w.b, ws->rw.b);
            }
        }
    } else {
        mpz_mul(out->w.a, left.w.a, right.q);
        mpz_add(out->w.a, out->w.a, ws->rw.a);
        mpz_mul(out->w.b, left.w.b, right.q);
        mpz_add(out->w.b, out->w.b, ws->rw.b);
    }

    d2_seg_clear(&left);
    d2_seg_clear(&right);
}

static int d2_bs_pi_2lane(uint64_t digits,
                          uint64_t guard,
                          const char *path,
                          double *seconds,
                          unsigned long *terms_out) {
    double t0 = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);

    const long double digits_per_term =
        24.95589965765426673091470594098813130578L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 30.0L) /
                            digits_per_term) + 2UL;

    mpz_t z0, z1, zd, b0, b1, bd;
    mpz_inits(z0, z1, zd, b0, b1, bd, NULL);
    if (mpz_set_str(z0, D2_Z0, 10) != 0 ||
        mpz_set_str(z1, D2_Z1, 10) != 0 ||
        mpz_set_str(zd, D2_ZD, 10) != 0 ||
        mpz_set_str(b0, D2_B0, 10) != 0 ||
        mpz_set_str(b1, D2_B1, 10) != 0 ||
        mpz_set_str(bd, D2_BD, 10) != 0) {
        mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
        return 0;
    }

    d2_seg_t tree;
    d2_seg_init(&tree);
    d2_ws_t ws;
    d2_ws_init(&ws);
    d2_bs_build_2lane(&tree, 0UL, terms, z0, z1, zd, b0, b1, bd, &ws);

    mpf_t omega, root61, fz0, fz1, fzd, fw0, fw1, fq, fbd,
          z, W, K, tmp, inv_pi, pi;
    mpf_init2(omega, bits); mpf_init2(root61, bits);
    mpf_init2(fz0, bits); mpf_init2(fz1, bits);
    mpf_init2(fzd, bits); mpf_init2(fw0, bits);
    mpf_init2(fw1, bits); mpf_init2(fq, bits);
    mpf_init2(fbd, bits); mpf_init2(z, bits);
    mpf_init2(W, bits); mpf_init2(K, bits);
    mpf_init2(tmp, bits); mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    mpf_set_ui(root61, 61UL);
    mpf_sqrt(root61, root61);
    mpf_add_ui(omega, root61, 1UL);
    mpf_div_ui(omega, omega, 2UL);

    mpf_set_z(fz0, z0);
    mpf_set_z(fz1, z1);
    mpf_mul(fz1, fz1, omega);
    mpf_add(z, fz0, fz1);
    mpf_set_z(fzd, zd);
    mpf_div(z, z, fzd);

    mpf_set_z(fw0, tree.w.a);
    mpf_set_z(fw1, tree.w.b);
    mpf_mul(fw1, fw1, omega);
    mpf_add(W, fw0, fw1);
    mpf_set_z(fq, tree.q);
    mpf_div(W, W, fq);
    mpf_set_z(fbd, bd);
    mpf_div(W, W, fbd);

    mpf_set_ui(K, 427UL);
    mpf_sqrt(K, K);
    mpf_div_ui(K, K, 6UL);
    mpf_ui_sub(tmp, 1UL, z);
    mpf_sqrt(tmp, tmp);
    mpf_mul(K, K, tmp);

    mpf_mul(inv_pi, K, W);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0;
    if (terms_out) *terms_out = terms;

    mpf_clears(omega, root61, fz0, fz1, fzd, fw0, fw1, fq, fbd,
               z, W, K, tmp, inv_pi, pi, NULL);
    d2_ws_clear(&ws);
    d2_seg_clear(&tree);
    mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
    return ok;
}

static int run_parallel_case(uint64_t digits) {
    const uint64_t guard = 192ULL;
    char ref_path[256], d2_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/d2parallel_ref_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/d2parallel_%llu.txt",
             (unsigned long long)digits);

    double tc = 0.0, td = 0.0;
    unsigned long nc = 0UL, nd = 0UL;
    if (!chud_bs_pi(digits, guard, ref_path, &tc, &nc) ||
        !d2_bs_pi_2lane(digits, guard, d2_path, &td, &nd)) return 0;
    int eq = files_equal(ref_path, d2_path);

    printf("%llu,%.9f,%lu,%.9f,%lu,%.6f,%lu,%d\n",
           (unsigned long long)digits,
           tc, nc, td, nd, td / tc,
           (unsigned long)d2_parallel_threshold_bits, eq);
    fflush(stdout);
    return eq;
}

int main(int argc, char **argv) {
    const char *e = getenv("D2_PARALLEL_THRESHOLD_BITS");
    if (e && *e) {
        unsigned long long v = strtoull(e, NULL, 10);
        if (v > 0ULL) d2_parallel_threshold_bits = (mp_bitcnt_t)v;
    }

    omp_set_dynamic(0);
    omp_set_num_threads(2);

    printf("digits,chud_bs_seconds,chud_terms,d2_2lane_seconds,d2_terms,d2_over_chud,parallel_threshold_bits,equal\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {100000ULL, 300000ULL, 1000000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i)
            if (!run_parallel_case(defaults[i])) return 1;
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !run_parallel_case(digits)) return 1;
    }
    return 0;
}
