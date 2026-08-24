#define main original_d2_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

#include <omp.h>

/* Exact two-lane D2 execution in the integral basis
 *   omega = (1 + sqrt(61))/2, omega^2 = omega + 15.
 *
 * A ring product
 *   (a+b*omega)(c+d*omega)
 * has two output coordinates. For large operands they are assigned to two
 * independent OpenMP workers:
 *
 *   lane A: ac + 15 bd
 *   lane B: (a+b)(c+d) - ac
 *
 * The ac product is deliberately duplicated so neither output lane waits for
 * the other. This is exact integer arithmetic, not a floating conjugate-lane
 * approximation. Small operands retain the 3-product serial kernel.
 */

static const char *P2_Z0 = "-59818419102592333";
static const char *P2_Z1 = "13579278976889262";
static const char *P2_ZD = "609541351191872000";
static const char *P2_B0 = "320004750671";
static const char *P2_B1 = "-56552805760";
static const char *P2_BD = "766923569391";

static mp_bitcnt_t p2_threshold_bits = 65536UL;

typedef struct {
    mpz_t a;
    mpz_t b;
} p2_ring_t;

typedef struct {
    mpz_t q;
    p2_ring_t p;
    p2_ring_t w;
} p2_seg_t;

typedef struct {
    mpz_t num, den, g, weight0;
    mpz_t m0, m1, m2, t0, t1;
    p2_ring_t rw;
} p2_ws_t;

static void p2_ring_init(p2_ring_t *r) {
    mpz_init(r->a);
    mpz_init(r->b);
}

static void p2_ring_clear(p2_ring_t *r) {
    mpz_clear(r->a);
    mpz_clear(r->b);
}

static void p2_seg_init(p2_seg_t *s) {
    mpz_init(s->q);
    p2_ring_init(&s->p);
    p2_ring_init(&s->w);
}

static void p2_seg_clear(p2_seg_t *s) {
    mpz_clear(s->q);
    p2_ring_clear(&s->p);
    p2_ring_clear(&s->w);
}

static void p2_ws_init(p2_ws_t *ws) {
    mpz_inits(ws->num, ws->den, ws->g, ws->weight0,
              ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
    p2_ring_init(&ws->rw);
}

static void p2_ws_clear(p2_ws_t *ws) {
    p2_ring_clear(&ws->rw);
    mpz_clears(ws->num, ws->den, ws->g, ws->weight0,
               ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
}

static mp_bitcnt_t p2_max_bits(const p2_ring_t *x, const p2_ring_t *y) {
    mp_bitcnt_t m = 0, v;
    v = mpz_sizeinbase(x->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(x->b, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->b, 2); if (v > m) m = v;
    return m;
}

static void p2_ring_mul_serial(p2_ring_t *out,
                               const p2_ring_t *x,
                               const p2_ring_t *y,
                               p2_ws_t *ws) {
    mpz_mul(ws->m0, x->a, y->a);
    mpz_mul(ws->m1, x->b, y->b);
    mpz_add(ws->t0, x->a, x->b);
    mpz_add(ws->t1, y->a, y->b);
    mpz_mul(ws->m2, ws->t0, ws->t1);
    mpz_set(out->a, ws->m0);
    mpz_addmul_ui(out->a, ws->m1, 15UL);
    mpz_sub(out->b, ws->m2, ws->m0);
}

static void p2_ring_mul_2lane(p2_ring_t *out,
                              const p2_ring_t *x,
                              const p2_ring_t *y,
                              p2_ws_t *ws) {
    if (p2_max_bits(x, y) < p2_threshold_bits) {
        p2_ring_mul_serial(out, x, y, ws);
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

static void p2_leaf(p2_seg_t *out,
                    unsigned long n,
                    const mpz_t z0,
                    const mpz_t z1,
                    const mpz_t zd,
                    const mpz_t b0,
                    const mpz_t b1,
                    const mpz_t bd,
                    p2_ws_t *ws) {
    mpz_set_ui(ws->num, 6UL * n + 1UL);
    mpz_mul_ui(ws->num, ws->num, 2UL * n + 1UL);
    mpz_mul_ui(ws->num, ws->num, 6UL * n + 5UL);

    mpz_set_ui(ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, n + 1UL);
    mpz_mul_ui(ws->den, ws->den, 72UL);
    mpz_mul(ws->den, ws->den, zd);

    mpz_gcd(ws->g, ws->num, ws->den);
    mpz_divexact(ws->num, ws->num, ws->g);
    mpz_divexact(out->q, ws->den, ws->g);

    mpz_mul(out->p.a, z0, ws->num);
    mpz_mul(out->p.b, z1, ws->num);

    mpz_mul_ui(ws->weight0, bd, 6UL * n);
    mpz_add(ws->weight0, ws->weight0, b0);
    mpz_mul(out->w.a, ws->weight0, out->q);
    mpz_mul(out->w.b, b1, out->q);
}

static void p2_build(p2_seg_t *out,
                     unsigned long a,
                     unsigned long b,
                     const mpz_t z0,
                     const mpz_t z1,
                     const mpz_t zd,
                     const mpz_t b0,
                     const mpz_t b1,
                     const mpz_t bd,
                     p2_ws_t *ws) {
    if (b - a == 1UL) {
        p2_leaf(out, a, z0, z1, zd, b0, b1, bd, ws);
        return;
    }

    unsigned long m = a + (b - a) / 2UL;
    p2_seg_t left, right;
    p2_seg_init(&left);
    p2_seg_init(&right);
    p2_build(&left, a, m, z0, z1, zd, b0, b1, bd, ws);
    p2_build(&right, m, b, z0, z1, zd, b0, b1, bd, ws);

    mpz_mul(out->q, left.q, right.q);
    p2_ring_mul_2lane(&out->p, &left.p, &right.p, ws);
    p2_ring_mul_2lane(&ws->rw, &left.p, &right.w, ws);

    if (mpz_sizeinbase(right.q, 2) >= p2_threshold_bits) {
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

    p2_seg_clear(&left);
    p2_seg_clear(&right);
}

static int p2_pi(uint64_t digits,
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
    if (mpz_set_str(z0, P2_Z0, 10) != 0 ||
        mpz_set_str(z1, P2_Z1, 10) != 0 ||
        mpz_set_str(zd, P2_ZD, 10) != 0 ||
        mpz_set_str(b0, P2_B0, 10) != 0 ||
        mpz_set_str(b1, P2_B1, 10) != 0 ||
        mpz_set_str(bd, P2_BD, 10) != 0) {
        mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
        return 0;
    }

    p2_seg_t tree;
    p2_seg_init(&tree);
    p2_ws_t ws;
    p2_ws_init(&ws);
    p2_build(&tree, 0UL, terms, z0, z1, zd, b0, b1, bd, &ws);

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
    p2_ws_clear(&ws);
    p2_seg_clear(&tree);
    mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
    return ok;
}

static int run_case_2lane(uint64_t digits) {
    const uint64_t guard = 192ULL;
    char ref_path[256], d2_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/d2parallel_ref_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/d2parallel_%llu.txt",
             (unsigned long long)digits);

    double tc = 0.0, td = 0.0;
    unsigned long nc = 0UL, nd = 0UL;
    if (!chud_bs_pi(digits, guard, ref_path, &tc, &nc) ||
        !p2_pi(digits, guard, d2_path, &td, &nd)) return 0;
    int eq = files_equal(ref_path, d2_path);
    printf("%llu,%.9f,%lu,%.9f,%lu,%.6f,%lu,%d\n",
           (unsigned long long)digits,
           tc, nc, td, nd, td / tc,
           (unsigned long)p2_threshold_bits, eq);
    fflush(stdout);
    return eq;
}

int main(int argc, char **argv) {
    const char *e = getenv("D2_PARALLEL_THRESHOLD_BITS");
    if (e && *e) {
        unsigned long long v = strtoull(e, NULL, 10);
        if (v > 0ULL) p2_threshold_bits = (mp_bitcnt_t)v;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(2);

    printf("digits,chud_bs_seconds,chud_terms,d2_2lane_seconds,d2_terms,d2_over_chud,parallel_threshold_bits,equal\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {100000ULL, 300000ULL, 1000000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i)
            if (!run_case_2lane(defaults[i])) return 1;
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !run_case_2lane(digits)) return 1;
    }
    return 0;
}
