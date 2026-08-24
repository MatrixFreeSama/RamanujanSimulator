#define main original_d2_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/*
 * Accelerated D2 only.
 *
 * Instead of the large-coefficient x=H2/J basis, use the integral basis
 *   omega = (1 + sqrt(61))/2,  omega^2 = omega + 15.
 *
 * For the selected D=-427 branch the two algebraic quantities needed by the
 * explicit formula become
 *
 *   z = (Z0 + Z1*omega)/ZD
 *   beta = 1-alpha = (B0 + B1*omega)/BD
 *
 * with the small exact coefficients below. This keeps the same binary-
 * splitting/product-tree skeleton as Chudnovsky while avoiding the large
 * H1/H2 coefficients in every quadratic-ring multiplication.
 *
 * The product-tree scratch integers are allocated once per complete solve and
 * reused at every leaf/merge. This removes a large number of GMP init/clear
 * calls without changing the arithmetic or the tree topology.
 */

static const char *D2_Z0 = "-59818419102592333";
static const char *D2_Z1 = "13579278976889262";
static const char *D2_ZD = "609541351191872000";

static const char *D2_B0 = "320004750671";
static const char *D2_B1 = "-56552805760";
static const char *D2_BD = "766923569391";

typedef struct {
    mpz_t a;
    mpz_t b;
} d2_ring_t;

typedef struct {
    mpz_t q;
    d2_ring_t p;
    d2_ring_t w;
} d2_seg_t;

typedef struct {
    mpz_t num;
    mpz_t den;
    mpz_t g;
    mpz_t weight0;
    mpz_t m0;
    mpz_t m1;
    mpz_t m2;
    mpz_t t0;
    mpz_t t1;
    d2_ring_t rw;
} d2_ws_t;

static void d2_ring_init(d2_ring_t *r) {
    mpz_init(r->a);
    mpz_init(r->b);
}

static void d2_ring_clear(d2_ring_t *r) {
    mpz_clear(r->a);
    mpz_clear(r->b);
}

static void d2_seg_init(d2_seg_t *s) {
    mpz_init(s->q);
    d2_ring_init(&s->p);
    d2_ring_init(&s->w);
}

static void d2_seg_clear(d2_seg_t *s) {
    mpz_clear(s->q);
    d2_ring_clear(&s->p);
    d2_ring_clear(&s->w);
}

static void d2_ws_init(d2_ws_t *ws) {
    mpz_inits(ws->num, ws->den, ws->g, ws->weight0,
              ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
    d2_ring_init(&ws->rw);
}

static void d2_ws_clear(d2_ws_t *ws) {
    d2_ring_clear(&ws->rw);
    mpz_clears(ws->num, ws->den, ws->g, ws->weight0,
               ws->m0, ws->m1, ws->m2, ws->t0, ws->t1, NULL);
}

/*
 * (a+b*w)(c+d*w), w^2=w+15.
 * m0=ac, m1=bd, m2=(a+b)(c+d), hence
 *   constant = m0 + 15*m1
 *   w-coeff  = m2 - m0.
 */
static void d2_ring_mul(d2_ring_t *out,
                        const d2_ring_t *x,
                        const d2_ring_t *y,
                        d2_ws_t *ws) {
    mpz_mul(ws->m0, x->a, y->a);
    mpz_mul(ws->m1, x->b, y->b);
    mpz_add(ws->t0, x->a, x->b);
    mpz_add(ws->t1, y->a, y->b);
    mpz_mul(ws->m2, ws->t0, ws->t1);

    mpz_set(out->a, ws->m0);
    mpz_addmul_ui(out->a, ws->m1, 15UL);
    mpz_sub(out->b, ws->m2, ws->m0);
}

static void d2_leaf(d2_seg_t *out,
                    unsigned long n,
                    const mpz_t z0,
                    const mpz_t z1,
                    const mpz_t zd,
                    const mpz_t b0,
                    const mpz_t b1,
                    const mpz_t bd,
                    d2_ws_t *ws) {
    /* r_n = A_n*z / (72*(n+1)^3). */
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

    /* Segment invariant: W = w/(BD*q), with W_leaf=beta+6n. */
    mpz_mul_ui(ws->weight0, bd, 6UL * n);
    mpz_add(ws->weight0, ws->weight0, b0);
    mpz_mul(out->w.a, ws->weight0, out->q);
    mpz_mul(out->w.b, b1, out->q);
}

static void d2_bs_build(d2_seg_t *out,
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
    d2_bs_build(&left, a, m, z0, z1, zd, b0, b1, bd, ws);
    d2_bs_build(&right, m, b, z0, z1, zd, b0, b1, bd, ws);

    mpz_mul(out->q, left.q, right.q);
    d2_ring_mul(&out->p, &left.p, &right.p, ws);

    /* W = W_L + P_L*W_R. */
    d2_ring_mul(&ws->rw, &left.p, &right.w, ws);
    mpz_mul(out->w.a, left.w.a, right.q);
    mpz_add(out->w.a, out->w.a, ws->rw.a);
    mpz_mul(out->w.b, left.w.b, right.q);
    mpz_add(out->w.b, out->w.b, ws->rw.b);

    d2_seg_clear(&left);
    d2_seg_clear(&right);
}

static int d2_bs_pi(uint64_t digits,
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
    d2_bs_build(&tree, 0UL, terms, z0, z1, zd, b0, b1, bd, &ws);

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

    /* z = (Z0 + Z1*omega)/ZD. */
    mpf_set_z(fz0, z0);
    mpf_set_z(fz1, z1);
    mpf_mul(fz1, fz1, omega);
    mpf_add(z, fz0, fz1);
    mpf_set_z(fzd, zd);
    mpf_div(z, z, fzd);

    /* W = (w.a + w.b*omega)/(BD*q). */
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

static int run_d2_bs_case(uint64_t digits) {
    const uint64_t guard = 192ULL;
    char ref_path[256], d2_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/d2bs_ref_%llu.txt",
             (unsigned long long)digits);
    snprintf(d2_path, sizeof(d2_path), "build/d2bs_%llu.txt",
             (unsigned long long)digits);

    double tref = 0.0, td2 = 0.0;
    unsigned long nref = 0UL, nd2 = 0UL;
    if (!chud_bs_pi(digits, guard, ref_path, &tref, &nref) ||
        !d2_bs_pi(digits, guard, d2_path, &td2, &nd2)) return 0;

    int eq = files_equal(ref_path, d2_path);
    printf("%llu,%.9f,%lu,%.9f,%lu,%.6f,%d\n",
           (unsigned long long)digits,
           tref, nref,
           td2, nd2,
           td2 / tref,
           eq);
    fflush(stdout);
    return eq;
}

int main(int argc, char **argv) {
    printf("digits,chud_bs_seconds,chud_terms,d2_bs_seconds,d2_terms,d2_over_chud_bs,equal\n");
    if (argc <= 1) {
        const uint64_t defaults[] = {1000ULL, 10000ULL, 30000ULL, 100000ULL, 200000ULL};
        for (size_t i = 0; i < sizeof(defaults) / sizeof(defaults[0]); ++i)
            if (!run_d2_bs_case(defaults[i])) return 1;
        return 0;
    }
    for (int i = 1; i < argc; ++i) {
        uint64_t digits = strtoull(argv[i], NULL, 10);
        if (digits == 0ULL || !run_d2_bs_case(digits)) return 1;
    }
    return 0;
}
