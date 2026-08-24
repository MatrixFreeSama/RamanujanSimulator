#define main original_d2_benchmark_main
#include "benchmark_original_d2_basic_vs_chud.c"
#undef main

/*
 * Accelerated D2 only.  This is the class-number-2 formula with the same
 * product-tree idea used by Chudnovsky binary splitting.
 *
 * Let x = H2/J.  Then x^2 + H1*x + H2 = 0 and z = 1728*x/H2.
 * The exact linear weight also lies in Q(x):
 *
 *   beta = 1-alpha
 *        = (166325935146432054371646921561600 - 66389*x)
 *          / 1799585700948647322182745887577600.
 *
 * Therefore the tree can carry one weighted sum W directly, rather than
 * separate F and theta(F) accumulators.  Each merge carries only Q, P and W,
 * matching the P/Q/T structure of ordinary binary splitting as closely as a
 * quadratic coefficient ring permits.
 */

static const char *D2_BETA0 = "166325935146432054371646921561600";
static const char *D2_BETA1 = "-66389";
static const char *D2_BETAD = "1799585700948647322182745887577600";

typedef struct {
    mpz_t a;
    mpz_t b;
} d2_ring_t;

typedef struct {
    mpz_t q;
    d2_ring_t p;
    d2_ring_t w;
} d2_seg_t;

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

/* Three-large-product multiplication in Z[x]/(x^2+H1*x+H2). */
static void d2_ring_mul(d2_ring_t *out,
                        const d2_ring_t *x,
                        const d2_ring_t *y,
                        const mpz_t h1,
                        const mpz_t h2,
                        mpz_t m0,
                        mpz_t m1,
                        mpz_t m2,
                        mpz_t t0,
                        mpz_t t1) {
    mpz_mul(m0, x->a, y->a);
    mpz_mul(m1, x->b, y->b);
    mpz_add(t0, x->a, x->b);
    mpz_add(t1, y->a, y->b);
    mpz_mul(m2, t0, t1);

    mpz_mul(t0, h2, m1);
    mpz_sub(out->a, m0, t0);

    mpz_sub(m2, m2, m0);
    mpz_sub(m2, m2, m1);
    mpz_mul(t0, h1, m1);
    mpz_sub(out->b, m2, t0);
}

static void d2_leaf(d2_seg_t *out,
                    unsigned long n,
                    const mpz_t h2,
                    const mpz_t beta0,
                    const mpz_t beta1,
                    const mpz_t betad) {
    mpz_t num, den, g, weight0;
    mpz_inits(num, den, g, weight0, NULL);

    /* r_n = A_n*1728*x / (72*(n+1)^3*H2). */
    mpz_set_ui(num, 6UL * n + 1UL);
    mpz_mul_ui(num, num, 2UL * n + 1UL);
    mpz_mul_ui(num, num, 6UL * n + 5UL);
    mpz_mul_ui(num, num, 1728UL);

    mpz_set_ui(den, n + 1UL);
    mpz_mul_ui(den, den, n + 1UL);
    mpz_mul_ui(den, den, n + 1UL);
    mpz_mul_ui(den, den, 72UL);
    mpz_mul(den, den, h2);

    mpz_gcd(g, num, den);
    mpz_divexact(num, num, g);
    mpz_divexact(out->q, den, g);

    mpz_set_ui(out->p.a, 0UL);
    mpz_set(out->p.b, num);

    /*
     * W_leaf = beta + 6n.
     * The segment invariant is W = w / (BETAD*q), hence multiply the exact
     * ring weight by q at the leaf.
     */
    mpz_mul_ui(weight0, betad, 6UL * n);
    mpz_add(weight0, weight0, beta0);
    mpz_mul(out->w.a, weight0, out->q);
    mpz_mul(out->w.b, beta1, out->q);

    mpz_clears(num, den, g, weight0, NULL);
}

static void d2_bs_build(d2_seg_t *out,
                        unsigned long a,
                        unsigned long b,
                        const mpz_t h1,
                        const mpz_t h2,
                        const mpz_t beta0,
                        const mpz_t beta1,
                        const mpz_t betad) {
    if (b - a == 1UL) {
        d2_leaf(out, a, h2, beta0, beta1, betad);
        return;
    }

    unsigned long m = a + (b - a) / 2UL;
    d2_seg_t left, right;
    d2_seg_init(&left);
    d2_seg_init(&right);
    d2_bs_build(&left, a, m, h1, h2, beta0, beta1, betad);
    d2_bs_build(&right, m, b, h1, h2, beta0, beta1, betad);

    mpz_t m0, m1, m2, t0, t1;
    mpz_inits(m0, m1, m2, t0, t1, NULL);
    d2_ring_t rw;
    d2_ring_init(&rw);

    mpz_mul(out->q, left.q, right.q);
    d2_ring_mul(&out->p, &left.p, &right.p,
                h1, h2, m0, m1, m2, t0, t1);

    /* W = W_L + P_L*W_R. */
    d2_ring_mul(&rw, &left.p, &right.w,
                h1, h2, m0, m1, m2, t0, t1);
    mpz_mul(out->w.a, left.w.a, right.q);
    mpz_add(out->w.a, out->w.a, rw.a);
    mpz_mul(out->w.b, left.w.b, right.q);
    mpz_add(out->w.b, out->w.b, rw.b);

    d2_ring_clear(&rw);
    mpz_clears(m0, m1, m2, t0, t1, NULL);
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

    mpz_t h1, h2, beta0, beta1, betad, disc, tmpz;
    mpz_inits(h1, h2, beta0, beta1, betad, disc, tmpz, NULL);
    if (mpz_set_str(h1, D2_H1, 10) != 0 ||
        mpz_set_str(h2, D2_H2, 10) != 0 ||
        mpz_set_str(beta0, D2_BETA0, 10) != 0 ||
        mpz_set_str(beta1, D2_BETA1, 10) != 0 ||
        mpz_set_str(betad, D2_BETAD, 10) != 0) {
        mpz_clears(h1, h2, beta0, beta1, betad, disc, tmpz, NULL);
        return 0;
    }

    d2_seg_t tree;
    d2_seg_init(&tree);
    d2_bs_build(&tree, 0UL, terms, h1, h2, beta0, beta1, betad);

    mpz_mul(disc, h1, h1);
    mpz_mul_ui(tmpz, h2, 4UL);
    mpz_sub(disc, disc, tmpz);

    mpf_t fh1, fh2, fdisc, root, x, fw0, fw1, fq, fbd,
          W, y, z, K, tmp, inv_pi, pi;
    mpf_init2(fh1, bits); mpf_init2(fh2, bits);
    mpf_init2(fdisc, bits); mpf_init2(root, bits);
    mpf_init2(x, bits); mpf_init2(fw0, bits);
    mpf_init2(fw1, bits); mpf_init2(fq, bits);
    mpf_init2(fbd, bits); mpf_init2(W, bits);
    mpf_init2(y, bits); mpf_init2(z, bits);
    mpf_init2(K, bits); mpf_init2(tmp, bits);
    mpf_init2(inv_pi, bits); mpf_init2(pi, bits);

    /* Stable small embedding x = -2H2/(H1+sqrt(H1^2-4H2)). */
    mpf_set_z(fh1, h1);
    mpf_set_z(fh2, h2);
    mpf_set_z(fdisc, disc);
    mpf_sqrt(root, fdisc);
    mpf_add(tmp, fh1, root);
    mpf_mul_ui(x, fh2, 2UL);
    mpf_neg(x, x);
    mpf_div(x, x, tmp);

    /* W = (w.a + w.b*x)/(BETAD*q). */
    mpf_set_z(fw0, tree.w.a);
    mpf_set_z(fw1, tree.w.b);
    mpf_mul(fw1, fw1, x);
    mpf_add(W, fw0, fw1);
    mpf_set_z(fq, tree.q);
    mpf_div(W, W, fq);
    mpf_set_z(fbd, betad);
    mpf_div(W, W, fbd);

    mpf_div(y, x, fh2);
    mpf_mul_ui(z, y, 1728UL);

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

    mpf_clears(fh1, fh2, fdisc, root, x, fw0, fw1, fq, fbd,
               W, y, z, K, tmp, inv_pi, pi, NULL);
    d2_seg_clear(&tree);
    mpz_clears(h1, h2, beta0, beta1, betad, disc, tmpz, NULL);
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
