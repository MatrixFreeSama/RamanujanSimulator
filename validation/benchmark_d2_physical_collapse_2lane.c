#define main hc_serial_benchmark_main
#include "benchmark_d2_physical_collapse.c"
#undef main

#include <omp.h>

/*
 * Same physical-embedding collapse, but retain the proven two-coordinate
 * parallel schedule inside each exact lower-tree chunk. Only the few largest
 * upper merges are replaced by scalar physical-branch accumulation.
 */

static mp_bitcnt_t hc2_threshold_bits = 4096UL;

static mp_bitcnt_t hc2_max_bits(const hc_ring_t *x, const hc_ring_t *y) {
    mp_bitcnt_t m = 0, v;
    v = mpz_sizeinbase(x->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(x->b, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->a, 2); if (v > m) m = v;
    v = mpz_sizeinbase(y->b, 2); if (v > m) m = v;
    return m;
}

static void hc2_ring_mul(hc_ring_t *out,
                         const hc_ring_t *x,
                         const hc_ring_t *y,
                         hc_ws_t *ws) {
    if (hc2_max_bits(x, y) < hc2_threshold_bits) {
        hc_ring_mul(out, x, y, ws);
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

static void hc2_build(hc_seg_t *out,
                      unsigned long a,
                      unsigned long b,
                      const mpz_t z0,
                      const mpz_t z1,
                      const mpz_t zd,
                      const mpz_t b0,
                      const mpz_t b1,
                      const mpz_t bd,
                      hc_ws_t *ws) {
    if (b - a == 1UL) {
        hc_leaf(out, a, z0, z1, zd, b0, b1, bd, ws);
        return;
    }

    unsigned long m = a + (b - a) / 2UL;
    hc_seg_t left, right;
    hc_seg_init(&left);
    hc_seg_init(&right);
    hc2_build(&left, a, m, z0, z1, zd, b0, b1, bd, ws);
    hc2_build(&right, m, b, z0, z1, zd, b0, b1, bd, ws);

    mpz_mul(out->q, left.q, right.q);
    hc2_ring_mul(&out->p, &left.p, &right.p, ws);
    hc2_ring_mul(&ws->rw, &left.p, &right.w, ws);

    if (mpz_sizeinbase(right.q, 2) >= hc2_threshold_bits) {
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

    hc_seg_clear(&left);
    hc_seg_clear(&right);
}

static int hc2_pi(uint64_t digits,
                  uint64_t guard,
                  unsigned long chunk_terms,
                  const char *path,
                  double *seconds,
                  unsigned long *terms_out) {
    double t0s = now_seconds();
    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    const long double digits_per_term =
        24.95589965765426673091470594098813130578L;
    unsigned long terms =
        (unsigned long)ceill(((long double)digits + (long double)guard + 30.0L) /
                            digits_per_term) + 2UL;

    mpz_t z0, z1, zd, b0, b1, bd;
    mpz_inits(z0, z1, zd, b0, b1, bd, NULL);
    if (mpz_set_str(z0, HC_Z0, 10) != 0 ||
        mpz_set_str(z1, HC_Z1, 10) != 0 ||
        mpz_set_str(zd, HC_ZD, 10) != 0 ||
        mpz_set_str(b0, HC_B0, 10) != 0 ||
        mpz_set_str(b1, HC_B1, 10) != 0 ||
        mpz_set_str(bd, HC_BD, 10) != 0) {
        mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
        return 0;
    }

    mpf_t omega, root61, fbd, p_total, w_total, p_chunk, w_chunk,
          tmp0, tmp1, fq, tmp2, z, K, inv_pi, pi;
    mpf_init2(omega, bits); mpf_init2(root61, bits);
    mpf_init2(fbd, bits); mpf_init2(p_total, bits);
    mpf_init2(w_total, bits); mpf_init2(p_chunk, bits);
    mpf_init2(w_chunk, bits); mpf_init2(tmp0, bits);
    mpf_init2(tmp1, bits); mpf_init2(fq, bits);
    mpf_init2(tmp2, bits); mpf_init2(z, bits);
    mpf_init2(K, bits); mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    mpf_set_ui(root61, 61UL);
    mpf_sqrt(root61, root61);
    mpf_add_ui(omega, root61, 1UL);
    mpf_div_ui(omega, omega, 2UL);
    mpf_set_z(fbd, bd);
    mpf_set_ui(p_total, 1UL);
    mpf_set_ui(w_total, 0UL);

    hc_ws_t ws;
    hc_ws_init(&ws);
    for (unsigned long a = 0UL; a < terms; a += chunk_terms) {
        unsigned long b = a + chunk_terms;
        if (b > terms) b = terms;
        hc_seg_t seg;
        hc_seg_init(&seg);
        hc2_build(&seg, a, b, z0, z1, zd, b0, b1, bd, &ws);
        hc_project(&seg, omega, fbd, p_chunk, w_chunk, tmp0, tmp1, fq);
        mpf_mul(tmp2, p_total, w_chunk);
        mpf_add(w_total, w_total, tmp2);
        mpf_mul(p_total, p_total, p_chunk);
        hc_seg_clear(&seg);
    }
    hc_ws_clear(&ws);

    mpf_set_z(tmp0, z0);
    mpf_set_z(tmp1, z1);
    mpf_mul(tmp1, tmp1, omega);
    mpf_add(z, tmp0, tmp1);
    mpf_set_z(fq, zd);
    mpf_div(z, z, fq);

    mpf_set_ui(K, 427UL);
    mpf_sqrt(K, K);
    mpf_div_ui(K, K, 6UL);
    mpf_ui_sub(tmp0, 1UL, z);
    mpf_sqrt(tmp0, tmp0);
    mpf_mul(K, K, tmp0);
    mpf_mul(inv_pi, K, w_total);
    mpf_ui_div(pi, 1UL, inv_pi);

    int ok = write_pi(path, pi, digits);
    if (seconds) *seconds = now_seconds() - t0s;
    if (terms_out) *terms_out = terms;

    mpf_clears(omega, root61, fbd, p_total, w_total, p_chunk, w_chunk,
               tmp0, tmp1, fq, tmp2, z, K, inv_pi, pi, NULL);
    mpz_clears(z0, z1, zd, b0, b1, bd, NULL);
    return ok;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s DIGITS CHUNK_TERMS\n", argv[0]);
        return 2;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(2);

    uint64_t digits = strtoull(argv[1], NULL, 10);
    unsigned long chunk_terms = strtoul(argv[2], NULL, 10);
    const uint64_t guard = 256ULL;
    char ref_path[256], hc_path[256];
    snprintf(ref_path, sizeof(ref_path), "build/hc2_ref_%llu_%lu.txt",
             (unsigned long long)digits, chunk_terms);
    snprintf(hc_path, sizeof(hc_path), "build/hc2_%llu_%lu.txt",
             (unsigned long long)digits, chunk_terms);

    double tc = 0.0, th = 0.0;
    unsigned long nc = 0UL, nh = 0UL;
    printf("digits,chunk_terms,chud_seconds,chud_terms,collapse2_seconds,d2_terms,collapse2_over_chud,equal\n");
    if (!chud_bs_pi(digits, guard, ref_path, &tc, &nc) ||
        !hc2_pi(digits, guard, chunk_terms, hc_path, &th, &nh)) return 1;
    int eq = files_equal(ref_path, hc_path);
    printf("%llu,%lu,%.9f,%lu,%.9f,%lu,%.6f,%d\n",
           (unsigned long long)digits, chunk_terms,
           tc, nc, th, nh, th / tc, eq);
    return eq ? 0 : 1;
}
