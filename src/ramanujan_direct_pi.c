#include "ramanujan_c_common.h"

#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char *argv0) {
    fprintf(stderr,
            "Usage: %s [--digits N] [--guard N] [--output FILE] [--quiet]\n",
            argv0);
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

/*
 * Direct upstream transform. The Newton loop stops when the rounded mpf
 * iterate no longer changes, so newton_cap is a safety cap rather than a
 * mandatory iteration count. The transport equations are identical to the
 * reference rj_transform_step path.
 */
static int direct_transform_step(const RjPoly *poly, int p, int newton_cap,
                                 const RjState *in, RjState *out,
                                 int *out_newton_iters) {
    mp_bitcnt_t bits = mpf_get_prec(in->z);
    mpf_t y, oldy, P, Px, Py, Pxx, Pxy, Pyy, delta;
    mpf_t r, s, M, L, tmp, tmp2, den;
    mpf_init2(y, bits);       mpf_init2(oldy, bits);
    mpf_init2(P, bits);       mpf_init2(Px, bits);
    mpf_init2(Py, bits);      mpf_init2(Pxx, bits);
    mpf_init2(Pxy, bits);     mpf_init2(Pyy, bits);
    mpf_init2(delta, bits);   mpf_init2(r, bits);
    mpf_init2(s, bits);       mpf_init2(M, bits);
    mpf_init2(L, bits);       mpf_init2(tmp, bits);
    mpf_init2(tmp2, bits);    mpf_init2(den, bits);

    /* Small branch predictor: y0 = x^p / 1728^(p-1). */
    mpf_pow_ui(y, in->z, (unsigned long)p);
    mpf_set_ui(den, 1728);
    mpf_pow_ui(den, den, (unsigned long)(p - 1));
    mpf_div(y, y, den);

    int used = 0;
    for (int it = 0; it < newton_cap; ++it) {
        rj_eval_poly_all(poly, in->z, y, P, Px, Py, Pxx, Pxy, Pyy);
        if (mpf_sgn(Py) == 0) goto fail;
        mpf_div(delta, P, Py);
        mpf_set(oldy, y);
        mpf_sub(y, y, delta);
        used = it + 1;
        if (mpf_cmp(y, oldy) == 0) break;
    }

    rj_eval_poly_all(poly, in->z, y, P, Px, Py, Pxx, Pxy, Pyy);
    if (mpf_sgn(Py) == 0 || mpf_sgn(y) == 0) goto fail;

    /* r = -Px/Py */
    mpf_div(r, Px, Py);
    mpf_neg(r, r);

    /* s = -(Pxx + 2 Pxy r + Pyy r^2)/Py */
    mpf_mul(tmp, Pxy, r);
    mpf_mul_ui(tmp, tmp, 2);
    mpf_add(tmp, tmp, Pxx);
    mpf_mul(tmp2, r, r);
    mpf_mul(tmp2, tmp2, Pyy);
    mpf_add(tmp, tmp, tmp2);
    mpf_div(s, tmp, Py);
    mpf_neg(s, s);

    /* M = (1/p)(x/y)r sqrt((1-x)/(1-y)) */
    mpf_div(M, in->z, y);
    mpf_mul(M, M, r);
    mpf_div_ui(M, M, (unsigned long)p);
    mpf_ui_sub(tmp, 1, in->z);
    mpf_ui_sub(tmp2, 1, y);
    mpf_div(tmp, tmp, tmp2);
    mpf_sqrt(tmp, tmp);
    mpf_mul(M, M, tmp);
    if (mpf_sgn(M) == 0 || mpf_sgn(r) == 0) goto fail;

    /* L = 1-xr/y+xs/r-x/(2(1-x))+xr/(2(1-y)) */
    mpf_set_ui(L, 1);
    mpf_mul(tmp, in->z, r); mpf_div(tmp, tmp, y); mpf_sub(L, L, tmp);
    mpf_mul(tmp, in->z, s); mpf_div(tmp, tmp, r); mpf_add(L, L, tmp);
    mpf_ui_sub(tmp2, 1, in->z); mpf_mul_ui(tmp2, tmp2, 2);
    mpf_div(tmp, in->z, tmp2); mpf_sub(L, L, tmp);
    mpf_ui_sub(tmp2, 1, y); mpf_mul_ui(tmp2, tmp2, 2);
    mpf_mul(tmp, in->z, r); mpf_div(tmp, tmp, tmp2); mpf_add(L, L, tmp);

    /* u'=(u-vL)/M ; v'=v*x*r/(y*M) */
    mpf_mul(tmp, in->v, L);
    mpf_sub(tmp, in->u, tmp);
    mpf_div(out->u, tmp, M);
    mpf_mul(tmp, in->v, in->z);
    mpf_mul(tmp, tmp, r);
    mpf_mul(tmp2, y, M);
    mpf_div(out->v, tmp, tmp2);
    mpf_set(out->z, y);

    if (out_newton_iters) *out_newton_iters = used;
    mpf_clear(y); mpf_clear(oldy); mpf_clear(P); mpf_clear(Px); mpf_clear(Py);
    mpf_clear(Pxx); mpf_clear(Pxy); mpf_clear(Pyy); mpf_clear(delta);
    mpf_clear(r); mpf_clear(s); mpf_clear(M); mpf_clear(L);
    mpf_clear(tmp); mpf_clear(tmp2); mpf_clear(den);
    return 1;

fail:
    mpf_clear(y); mpf_clear(oldy); mpf_clear(P); mpf_clear(Px); mpf_clear(Py);
    mpf_clear(Pxx); mpf_clear(Pxy); mpf_clear(Pyy); mpf_clear(delta);
    mpf_clear(r); mpf_clear(s); mpf_clear(M); mpf_clear(L);
    mpf_clear(tmp); mpf_clear(tmp2); mpf_clear(den);
    return 0;
}

int main(int argc, char **argv) {
    uint64_t digits = 1000;
    uint64_t guard = 140;
    const char *out_path = NULL;
    int quiet = 0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--digits") == 0 && i + 1 < argc) {
            digits = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--guard") == 0 && i + 1 < argc) {
            guard = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--quiet") == 0) {
            quiet = 1;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (digits == 0 || digits > (uint64_t)SIZE_MAX - 2) {
        fprintf(stderr, "invalid digit count\n");
        return 2;
    }

    mp_bitcnt_t bits = rj_bits_for_decimal(digits, guard);
    RjPoly p2 = {0}, p5 = {0};
    if (!rj_get_z_poly(2, &p2) || !rj_get_z_poly(5, &p5)) {
        fprintf(stderr, "cannot initialize z-polynomials\n");
        return 1;
    }

    RjState a, b;
    rj_state_init(&a, bits);
    rj_state_init(&b, bits);
    rj_initial_state(&a);
    RjState *cur = &a, *nxt = &b;

    /*
     * The observer itself can sum more than one hypergeometric term, so this
     * is a performance stop, not a mathematical validity gate. Once |z| is
     * below the requested working decimal scale, further modular layers only
     * make the observer unnecessarily over-converged.
     */
    const long double target_depth = (long double)digits + (long double)guard;
    size_t layers = 0;
    int total_newton = 0;

    while (minus_log10_abs(cur->z) < target_depth &&
           layers < RJ_SEQUENCE_10K4_LEN) {
        int p = RJ_SEQUENCE_10K4[layers];
        RjPoly *poly = (p == 2) ? &p2 : ((p == 5) ? &p5 : NULL);
        int used = 0;
        if (!poly || !direct_transform_step(poly, p, 64, cur, nxt, &used)) {
            fprintf(stderr, "direct transform failed at layer %zu\n", layers + 1);
            return 1;
        }
        total_newton += used;
        RjState *tmp_state = cur; cur = nxt; nxt = tmp_state;
        ++layers;
    }

    mpf_t F, thetaF, inv_pi, pi;
    mpf_init2(F, bits);
    mpf_init2(thetaF, bits);
    mpf_init2(inv_pi, bits);
    mpf_init2(pi, bits);

    if (!rj_hyper_F_theta(cur->z, (unsigned long)(digits + guard),
                          F, thetaF, 10000000UL)) {
        fprintf(stderr, "hypergeometric observer failed\n");
        return 1;
    }
    mpf_mul(inv_pi, cur->u, F);
    mpf_mul(pi, cur->v, thetaF);
    mpf_add(inv_pi, inv_pi, pi);
    mpf_ui_div(pi, 1, inv_pi);

    FILE *fp = stdout;
    if (out_path) {
        fp = fopen(out_path, "wb");
        if (!fp) {
            perror("fopen");
            return 1;
        }
    }
    if (!rj_write_decimal_mpf(fp, pi, (size_t)digits)) {
        fprintf(stderr, "decimal serialization failed\n");
        return 1;
    }
    fputc('\n', fp);
    if (out_path) fclose(fp);

    if (!quiet) {
        fprintf(stderr,
                "DIRECT_PI digits=%llu layers=%zu newton=%d depth=%.6Lf bits=%lu\n",
                (unsigned long long)digits, layers, total_newton,
                minus_log10_abs(cur->z), (unsigned long)bits);
    }

    mpf_clear(F); mpf_clear(thetaF); mpf_clear(inv_pi); mpf_clear(pi);
    rj_state_clear(&a); rj_state_clear(&b);
    rj_poly_clear(&p2); rj_poly_clear(&p5);
    return 0;
}
