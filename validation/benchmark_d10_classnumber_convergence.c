#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Class-number-10 convergence probe for the principal CM point of
 * discriminant Delta = -13843.
 *
 * H_D(j) has degree 10.  With y = 1/j, the reciprocal class polynomial is
 *
 *   R_D(y) = 1 + h9*y + h8*y^2 + ... + h0*y^10.
 *
 * The desired principal real root is the tiny negative root near -1/h9.
 * The Ramanujan/Chudnovsky level-1 direct-series parameter is z = 1728*y.
 * This program determines z from the integer class polynomial only.  It does
 * not use pi to construct z.  The asymptotic decimal gain of one series term
 * is -log10(|z|).
 *
 * This file intentionally certifies only the singular-modulus/convergence
 * side.  A standalone D10 pi evaluator additionally needs an independently
 * algebraic representation of the companion coefficient (alpha/beta).
 */

static const char *HD_RECIP[11] = {
    "1",
    "33685656633958669107888967803661931714647776290935715365741200718405786535119228031508290986188379856550182868818065771167002034269118770355222717680986890240000",
    "-373347885499824232804010412882330497087294334996397135698951537102944942755876372620309978804222318847935128447068247937919760067465510041830754488642742834429952000000",
    "3095576677826745001791647813730045774163789379693210329707194097504345347048291101401257306468618558671368494783528572132060746016311647242448231212195585331759677440000000000",
    "-257203233837035525461655402777787102338440698835254775521328981264810710382429441484236329698930431430416431736284956301679573229169891425668900992511329873883243216896000000000000",
    "70847587572824648936525224372543186878973589986261303386779554786605638312391735378070425335097090463401517210305516918216139834018597205809483527019979311359260098560000000000000000000",
    "-284990521388995246316162827656061115574939004181220209497083465422269713552869415073979491221639144085206093288268428426324235462148160575299769337129153682332826491420672000000000000000000",
    "418369010729471937388837995958372278400439494046430042571202929122465271027869228122325714946016657143974277552123652771875407817385041791004748472603827319008838266388480000000000000000",
    "-217350279991328190761985550176942873049694056194412714001325339136414424489817793868040748634553400591278443383092632547675780889787677168496425167768288712702245816238080000000000000000000000000",
    "-16384491475280569857609606791335695552666625484765086387174938658663813020359754659184361553122188956373339449603305081745136824483583488172583077501846833379783380828160000000000000000000000000000",
    "-370379678802534904801278783784712516669511657545292686541317920733750417605390874497797114616248919946333551664716964637548220313449661671973306049203399881120521650176000000000000000000000000000000"
};

static int load_coeffs(mpf_t c[11], mp_bitcnt_t bits) {
    mpz_t z;
    mpz_init(z);
    for (int i = 0; i <= 10; ++i) {
        mpf_init2(c[i], bits);
        if (mpz_set_str(z, HD_RECIP[i], 10) != 0) {
            mpz_clear(z);
            return 0;
        }
        mpf_set_z(c[i], z);
    }
    mpz_clear(z);
    return 1;
}

static void eval_poly_and_derivative(const mpf_t y, mpf_t c[11],
                                     mpf_t p, mpf_t dp, mpf_t tmp) {
    mpf_set(p, c[10]);
    mpf_set_ui(dp, 0UL);
    for (int k = 9; k >= 0; --k) {
        mpf_mul(dp, dp, y);
        mpf_add(dp, dp, p);
        mpf_mul(tmp, p, y);
        mpf_add(p, tmp, c[k]);
    }
}

int main(void) {
    const mp_bitcnt_t bits = 4096;
    mpf_t c[11];
    if (!load_coeffs(c, bits)) return 2;

    mpf_t y, p, dp, step, tmp, z, az;
    mpf_init2(y, bits);    mpf_init2(p, bits);
    mpf_init2(dp, bits);   mpf_init2(step, bits);
    mpf_init2(tmp, bits);  mpf_init2(z, bits);
    mpf_init2(az, bits);

    /* Principal root: y ~ -1/h9. */
    mpf_set_si(y, -1L);
    mpf_div(y, y, c[1]);

    for (int iter = 0; iter < 24; ++iter) {
        eval_poly_and_derivative(y, c, p, dp, tmp);
        mpf_div(step, p, dp);
        mpf_sub(y, y, step);
    }

    mpf_mul_ui(z, y, 1728UL);
    mpf_abs(az, z);

    mp_exp_t exp2 = 0;
    double mant = mpf_get_d_2exp(&exp2, az);
    long double decimal_depth =
        -(log10l((long double)mant) + (long double)exp2 * log10l(2.0L));

    printf("Delta=-13843\n");
    printf("class_number=10\n");
    printf("algebraic_degree_of_j=10\n");
    printf("z=\n");
    mpf_out_str(stdout, 10, 90, z);
    putchar('\n');
    printf("digits_per_term=%.18Lf\n", decimal_depth);
    printf("D2_reference_digits_per_term=24.955899657654266731\n");
    printf("Chud_reference_digits_per_term=14.181647462725477656\n");
    printf("gain_over_D2=%.9Lf\n", decimal_depth / 24.955899657654266731L);
    printf("gain_over_Chud=%.9Lf\n", decimal_depth / 14.181647462725477656L);

    for (unsigned long digits = 1000UL; digits <= 1000000UL; digits *= 10UL) {
        unsigned long terms = (unsigned long)ceill(((long double)digits + 32.0L) / decimal_depth) + 1UL;
        printf("predicted_terms_%lu=%lu\n", digits, terms);
    }

    for (int i = 0; i <= 10; ++i) mpf_clear(c[i]);
    mpf_clears(y, p, dp, step, tmp, z, az, NULL);
    return 0;
}
