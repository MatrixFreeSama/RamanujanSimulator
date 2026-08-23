#include "ramanujan_c_common.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static long double minus_log10_abs(const mpf_t z){
    if(mpf_sgn(z)==0)return INFINITY;
    mpf_t a;mpf_init2(a,mpf_get_prec(z));mpf_abs(a,z);mp_exp_t e2=0;double m=mpf_get_d_2exp(&e2,a);mpf_clear(a);
    return -(log10l(fabsl((long double)m))+(long double)e2*log10l(2.0L));
}

static int run_case(int p,int layers,unsigned long dps){
    RjPoly Phi={0},P={0}; int target_hi=(p<=3)?50:60; if(!rj_generate_phi_prime(p,target_hi,&Phi))return 0; if(!rj_z_polynomial_from_phi(&Phi,&P)){rj_poly_clear(&Phi);return 0;}
    mp_bitcnt_t bits=rj_bits_for_decimal(dps,80);RjState a,b;rj_state_init(&a,bits);rj_state_init(&b,bits);rj_initial_state(&a);RjState *cur=&a,*nxt=&b;
    printf("\n=== C nested modular branch p=%d ===\n",p);printf("z-polynomial degree: %d %d\n",P.degree_x,P.degree_y);
    printf("layer 0: -log10(|z|) = %.24Lg\n",minus_log10_abs(cur->z));
    for(int k=1;k<=layers;k++){
        if(!rj_transform_step(&P,p,14,cur,nxt)){fprintf(stderr,"transform failure p=%d layer=%d\n",p,k);return 0;}
        RjState *sw=cur;cur=nxt;nxt=sw;printf("layer %d: -log10(|z|) = %.24Lg\n",k,minus_log10_abs(cur->z));
    }
    rj_state_clear(&a);rj_state_clear(&b);rj_poly_clear(&P);rj_poly_clear(&Phi);return 1;
}

int main(void){
    if(!run_case(2,4,500)||!run_case(3,3,800)||!run_case(5,3,2600))return 1;
    long double R0=3.0L*log10l(53360.0L);long double A=R0+log10l(1728.0L);
    puts("\n=== 1000-layer asymptotic digits-per-final-series-term scale ===");
    for(int p=2;p<=5;p+=(p==2?1:2)){
        if (p == 4) continue;
        long double log10R = log10l(A) + 1000.0L * log10l((long double)p);
        long double e = floorl(log10R);
        long double m = powl(10.0L, log10R - e);
        printf("p=%d: ~ %.18Lfe%.0Lf digits\n",p,m,e);
    }
    return 0;
}
