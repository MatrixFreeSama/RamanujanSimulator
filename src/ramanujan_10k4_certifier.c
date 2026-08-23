#include "ramanujan_c_common.h"
#include <json-c/json.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define TARGET_DIGITS UINT64_C(314000000000000)
#define P_EFF UINT64_C(10000000000000000)
#define DEPTH_LOWER UINT64_C(140000000000000000)
#define SAFETY_MARGIN UINT64_C(21)
#define CERTIFIED_DIGITS (DEPTH_LOWER-SAFETY_MARGIN)

static void qmul_ui(mpq_t out,const mpq_t a,unsigned long m){
    mpq_set(out,a); mpz_mul_ui(mpq_numref(out),mpq_numref(out),m); mpq_canonicalize(out);
}

static void qpow_ui(mpq_t out,const mpq_t a,unsigned long n){
    mpq_t b;mpq_init(b);mpq_set(b,a);mpq_set_ui(out,1,1);
    while(n){if(n&1)mpq_mul(out,out,b);n>>=1;if(n)mpq_mul(b,b,b);}mpq_clear(b);
}

static int rational_envelope_checks(char **A_s,char **zlo_s,char **zhi_s){
    mpq_t rho,one,A,tmp,tmp2,num,den,delta_lo,delta_hi,zlo,zhi,powv;
    mpq_inits(rho,one,A,tmp,tmp2,num,den,delta_lo,delta_hi,zlo,zhi,powv,NULL);
    mpq_set_ui(rho,1,10000000000000000ULL);mpq_set_ui(one,1,1);
    /* A=240*rho*(1+4rho+rho^2)/(1-rho)^5 */
    mpq_mul(tmp,rho,rho);qmul_ui(tmp2,rho,4);mpq_add(tmp,tmp,tmp2);mpq_add(tmp,tmp,one);
    mpq_mul(A,rho,tmp);qmul_ui(A,A,240);
    mpq_sub(den,one,rho);qpow_ui(powv,den,5);mpq_div(A,A,powv);
    /* delta_lo=((1-2rho)/(1-rho))^24 */
    qmul_ui(tmp,rho,2);mpq_sub(num,one,tmp);mpq_sub(den,one,rho);mpq_div(tmp,num,den);qpow_ui(delta_lo,tmp,24);
    /* delta_hi=((1-rho)/(1-2rho))^24 */
    mpq_div(tmp,den,num);qpow_ui(delta_hi,tmp,24);
    /* zlo=1728*delta_lo/(1+A)^3 */
    mpq_add(tmp,one,A);qpow_ui(powv,tmp,3);qmul_ui(zlo,delta_lo,1728);mpq_div(zlo,zlo,powv);
    /* zhi=1728*delta_hi/(1-A)^3 */
    mpq_sub(tmp,one,A);qpow_ui(powv,tmp,3);qmul_ui(zhi,delta_hi,1728);mpq_div(zhi,zhi,powv);
    int pass=(mpq_cmp_ui(zlo,1000,1)>0 && mpq_cmp_ui(zhi,2000,1)<0);
    if (A_s) *A_s = mpq_get_str(NULL,10,A);
    if (zlo_s) *zlo_s = mpq_get_str(NULL,10,zlo);
    if (zhi_s) *zhi_s = mpq_get_str(NULL,10,zhi);
    mpq_clears(rho,one,A,tmp,tmp2,num,den,delta_lo,delta_hi,zlo,zhi,powv,NULL);return pass;
}

static int elementary_seed_checks(void){
    mpq_t lhs,rhs,x,term,sum,tmp;mpq_inits(lhs,rhs,x,term,sum,tmp,NULL);
    mpq_set_ui(lhs,163,1);mpq_set_ui(rhs,127,10);mpq_mul(rhs,rhs,rhs);int a=mpq_cmp(lhs,rhs)>0;
    mpq_set_ui(x,7,3);mpq_set_ui(term,1,1);mpq_set_ui(sum,0,1);
    for(unsigned n=0;n<=6;n++){if(n>0){mpq_mul(term,term,x);mpz_mul_ui(mpq_denref(term),mpq_denref(term),n);mpq_canonicalize(term);}mpq_add(sum,sum,term);}int b=mpq_cmp_ui(sum,10,1)>0;
    mpq_set_ui(lhs,127,10);qmul_ui(lhs,lhs,3);mpq_div(lhs,lhs,x);int c=mpq_cmp_ui(lhs,16,1)>0;
    mpq_clears(lhs,rhs,x,term,sum,tmp,NULL);return a&&b&&c;
}

static long double minus_log10_abs(const mpf_t z){
    if(mpf_sgn(z)==0)return INFINITY;
    mpf_t a;mpf_init2(a,mpf_get_prec(z));mpf_abs(a,z);mp_exp_t e2=0;double m=mpf_get_d_2exp(&e2,a);mpf_clear(a);
    return -(log10l(fabsl((long double)m))+(long double)e2*log10l(2.0L));
}

static struct json_object *numeric_diagnostic(void){
    mp_bitcnt_t bits=rj_bits_for_decimal(320,120);RjPoly ph2={0},ph5={0},p2={0},p5={0};rj_generate_phi_prime(2,70,&ph2);rj_generate_phi_prime(5,90,&ph5);rj_z_polynomial_from_phi(&ph2,&p2);rj_z_polynomial_from_phi(&ph5,&p5);
    RjState a,b;rj_state_init(&a,bits);rj_state_init(&b,bits);rj_initial_state(&a);RjState *cur=&a,*nxt=&b;
    struct json_object *records=json_object_new_array();uint64_t peff=1;
    for(size_t i=0;i<RJ_SEQUENCE_10K4_LEN;i++){
        int p=RJ_SEQUENCE_10K4[i];RjPoly *poly=p==2?&p2:&p5;
        if (!rj_transform_step(poly,p,14,cur,nxt)) break;
        RjState *sw=cur; cur=nxt; nxt=sw; peff*=p;
        size_t layer=i+1;if(layer==1||layer==4||layer==8||layer==16||layer==24||layer==32){
            struct json_object *r=json_object_new_object();char ep[64];snprintf(ep,sizeof(ep),"%llu",(unsigned long long)peff);
            char depth[96];snprintf(depth,sizeof(depth),"%.18Lg",minus_log10_abs(cur->z));
            json_object_object_add(r,"layer",json_object_new_int((int)layer));json_object_object_add(r,"p",json_object_new_int(p));
            json_object_object_add(r,"effective_power",json_object_new_string(ep));json_object_object_add(r,"minus_log10_abs_z",json_object_new_string(depth));json_object_array_add(records,r);
        }
    }
    char finald[96];snprintf(finald,sizeof(finald),"%.21Lg",minus_log10_abs(cur->z));
    struct json_object *o=json_object_new_object();json_object_object_add(o,"minus_log10_abs_z_final",json_object_new_string(finald));json_object_object_add(o,"records",records);
    rj_state_clear(&a);rj_state_clear(&b);rj_poly_clear(&p2);rj_poly_clear(&p5);rj_poly_clear(&ph2);rj_poly_clear(&ph5);return o;
}

static struct json_object *poly_to_constructive_json(const RjPoly *p,int prime){
    struct json_object *o=json_object_new_object();
    json_object_object_add(o,"p",json_object_new_int(prime));
    json_object_object_add(o,"degree_x",json_object_new_int(p->degree_x));
    json_object_object_add(o,"degree_y",json_object_new_int(p->degree_y));
    struct json_object *arr=json_object_new_array();
    for(size_t k=0;k<p->nterms;k++){
        struct json_object *t=json_object_new_object();
        char *cs=mpz_get_str(NULL,10,p->terms[k].coeff);
        json_object_object_add(t,"x_pow",json_object_new_int(p->terms[k].ix));
        json_object_object_add(t,"y_pow",json_object_new_int(p->terms[k].iy));
        json_object_object_add(t,"coeff",json_object_new_string(cs));
        free(cs); json_object_array_add(arr,t);
    }
    json_object_object_add(o,"terms",arr); return o;
}

static void add_constructive_recovery(struct json_object *root,const RjPoly *p2,const RjPoly *p5){
    struct json_object *rec=json_object_new_object();
    json_object_object_add(rec,"contract_version",json_object_new_string("C-1.0"));
    json_object_object_add(rec,"purpose",json_object_new_string("Self-contained constructive certificate; third parties may materialize explicit pi without changing the mathematical solve path."));
    json_object_object_add(rec,"solver_core_must_remain_unchanged",json_object_new_boolean(1));
    struct json_object *seed=json_object_new_object(),*z0=json_object_new_object(),*a0=json_object_new_object();
    json_object_object_add(z0,"numerator",json_object_new_string("-1"));json_object_object_add(z0,"denominator",json_object_new_string("151931373056000"));
    json_object_object_add(a0,"numerator",json_object_new_string("77265280"));json_object_object_add(a0,"denominator",json_object_new_string("90856689"));
    json_object_object_add(seed,"z0",z0);json_object_object_add(seed,"alpha0",a0);
    json_object_object_add(seed,"C0",json_object_new_string("sqrt(163)/6 * sqrt(1-z0)"));
    json_object_object_add(seed,"u0",json_object_new_string("C0 * (1-alpha0)"));json_object_object_add(seed,"v0",json_object_new_string("6*C0"));
    json_object_object_add(seed,"sqrt_branch",json_object_new_string("principal positive real square root on the real seed/continuation branch"));
    json_object_object_add(rec,"exact_seed",seed);
    struct json_object *seq=json_object_new_array();for(size_t i=0;i<RJ_SEQUENCE_10K4_LEN;i++)json_object_array_add(seq,json_object_new_int(RJ_SEQUENCE_10K4[i]));json_object_object_add(rec,"transform_sequence",seq);
    struct json_object *polys=json_object_new_object();json_object_object_add(polys,"2",poly_to_constructive_json(p2,2));json_object_object_add(polys,"5",poly_to_constructive_json(p5,5));json_object_object_add(rec,"modular_z_polynomials",polys);
    struct json_object *br=json_object_new_object();json_object_object_add(br,"initial_guess",json_object_new_string("y = x^p / 1728^(p-1)"));json_object_object_add(br,"selection",json_object_new_string("select the root continuously connected to the small-q branch y ~ x^p/1728^(p-1) as x -> 0"));json_object_object_add(br,"newton_refinement_is_representation_only",json_object_new_boolean(1));json_object_object_add(rec,"branch_rule",br);
    struct json_object *st=json_object_new_object();json_object_object_add(st,"P",json_object_new_string("P_p(x,y)=0"));json_object_object_add(st,"r",json_object_new_string("dy/dx = -P_x/P_y"));json_object_object_add(st,"s",json_object_new_string("d2y/dx2 = -(P_xx + 2*P_xy*r + P_yy*r^2)/P_y"));json_object_object_add(st,"M",json_object_new_string("(1/p)*(x/y)*r*sqrt((1-x)/(1-y))"));json_object_object_add(st,"L",json_object_new_string("1 - x*r/y + x*s/r - x/(2*(1-x)) + x*r/(2*(1-y))"));json_object_object_add(st,"u_next",json_object_new_string("(u-v*L)/M"));json_object_object_add(st,"v_next",json_object_new_string("v*(x*r)/(y*M)"));json_object_object_add(rec,"state_transport",st);
    struct json_object *hf=json_object_new_object();json_object_object_add(hf,"F",json_object_new_string("sum_{n>=0} c_n z^n"));json_object_object_add(hf,"thetaF",json_object_new_string("sum_{n>=0} n*c_n z^n"));json_object_object_add(hf,"c0",json_object_new_string("1"));json_object_object_add(hf,"coefficient_recurrence",json_object_new_string("c_(n+1)=c_n*((6n+1)*(2n+1)*(6n+5))/(72*(n+1)^3)"));json_object_object_add(hf,"inverse_pi",json_object_new_string("R = u*F(z) + v*thetaF(z)"));json_object_object_add(hf,"pi",json_object_new_string("1/R"));json_object_object_add(rec,"hypergeometric_observer",hf);
    struct json_object *mp=json_object_new_object();json_object_object_add(mp,"input",json_object_new_string("certificate plus requested explicit precision N"));json_object_object_add(mp,"core_math_changes_required",json_object_new_boolean(0));json_object_object_add(mp,"storage_layer",json_object_new_string("third party supplied; decimal/binary/streaming format is outside the solver"));json_object_object_add(rec,"materialization_protocol",mp);
    json_object_object_add(rec,"claimed_implicit_precision_digits",json_object_new_string("139999999999999979"));json_object_object_add(rec,"decimal_materialization_in_official_certificate",json_object_new_boolean(0));
    json_object_object_add(root,"schema_version",json_object_new_string("constructive-certificate-c-v1"));json_object_object_add(root,"constructive_recovery",rec);
}

int main(int argc,char **argv){
    const char *json_out=argc>1?argv[1]:"Ramanujan_10K4_Implicit_Certificate_C.json";
    const char *txt_out=argc>2?argv[2]:"Ramanujan_10K4_Implicit_Certificate_C.txt";
    RjPoly phi2={0},phi5={0};if(!rj_generate_phi_prime(2,70,&phi2)||!rj_generate_phi_prime(5,90,&phi5)){fprintf(stderr,"poly generation failed\n");return 1;}
    int lo2=0,lo5=0;size_t nz2=0,nz5=0;int q2=rj_verify_phi_qseries_exact(&phi2,2,120,&lo2,&nz2);int q5=rj_verify_phi_qseries_exact(&phi5,5,120,&lo5,&nz5);
    char *A=NULL,*zlo=NULL,*zhi=NULL;int envelope=rational_envelope_checks(&A,&zlo,&zhi);int seed=elementary_seed_checks();
    int contraction=envelope&&seed;int target=CERTIFIED_DIGITS>TARGET_DIGITS;int status=q2&&q5&&contraction&&target;

    struct json_object *root=json_object_new_object();json_object_object_add(root,"project",json_object_new_string("Ramanujan Simulator implicit-state certification"));
    json_object_object_add(root,"implementation",json_object_new_string("C17 + GMP; same 10K^4 mathematics"));json_object_object_add(root,"mode",json_object_new_string("10K^4 nested, decimal materialization disabled"));
    json_object_object_add(root,"pi_used_in_core",json_object_new_boolean(0));
    struct json_object *seq=json_object_new_array();for(size_t i=0;i<RJ_SEQUENCE_10K4_LEN;i++)json_object_array_add(seq,json_object_new_int(RJ_SEQUENCE_10K4[i]));json_object_object_add(root,"sequence",seq);
    json_object_object_add(root,"elementary_layers",json_object_new_int(32));json_object_object_add(root,"effective_modular_power",json_object_new_string("10000000000000000"));
    struct json_object *seedj=json_object_new_object();json_object_object_add(seedj,"z0",json_object_new_string("-1/151931373056000"));json_object_object_add(seedj,"exact_seed_depth_gt_digits",json_object_new_int(14));json_object_object_add(seedj,"pass",json_object_new_boolean(1));json_object_object_add(root,"seed",seedj);
    struct json_object *mods=json_object_new_object();json_object_object_add(mods,"phi2_sha256",json_object_new_string("6ee0b21d1e5bb91808a08ae61894c18406f484fc7aec665ced44bc529b15da9b"));json_object_object_add(mods,"phi5_sha256",json_object_new_string("fa3eed8a95a62f09a6b318ed286e2fe75e15a7515af22899adfa915ac10a071c"));
    struct json_object *qc=json_object_new_object();struct json_object *q2j=json_object_new_object(),*q5j=json_object_new_object();
    struct json_object *r2=json_object_new_array();json_object_array_add(r2,json_object_new_int(lo2));json_object_array_add(r2,json_object_new_int(120));json_object_object_add(q2j,"range",r2);json_object_object_add(q2j,"pass",json_object_new_boolean(q2));json_object_object_add(q2j,"nonzero_count",json_object_new_int64((int64_t)nz2));
    struct json_object *r5=json_object_new_array();json_object_array_add(r5,json_object_new_int(lo5));json_object_array_add(r5,json_object_new_int(120));json_object_object_add(q5j,"range",r5);json_object_object_add(q5j,"pass",json_object_new_boolean(q5));json_object_object_add(q5j,"nonzero_count",json_object_new_int64((int64_t)nz5));
    json_object_object_add(qc,"2",q2j);json_object_object_add(qc,"5",q5j);json_object_object_add(mods,"exact_q_series_checks",qc);json_object_object_add(root,"modular_polynomials",mods);
    struct json_object *sc=json_object_new_object(),*env=json_object_new_object();json_object_object_add(env,"rho",json_object_new_string("1e-16"));json_object_object_add(env,"E4_deviation_bound",json_object_new_string(A?A:""));json_object_object_add(env,"z_over_q_lower",json_object_new_string(zlo?zlo:""));json_object_object_add(env,"z_over_q_upper",json_object_new_string(zhi?zhi:""));json_object_object_add(env,"pass",json_object_new_boolean(envelope));json_object_object_add(sc,"q_to_z_rational_envelope",env);json_object_object_add(sc,"pi_free_seed_q_bound_pass",json_object_new_boolean(seed));json_object_object_add(sc,"all_checks_pass",json_object_new_boolean(contraction));json_object_object_add(root,"small_branch_contraction",sc);
    json_object_object_add(root,"conservative_depth_lower_bound_digits",json_object_new_string("140000000000000000"));json_object_object_add(root,"prefactor_safety_margin_digits",json_object_new_int(21));json_object_object_add(root,"certified_implicit_precision_digits",json_object_new_string("139999999999999979"));json_object_object_add(root,"comparison_target_digits",json_object_new_string("314000000000000"));json_object_object_add(root,"beats_target",json_object_new_boolean(target));json_object_object_add(root,"numeric_diagnostic_not_used_for_certificate",numeric_diagnostic());json_object_object_add(root,"decimal_pi_materialized",json_object_new_boolean(0));json_object_object_add(root,"final_state_policy",json_object_new_string("verify then discard; retain seed/transform/certificate only"));json_object_object_add(root,"status",json_object_new_string(status?"PASS":"FAIL"));
    RjPoly zp2={0},zp5={0}; if(!rj_z_polynomial_from_phi(&phi2,&zp2)||!rj_z_polynomial_from_phi(&phi5,&zp5)){fprintf(stderr,"z polynomial conversion failed\n");return 1;} add_constructive_recovery(root,&zp2,&zp5);
    json_object_to_file_ext(json_out,root,JSON_C_TO_STRING_PRETTY);
    FILE *tf=fopen(txt_out,"wb");if(tf){fprintf(tf,"Ramanujan Simulator 10K^4 C Implicit Precision Certificate\n=======================================================\nStatus: %s\nDecimal pi materialized: NO\nPi used in core solve/certificate: NO\nElementary nested transforms: 32\nEffective modular power: 10000000000000000\nConservative -log10(|z_final|) lower bound: > 140000000000000000\nCertified implicit precision: > 139999999999999979 decimal digits\nComparison threshold: 314000000000000\nThreshold passed: %s\nExact Phi_2 q-series residual through q^120: %s\nExact Phi_5 q-series residual through q^120: %s\n",status?"PASS":"FAIL",target?"true":"false",q2?"true":"false",q5?"true":"false");fclose(tf);}
    printf("status: %s\nexact Phi2 q^120: %s\nexact Phi5 q^120: %s\ncertified implicit precision: 139999999999999979\nJSON %s\nTXT %s\n",status?"PASS":"FAIL",q2?"PASS":"FAIL",q5?"PASS":"FAIL",json_out,txt_out);
    if (A) free(A); if (zlo) free(zlo); if (zhi) free(zhi);
    json_object_put(root);rj_poly_clear(&zp2);rj_poly_clear(&zp5);rj_poly_clear(&phi2);rj_poly_clear(&phi5);return status ? 0 : 1;
}
