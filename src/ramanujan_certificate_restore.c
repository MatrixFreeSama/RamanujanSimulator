#include "ramanujan_c_common.h"
#include <json-c/json.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static void usage(const char *argv0) {
    fprintf(stderr,"Usage: %s CERTIFICATE [--digits N] [--guard N] [--output FILE]\n",argv0);
}

static const char *jstr(struct json_object *o, const char *key) {
    struct json_object *v=NULL;
    if(!json_object_object_get_ex(o,key,&v)) return NULL;
    return json_object_get_string(v);
}

static int poly_from_json(struct json_object *node, RjPoly *out) {
    struct json_object *terms=NULL,*dx=NULL,*dy=NULL;
    if(!json_object_object_get_ex(node,"terms",&terms) || !json_object_is_type(terms,json_type_array)) return 0;
    if(!json_object_object_get_ex(node,"degree_x",&dx) || !json_object_object_get_ex(node,"degree_y",&dy)) return 0;
    size_t n=json_object_array_length(terms);
    memset(out,0,sizeof(*out)); out->degree_x=json_object_get_int(dx); out->degree_y=json_object_get_int(dy); out->nterms=n;
    out->terms=(RjTerm*)calloc(n,sizeof(RjTerm)); if(!out->terms) return 0;
    for(size_t k=0;k<n;k++) {
        struct json_object *t=json_object_array_get_idx(terms,k),*ix=NULL,*iy=NULL,*c=NULL;
        if(!json_object_object_get_ex(t,"x_pow",&ix)||!json_object_object_get_ex(t,"y_pow",&iy)||!json_object_object_get_ex(t,"coeff",&c)){rj_poly_clear(out);return 0;}
        out->terms[k].ix=json_object_get_int(ix); out->terms[k].iy=json_object_get_int(iy); mpz_init(out->terms[k].coeff);
        if(mpz_set_str(out->terms[k].coeff,json_object_get_string(c),10)!=0){rj_poly_clear(out);return 0;}
    }
    return 1;
}

static int seed_from_json(struct json_object *rec, RjState *s) {
    struct json_object *seed=NULL,*z0=NULL,*alpha0=NULL;
    if(!json_object_object_get_ex(rec,"exact_seed",&seed)) return 0;
    if(!json_object_object_get_ex(seed,"z0",&z0)||!json_object_object_get_ex(seed,"alpha0",&alpha0)) return 0;
    const char *zn=jstr(z0,"numerator"),*zd=jstr(z0,"denominator"),*an=jstr(alpha0,"numerator"),*ad=jstr(alpha0,"denominator");
    if(!zn||!zd||!an||!ad) return 0;
    mp_bitcnt_t bits=mpf_get_prec(s->z);
    mpz_t n,d; mpz_init(n);mpz_init(d);
    mpf_t a,C,tmp; mpf_init2(a,bits);mpf_init2(C,bits);mpf_init2(tmp,bits);
    mpz_set_str(n,zn,10);mpz_set_str(d,zd,10);mpf_set_z(s->z,n);mpf_set_z(tmp,d);mpf_div(s->z,s->z,tmp);
    mpz_set_str(n,an,10);mpz_set_str(d,ad,10);mpf_set_z(a,n);mpf_set_z(tmp,d);mpf_div(a,a,tmp);
    mpf_ui_sub(tmp,1,s->z);mpf_sqrt(tmp,tmp);mpf_set_ui(C,163);mpf_sqrt(C,C);mpf_div_ui(C,C,6);mpf_mul(C,C,tmp);
    mpf_ui_sub(tmp,1,a);mpf_mul(s->u,C,tmp);mpf_mul_ui(s->v,C,6);
    mpz_clear(n);mpz_clear(d);mpf_clear(a);mpf_clear(C);mpf_clear(tmp); return 1;
}

int main(int argc,char **argv) {
    if(argc<2){usage(argv[0]);return 2;}
    const char *cert_path=argv[1],*out_path=NULL; uint64_t digits=1000,guard=100;
    for(int i=2;i<argc;i++) {
        if(strcmp(argv[i],"--digits")==0 && i+1<argc) digits=strtoull(argv[++i],NULL,10);
        else if(strcmp(argv[i],"--guard")==0 && i+1<argc) guard=strtoull(argv[++i],NULL,10);
        else if(strcmp(argv[i],"--output")==0 && i+1<argc) out_path=argv[++i];
        else {usage(argv[0]);return 2;}
    }
    if(digits > (uint64_t)SIZE_MAX-2){fprintf(stderr,"digits too large for this build\n");return 2;}
    struct json_object *root=json_object_from_file(cert_path); if(!root){fprintf(stderr,"cannot parse certificate\n");return 1;}
    struct json_object *rec=NULL,*polys=NULL,*p2n=NULL,*p5n=NULL,*seq=NULL;
    if(!json_object_object_get_ex(root,"constructive_recovery",&rec)||!json_object_object_get_ex(rec,"modular_z_polynomials",&polys)||
       !json_object_object_get_ex(polys,"2",&p2n)||!json_object_object_get_ex(polys,"5",&p5n)||
       !json_object_object_get_ex(rec,"transform_sequence",&seq)) {fprintf(stderr,"certificate missing constructive fields\n");json_object_put(root);return 1;}
    RjPoly p2={0},p5={0}; if(!poly_from_json(p2n,&p2)||!poly_from_json(p5n,&p5)){fprintf(stderr,"bad polynomial table\n");json_object_put(root);return 1;}
    mp_bitcnt_t bits=rj_bits_for_decimal(digits,guard);
    RjState a,b; rj_state_init(&a,bits);rj_state_init(&b,bits); if(!seed_from_json(rec,&a)){fprintf(stderr,"bad seed\n");return 1;}
    RjState *cur=&a,*nxt=&b; size_t L=json_object_array_length(seq);
    for(size_t k=0;k<L;k++) {
        int p=json_object_get_int(json_object_array_get_idx(seq,k)); RjPoly *poly=p==2?&p2:(p==5?&p5:NULL);
        if(!poly||!rj_transform_step(poly,p,40,cur,nxt)){fprintf(stderr,"transform failed at layer %zu\n",k+1);return 1;}
        RjState *tmp=cur;cur=nxt;nxt=tmp;
    }
    mpf_t F,T,R,pi;mpf_init2(F,bits);mpf_init2(T,bits);mpf_init2(R,bits);mpf_init2(pi,bits);
    if(!rj_hyper_F_theta(cur->z,(unsigned long)(digits+guard),F,T,10000000UL)){fprintf(stderr,"hypergeometric observer failed\n");return 1;}
    mpf_mul(R,cur->u,F);mpf_mul(pi,cur->v,T);mpf_add(R,R,pi);mpf_ui_div(pi,1,R);
    FILE *fp=stdout; if(out_path){fp=fopen(out_path,"wb");if(!fp){perror("fopen");return 1;}}
    rj_write_decimal_mpf(fp,pi,(size_t)digits);fputc('\n',fp); if(out_path){fclose(fp);printf("WROTE %s (%llu digits after decimal)\n",out_path,(unsigned long long)digits);}
    mpf_clear(F);mpf_clear(T);mpf_clear(R);mpf_clear(pi);rj_state_clear(&a);rj_state_clear(&b);rj_poly_clear(&p2);rj_poly_clear(&p5);json_object_put(root);return 0;
}
