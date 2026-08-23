/* ---------- Exact q-series support for finite formal residual checks ---------- */
typedef struct { int lo,hi; mpq_t *v; } QSeries;
static int qs_init(QSeries *s,int lo,int hi) {
    long long span=(long long)hi-(long long)lo+1LL;
    if(span<=0 || span>1000000LL) return 0;
    s->lo=lo;s->hi=hi; size_t n=(size_t)span; s->v=(mpq_t*)malloc(n*sizeof(mpq_t)); if(!s->v)return 0;
    for(size_t i=0;i<n;i++){mpq_init(s->v[i]);mpq_set_ui(s->v[i],0,1);} return 1;
}
static void qs_clear(QSeries *s){ if(!s->v)return; size_t n=(size_t)(s->hi-s->lo+1);for(size_t i=0;i<n;i++)mpq_clear(s->v[i]);free(s->v);s->v=NULL; }
static mpq_ptr qs_at(QSeries *s,int e){ if(e<s->lo||e>s->hi)return NULL;return s->v[e-s->lo]; }
static mpq_srcptr qs_cat(const QSeries *s,int e){ if(e<s->lo||e>s->hi)return NULL;return s->v[e-s->lo]; }
static void qs_zero(QSeries *s){for(int e=s->lo;e<=s->hi;e++)mpq_set_ui(qs_at(s,e),0,1);}
static void qs_conv(QSeries *out,const QSeries *a,const QSeries *b){
    qs_zero(out); mpq_t t; mpq_init(t);
    for(int i=a->lo;i<=a->hi;i++){mpq_srcptr ai=qs_cat(a,i);if(!ai||mpq_sgn(ai)==0)continue;
        for(int j=b->lo;j<=b->hi;j++){int e=i+j;if(e<out->lo||e>out->hi)continue;mpq_srcptr bj=qs_cat(b,j);if(!bj||mpq_sgn(bj)==0)continue;
            mpq_mul(t,ai,bj);mpq_add(qs_at(out,e),qs_at(out,e),t);
        }} mpq_clear(t);
}
static unsigned long long sigma_u64(unsigned n,unsigned k){unsigned long long s=0;for(unsigned d=1;d*d<=n;d++)if(n%d==0){unsigned long long a=1,b=1;for(unsigned q=0;q<k;q++){a*=d;b*=n/d;}s+=a;if(d*d!=n)s+=b;}return s;}
static int qs_pow(QSeries *out,const QSeries *a,unsigned n,int lo,int hi){
    QSeries r={0},tmp={0};if(!qs_init(&r,lo,hi)||!qs_init(&tmp,lo,hi))return 0;mpq_set_ui(qs_at(&r,0),1,1);
    for(unsigned k=0;k<n;k++){qs_conv(&tmp,&r,a);for(int e=lo;e<=hi;e++)mpq_set(qs_at(&r,e),qs_at(&tmp,e));}
    for(int e=out->lo;e<=out->hi;e++){mpq_srcptr q=qs_cat(&r,e);if(q)mpq_set(qs_at(out,e),q);}qs_clear(&r);qs_clear(&tmp);return 1;
}
static int j_series_exact(int N,QSeries *out){
    QSeries E4={0},E6={0},E43={0},E62={0},Delta={0},dunit={0},inv={0},ratio={0};
    int H=N+3;
    if(!qs_init(&E4,0,H)||!qs_init(&E6,0,H)||!qs_init(&E43,0,H)||!qs_init(&E62,0,H)||!qs_init(&Delta,0,H)||!qs_init(&dunit,0,N+2)||!qs_init(&inv,0,N+2)||!qs_init(&ratio,0,N+2))return 0;
    mpq_set_ui(qs_at(&E4,0),1,1);mpq_set_ui(qs_at(&E6,0),1,1);
    for(int n=1;n<=H;n++){mpq_set_si(qs_at(&E4,n),(long)(240ULL*sigma_u64((unsigned)n,3)),1);mpq_set_si(qs_at(&E6,n),-(long)(504ULL*sigma_u64((unsigned)n,5)),1);}
    qs_pow(&E43,&E4,3,0,H);qs_pow(&E62,&E6,2,0,H);
    for(int n=0;n<=H;n++){mpq_sub(qs_at(&Delta,n),qs_at(&E43,n),qs_at(&E62,n));mpz_mul_ui(mpq_denref(qs_at(&Delta,n)), mpq_denref(qs_at(&Delta,n)), 1728); mpq_canonicalize(qs_at(&Delta,n));}
    for(int n=1;n<=N+3;n++){int e=n-1;if(e<=N+2)mpq_set(qs_at(&dunit,e),qs_at(&Delta,n));}
    if(mpq_cmp_ui(qs_at(&dunit,0),1,1)!=0){qs_clear(&E4);qs_clear(&E6);qs_clear(&E43);qs_clear(&E62);qs_clear(&Delta);qs_clear(&dunit);qs_clear(&inv);qs_clear(&ratio);return 0;}
    mpq_set_ui(qs_at(&inv,0),1,1);mpq_t sum,t;mpq_init(sum);mpq_init(t);
    for(int n=1;n<=N+2;n++){mpq_set_ui(sum,0,1);for(int k=1;k<=n;k++){mpq_mul(t,qs_at(&dunit,k),qs_at(&inv,n-k));mpq_add(sum,sum,t);}mpq_neg(qs_at(&inv,n),sum);}
    mpq_clear(sum);mpq_clear(t);qs_conv(&ratio,&E43,&inv);
    qs_zero(out);for(int e=-1;e<=N;e++){mpq_srcptr q=qs_cat(&ratio,e+1);if(q)mpq_set(qs_at(out,e),q);}
    qs_clear(&E4);qs_clear(&E6);qs_clear(&E43);qs_clear(&E62);qs_clear(&Delta);qs_clear(&dunit);qs_clear(&inv);qs_clear(&ratio);return 1;
}
static int qs_qsub(QSeries *out,const QSeries *in,int p){qs_zero(out);for(int e=in->lo;e<=in->hi;e++){int ep=e*p;if(ep<out->lo||ep>out->hi)continue;mpq_srcptr q=qs_cat(in,e);if(q)mpq_set(qs_at(out,ep),q);}return 1;}

typedef struct { int i,j; } RjPair;

static int pair_is_norm(RjPair q,int deg){ return q.i==0 && q.j==deg; }

static int poly_term_cmp_desc(const void *aa,const void *bb){
    const RjTerm *a=(const RjTerm*)aa,*b=(const RjTerm*)bb;
    if(a->ix!=b->ix) return b->ix-a->ix;
    return b->iy-a->iy;
}

static int rref_solve(mpq_t *mat,size_t rows,size_t ncols,mpq_t *sol){
    size_t rank=0; mpq_t piv,factor,tmp; mpq_init(piv);mpq_init(factor);mpq_init(tmp);
    for(size_t col=0;col<ncols;col++){
        size_t pr=rank;while(pr<rows && mpq_sgn(mat[pr*(ncols+1)+col])==0)pr++;
        if(pr==rows){mpq_clear(piv);mpq_clear(factor);mpq_clear(tmp);return 0;}
        if(pr!=rank){for(size_t c=0;c<=ncols;c++)mpq_swap(mat[pr*(ncols+1)+c],mat[rank*(ncols+1)+c]);}
        mpq_set(piv,mat[rank*(ncols+1)+col]);
        for(size_t c=col;c<=ncols;c++)mpq_div(mat[rank*(ncols+1)+c],mat[rank*(ncols+1)+c],piv);
        for(size_t r=0;r<rows;r++){
            if(r==rank || mpq_sgn(mat[r*(ncols+1)+col])==0)continue;
            mpq_set(factor,mat[r*(ncols+1)+col]);
            for(size_t c=col;c<=ncols;c++){
                mpq_mul(tmp,factor,mat[rank*(ncols+1)+c]);
                mpq_sub(mat[r*(ncols+1)+c],mat[r*(ncols+1)+c],tmp);
            }
        }
        rank++;
    }
    if(rank<ncols){mpq_clear(piv);mpq_clear(factor);mpq_clear(tmp);return 0;}
    for(size_t c=0;c<ncols;c++)mpq_set(sol[c],mat[c*(ncols+1)+ncols]);
    /* check remaining rows 0 = rhs */
    for(size_t r=ncols;r<rows;r++){
        int allzero=1;for(size_t c=0;c<ncols;c++)if(mpq_sgn(mat[r*(ncols+1)+c])!=0){allzero=0;break;}
        if(allzero && mpq_sgn(mat[r*(ncols+1)+ncols])!=0){mpq_clear(piv);mpq_clear(factor);mpq_clear(tmp);return 0;}
    }
    mpq_clear(piv);mpq_clear(factor);mpq_clear(tmp);return 1;
}

int rj_generate_phi_prime(int p,int target_hi,RjPoly *out){
    if (p < 2) return 0;
    int deg = p + 1;
    int lo = -deg * (p + 1) - 2;
    int internal_hi = target_hi + deg * p + deg + 24;
    QSeries Jfull={0},Ksrc={0},J={0},K={0};
    if(!qs_init(&Jfull,-1,internal_hi+8)||!qs_init(&Ksrc,-1,internal_hi/p+12)||!qs_init(&J,-1,internal_hi)||!qs_init(&K,-p,internal_hi))return 0;
    if(!j_series_exact(internal_hi+8,&Jfull)||!j_series_exact(internal_hi/p+12,&Ksrc)){qs_clear(&Jfull);qs_clear(&Ksrc);qs_clear(&J);qs_clear(&K);return 0;}
    for (int e = -1; e <= internal_hi; e++)
        mpq_set(qs_at(&J,e), qs_at(&Jfull,e));
    qs_qsub(&K,&Ksrc,p);
    QSeries *Jp=(QSeries*)calloc((size_t)deg+1,sizeof(QSeries));QSeries *Kp=(QSeries*)calloc((size_t)deg+1,sizeof(QSeries));
    if(!Jp||!Kp)return 0;
    for(int i=0;i<=deg;i++){qs_init(&Jp[i],lo,internal_hi);qs_init(&Kp[i],lo,internal_hi);}mpq_set_ui(qs_at(&Jp[0],0),1,1);mpq_set_ui(qs_at(&Kp[0],0),1,1);
    for(int i=1;i<=deg;i++){qs_conv(&Jp[i],&Jp[i-1],&J);qs_conv(&Kp[i],&Kp[i-1],&K);}
    size_t npairs=(size_t)(deg+1)*(deg+2)/2;RjPair *pairs=(RjPair*)malloc(npairs*sizeof(RjPair));QSeries *basis=(QSeries*)calloc(npairs,sizeof(QSeries));
    size_t qi=0;for(int i=0;i<=deg;i++)for(int j=i;j<=deg;j++){pairs[qi]=(RjPair){i,j};qs_init(&basis[qi],lo,target_hi);qs_conv(&basis[qi],&Jp[i],&Kp[j]);if(i!=j){QSeries t={0};qs_init(&t,lo,target_hi);qs_conv(&t,&Jp[j],&Kp[i]);for(int e=lo;e<=target_hi;e++)mpq_add(qs_at(&basis[qi],e),qs_at(&basis[qi],e),qs_at(&t,e));qs_clear(&t);}qi++;}
    size_t norm_idx=0;for(size_t q=0;q<npairs;q++)if(pair_is_norm(pairs[q],deg)){norm_idx=q;break;}
    size_t nunknown=npairs-1;size_t *umap=(size_t*)malloc(nunknown*sizeof(size_t));size_t u=0;for(size_t q=0;q<npairs;q++)if(q!=norm_idx)umap[u++]=q;
    size_t maxrows=(size_t)(target_hi-lo+1);mpq_t *mat=(mpq_t*)malloc(maxrows*(nunknown+1)*sizeof(mpq_t));size_t rows=0;
    for(size_t k=0;k<maxrows*(nunknown+1);k++)mpq_init(mat[k]);
    for(int e=lo;e<=target_hi;e++){
        int nonzero=mpq_sgn(qs_at(&basis[norm_idx],e))!=0;for(size_t c=0;c<nunknown&&!nonzero;c++)if(mpq_sgn(qs_at(&basis[umap[c]],e))!=0)nonzero=1;if(!nonzero)continue;
        for(size_t c=0;c<nunknown;c++)mpq_set(mat[rows*(nunknown+1)+c],qs_at(&basis[umap[c]],e));
        mpq_neg(mat[rows*(nunknown+1)+nunknown],qs_at(&basis[norm_idx],e));rows++;
    }
    mpq_t *sol=(mpq_t*)malloc(nunknown*sizeof(mpq_t));for(size_t c=0;c<nunknown;c++)mpq_init(sol[c]);int ok=rref_solve(mat,rows,nunknown,sol);
    if(ok){
        size_t cap=npairs*2;out->terms=(RjTerm*)calloc(cap,sizeof(RjTerm));out->nterms=0;out->degree_x=deg;out->degree_y=deg;
        for(size_t q=0;q<npairs;q++){
            mpq_t cq;mpq_init(cq);if(q==norm_idx)mpq_set_ui(cq,1,1);else{size_t cidx=0;while(cidx<nunknown&&umap[cidx]!=q)cidx++;mpq_set(cq,sol[cidx]);}
            if(mpq_sgn(cq)!=0){if(mpz_cmp_ui(mpq_denref(cq),1)!=0){ok=0;mpq_clear(cq);break;}RjPair pr=pairs[q];
                out->terms[out->nterms].ix=pr.i;out->terms[out->nterms].iy=pr.j;mpz_init(out->terms[out->nterms].coeff);mpz_set(out->terms[out->nterms].coeff,mpq_numref(cq));out->nterms++;
                if(pr.i!=pr.j){out->terms[out->nterms].ix=pr.j;out->terms[out->nterms].iy=pr.i;mpz_init(out->terms[out->nterms].coeff);mpz_set(out->terms[out->nterms].coeff,mpq_numref(cq));out->nterms++;}
            }mpq_clear(cq);
        }
        if(ok)qsort(out->terms,out->nterms,sizeof(RjTerm),poly_term_cmp_desc);else rj_poly_clear(out);
    }
    for (size_t c = 0; c < nunknown; c++) mpq_clear(sol[c]);
    free(sol);
    for (size_t k = 0; k < maxrows * (nunknown + 1); k++) mpq_clear(mat[k]);
    free(mat);
    free(umap);
    for (size_t q = 0; q < npairs; q++) qs_clear(&basis[q]);
    free(basis);
    free(pairs);
    for (int i = 0; i <= deg; i++) {
        qs_clear(&Jp[i]);
        qs_clear(&Kp[i]);
    }
    free(Jp);
    free(Kp);
    qs_clear(&Jfull);
    qs_clear(&Ksrc);
    qs_clear(&J);
    qs_clear(&K);
    return ok;
}
