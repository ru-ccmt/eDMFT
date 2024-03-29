// @Copyright 2007 Kristjan Haule
// 
struct Amp{
  double m, p;
};
struct Gmpr{
  double m, p, r;
};
struct Hxy{
  double bh, pr, pm, mr;
  Hxy(double bh_=0, double pr_=0, double pm_=0, double mr_=0) : bh(bh_), pr(pr_), pm(pm_), mr(mr_){};
};

struct Gmpr1 : public Gmpr
{
};
struct Gmpr2 : public Gmpr
{
};

struct Lxy{
  double bh, mr, pr, rm, rp, rr;
  Lxy(double bh_=0, double mr_=0, double pr_=0, double rm_=0, double rp_=0, double rr_=0) :
    bh(bh_), mr(mr_), pr(pr_), rm(rm_), rp(rp_), rr(rr_) {};
};

/************** Needed in Integrate2c1 ***************/
// To make sure that those functions are gonna be inlined
// everywhere they are written explicitely for each type.
// This is important for speed.
inline void Setv(Amp& A, const functionb<Amp>& f, int i)
{
  A.m = f[i].m;
  A.p = f[i].p;
}

inline void Seti(Amp& A, const functionb<Amp>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { A.m = A.p = 0; return;}
  A.m = f[ip.i].m + ip.p*(f[ip.i+1].m - f[ip.i].m);
  A.p = f[ip.i].p + ip.p*(f[ip.i+1].p - f[ip.i].p);
}

inline void Setv(Gmpr& G, const functionb<Gmpr>& f, int i)
{
  G.m = f[i].m;
  G.p = f[i].p;
  G.r = f[i].r;
}
inline void Seti(Gmpr& G, const functionb<Gmpr>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { G.m = G.p = G.r = 0; return;}
  G.m = f[ip.i].m + ip.p*(f[ip.i+1].m - f[ip.i].m);
  G.p = f[ip.i].p + ip.p*(f[ip.i+1].p - f[ip.i].p);
  G.r = f[ip.i].r + ip.p*(f[ip.i+1].r - f[ip.i].r);
}
inline void Setv(Gmpr1& G, const functionb<Gmpr1>& f, int i)
{
  G.m = f[i].m;
  G.p = f[i].p;
  G.r = f[i].r;
}
inline void Seti(Gmpr1& G, const functionb<Gmpr1>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { G.m = G.p = G.r = 0; return;}
  G.m = f[ip.i].m + ip.p*(f[ip.i+1].m - f[ip.i].m);
  G.p = f[ip.i].p + ip.p*(f[ip.i+1].p - f[ip.i].p);
  G.r = f[ip.i].r + ip.p*(f[ip.i+1].r - f[ip.i].r);
}
inline void Setv(Gmpr2& G, const functionb<Gmpr2>& f, int i)
{
  G.m = f[i].m;
  G.p = f[i].p;
  G.r = f[i].r;
}
inline void Seti(Gmpr2& G, const functionb<Gmpr2>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { G.m = G.p = G.r = 0; return;}
  G.m = f[ip.i].m + ip.p*(f[ip.i+1].m - f[ip.i].m);
  G.p = f[ip.i].p + ip.p*(f[ip.i+1].p - f[ip.i].p);
  G.r = f[ip.i].r + ip.p*(f[ip.i+1].r - f[ip.i].r);
}
/****************************************************/
/*************** Needed for Integrate2c1 ************/
inline void Setv(Hxy& H, const functionb<Hxy>& f, int i)
{
  H.bh = f[i].bh;
  H.pr = f[i].pr;
  H.pm = f[i].pm;
  H.mr = f[i].mr;
}
inline void Seti(Hxy& H, const functionb<Hxy>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { H.bh = H.pr = H.pm = H.mr = 0; return;}
  H.bh = f[ip.i].bh + ip.p*(f[ip.i+1].bh - f[ip.i].bh);
  H.pr = f[ip.i].pr + ip.p*(f[ip.i+1].pr - f[ip.i].pr);
  H.pm = f[ip.i].pm + ip.p*(f[ip.i+1].pm - f[ip.i].pm);
  H.mr = f[ip.i].mr + ip.p*(f[ip.i+1].mr - f[ip.i].mr);
}
inline void AddProduct(const Hxy& H, const Gmpr& G, /*Ri*/Hxy& sum, double dh)
{
  sum.bh += (H.pr * G.r - H.pm * G.m) * dh;
  sum.pr += (H.bh * G.r + H.pr * G.m) * dh;
  sum.mr +=  H.mr * G.p * dh;
}
inline void AddProduct(const Gmpr& G, const Hxy& H, /*Ri*/Hxy& sum, double dh)
{
  sum.bh += (H.pr * G.r - H.pm * G.m) * dh;
  sum.pr += (H.bh * G.r + H.pr * G.m) * dh;
  sum.mr +=  H.mr * G.p * dh;
}
inline void Setv(Lxy& L, const functionb<Lxy>& f, int i)
{
  L.bh = f[i].bh;
  L.mr = f[i].mr;
  L.pr = f[i].pr;
  L.rm = f[i].rm;
  L.rp = f[i].rp;
  L.rr = f[i].rr;
}
inline void Seti(Lxy& L, const functionb<Lxy>& f, const intpar& ip)
{
  if (ip.p<0 || ip.p>1) { L.bh=0; L.mr=0; L.pr=0; L.rm=0; L.rp=0; L.rr=0; return;}
  L.bh = f[ip.i].bh + ip.p*(f[ip.i+1].bh - f[ip.i].bh);
  L.mr = f[ip.i].mr + ip.p*(f[ip.i+1].mr - f[ip.i].mr);
  L.pr = f[ip.i].pr + ip.p*(f[ip.i+1].pr - f[ip.i].pr);
  L.rm = f[ip.i].rm + ip.p*(f[ip.i+1].rm - f[ip.i].rm);
  L.rp = f[ip.i].rp + ip.p*(f[ip.i+1].rp - f[ip.i].rp);
  L.rr = f[ip.i].rr + ip.p*(f[ip.i+1].rr - f[ip.i].rr);
}
inline void AddProduct(const Lxy& L, const Amp& A, Lxy& sum, double dh)
{
  double Amdh = A.m*dh, Apdh = A.p*dh;
  sum.rr += L.bh * Amdh;
  sum.rp += L.mr * Apdh;
  sum.rm += L.pr * Amdh;
  sum.pr += L.rm * Apdh;
  sum.mr += L.rp * Amdh;
  sum.bh += L.rr * Apdh;
}
inline void AddProduct(const Amp& A, const Lxy& L, Lxy& sum, double dh)
{
  double Amdh = A.m*dh, Apdh = A.p*dh;
  sum.rr += L.bh * Amdh;
  sum.rp += L.mr * Apdh;
  sum.rm += L.pr * Amdh;
  sum.pr += L.rm * Apdh;
  sum.mr += L.rp * Amdh;
  sum.bh += L.rr * Apdh;
}
/***************************************************/
/*********** Needed in the outside loop ************/
inline double Product(const Hxy& h, const /*Ri&*/Hxy& r)
{
  return h.bh * r.bh + h.pr * r.pr + h.mr * r.mr;
}
inline double Product(const Lxy& L1, const Lxy& L2)
{
  return L1.rr * L2.rr + L1.rp * L2.rp + L1.rm * L2.rm + L1.pr * L2.pr + L1.mr * L2.mr + L1.bh * L2.bh;
}
/***************************************************/
/******* Needed for MakeCombinedFunctions **********/
inline void Setv(Hxy& H, const Amp& A, const Gmpr& G)
{
  H.bh = A.p*G.m+A.m*G.p;
  H.pr = A.p*G.r;
  H.pm = A.p*G.m;
  H.mr = A.m*G.r;
}
inline void Setv(Hxy& H, const Gmpr& G, const Amp& A)
{
  H.bh = A.p*G.m+A.m*G.p;
  H.pr = A.p*G.r;
  H.pm = A.p*G.m;
  H.mr = A.m*G.r;
}
inline void Setv(Lxy& L, const Gmpr1& Gf, const Gmpr2& Gb)
{
  L.bh = Gf.m * Gb.p + Gf.p * Gb.m;
  L.mr = Gf.m * Gb.r;
  L.pr = Gf.p * Gb.r;
  L.rm = Gf.r * Gb.m;
  L.rp = Gf.r * Gb.p;
  L.rr = Gf.r * Gb.r;
}
inline void Setv(Lxy& L, const Gmpr2& Gb, const Gmpr1& Gf)
{
  L.bh = Gf.m * Gb.p + Gf.p * Gb.m;
  L.mr = Gf.m * Gb.r;
  L.pr = Gf.p * Gb.r;
  L.rm = Gf.r * Gb.m;
  L.rp = Gf.r * Gb.p;
  L.rr = Gf.r * Gb.r;
}

/***************************************************/



///******************* Just for debugging *********************/
//void PrintProduct(ostream& stream, const Hxy& H, const Gmpr& G)
//{
//  stream <<(H.pr * G.r - H.pm * G.m)<<" "<<(H.bh * G.r + H.pr * G.m)<<" "<<H.mr * G.p<<endl;
//}
//void PrintProduct(ostream& stream, const Gmpr& G, const Hxy& H)
//{
//  stream <<(H.pr * G.r - H.pm * G.m)<<" "<<(H.bh * G.r + H.pr * G.m)<<" "<<H.mr * G.p<<endl;
//}
//
