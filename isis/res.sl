require("gsl");
() = evalfile("conv");

define rcas(fn, k) {
  variable f, r, vb, b0, n, buf, k0, s, mlev;
  variable i0, i1, nb, tr, bs;

  f = fopen(fn, "r");

  vb = 0;
  b0 = -1;
  mlev = -1;
  while (1) {
    n = fgets(&buf, f);
    if (n == -1) break;
    if (strlen(buf) < 4) continue;
    buf = buf[[:-2]];
    s = strchop(strcompress(buf, " \t"), ' ', 0);
    if (s[0] == "NELE") {
      k0 = integer(s[2]);
      if (k0 == k) vb = 1;
      else vb = 0;
    } 
    if (s[0] == "ILEV") {
      i0 = integer(s[2]);
      if (i0 < 0) continue;
      if (mlev < i0) mlev = i0;
    }
    if (b0 < 0 and vb and s[0] == "IBLK") {
      b0 = integer(s[2]);
    }
  }
  mlev++;

  r = Double_Type[mlev,mlev];
  () = fseek(f, 0, SEEK_SET);
  vb = 0;
  bs = 0;
  while (1) {
    n = fgets(&buf, f);
    if (n == -1) break;
    if (strlen(buf) < 4) continue;
    buf = buf[[:-2]];
    s = strchop(strcompress(buf, " \t"), ' ', 0);
    if (s[0] == "NELE") {
      k0 = integer(s[2]);
      if (k0 == k) vb = 1;
      else vb = 0;
      bs = 0;
    }
    if (vb == 0) continue;
    if (s[0] == "ILEV") {
      i0 = integer(s[2]);
      if (i0 < 0) {
	vb = 0;
	continue;
      }
    }
    if (s[0] == "DENS") {
      bs = 1;
      continue;
    }
    if (s[0][0] == '-') {
      vb = 0;
      continue;
    }
    if (bs and isdigit(s[0])) {
      i1 = integer(s[0]) - b0;
      if (i1 < mlev) {
	nb = atof(s[1]);
	tr = atof(s[2]);
	r[i0,i1] = tr/nb;
      }
    }
  }

  () = fclose(f);

  for (i1 = 0; i1 < mlev; i1++) {
    for (i0 = 0; i0 < i1; i0++) {
      r[i1,i1] += r[i0,i1];
    }
  }

  return r;
}
      
define rre() {
  variable f, m, i, buf, s, b, c, n, i0, g0, wet, rfn, fn;

  if (_NARGS == 2) {
    rfn = ();
    fn = ();
  } else {
    fn = ();
    rfn = NULL;
  }

  f = fopen(fn, "r");

  m = 0;
  while (1) {
    n = fgets(&buf, f);
    if (n == -1) break;

    if (strlen(buf) < 4) continue;
    
    buf = strcompress(buf[[:-2]], " \t");
    if (isdigit(buf)) m++;
  }
  
  if (m == 0) return NULL, NULL, NULL, NULL;
    
  () = fseek(f, 0, SEEK_SET);
  b = Integer_Type[5,m];
  c = Double_Type[4,m];
  i = 0;
  wet = 0;
  while (1) {
    n = fgets(&buf, f);
    if (n == -1) break;

    if (strlen(buf) < 4) continue;
    
    buf = strcompress(buf[[:-2]], " \t");
    s = strchop(buf, ' ', 0);
    if (isdigit(buf)) {
      b[0,i] = integer(s[0]);
      b[1,i] = integer(s[1]);
      b[2,i] = integer(s[4]);
      b[3,i] = integer(s[5]);
      b[4,i] = integer(s[7]);
      c[0,i] = atof(s[9]);
      c[1,i] = atof(s[11-wet]);
      c[2,i] = atof(s[13-wet]);
      c[3,i] = atof(s[12-wet]);
      i++;
    } else if (s[0] == "ILEV") {
      i0 = integer(s[2]);
    } else if (s[0] == "JLEV") {
      g0 = integer(s[2]);
    }
  }

  () = fclose(f);

  if (rfn != NULL) {
    s = readcol(rfn, 3);
    c[2,*] = s;
  }

  return i0, g0, b, c;
} 

% voigt function with unit area.
define voigt(alpha, v) {
  variable a, b, c;
  variable v2, v3, fac1, fac2;
  variable p1,p2,p3,p4,p5,p6,p7;
  variable o1,o2,o3,o4,o5,o6,o7;
  variable q1,q2;
  variable r1,r2;
  variable H, n, w, vw, i, vb;

  a=Double_Type[8];
  b=Double_Type[8];
  c=Double_Type[8];

  a[1]=122.607931777104326;
  a[2]=214.382388694706425;
  a[3]=181.928533092181549;
  a[4]=93.155580458134410;
  a[5]=30.180142196210589;
  a[6]=5.912626209773153;
  a[7]=0.564189583562615;

  b[1]=122.607931773875350;
  b[2]=352.730625110963558;
  b[3]=457.334478783897737;
  b[4]=348.703917719495792;
  b[5]=170.354001821091472;
  b[6]=53.992906912940207;
  b[7]=10.479857114260399;

  c[1]=0.5641641;
  c[2]=0.8718681;
  c[3]=1.474395;
  c[4]=-19.57862;
  c[5]=802.4513;
  c[6]=-4850.316;
  c[7]=8031.468;
 
  n = length(v);
  H = Double_Type[n];
  vb = 2.5;
  if (alpha <= .001) {
    w = where(abs(v) >= vb);
    if (length(w) > 0) {
      v2   = v[w]* v[w];
      v3   = 1.0;
      fac1 = c[1];
      fac2 = c[1] * (v2 - 1.0);     
      for (i=1;i<=7;++i) {
	v3     = v3 * v2;
	fac1 = fac1 + c[i] / v3;
	fac2 = fac2 + c[i] / v3 * (v2 - i);
      }
      H[w] = exp(-v2) * (1. + fac2*alpha^2 * (1. - 2.*v2)) + fac1 * (alpha/v2);
    }
    w = where(abs(v) < vb);
  } else {
    w = [0:n-1];
  }
  if (length(w) > 0) {
    p1 = alpha;
    vw = v[w];
    o1 = -vw;
    p2 = (p1 * alpha + o1 * vw);
    o2 = (o1 * alpha - p1 * vw);
    p3 = (p2 * alpha + o2 * vw);
    o3 = (o2 * alpha - p2 * vw);
    p4 = (p3 * alpha + o3 * vw);
    o4 = (o3 * alpha - p3 * vw);
    p5 = (p4 * alpha + o4 * vw);
    o5 = (o4 * alpha - p4 * vw);
    p6 = (p5 * alpha + o5 * vw);
    o6 = (o5 * alpha - p5 * vw);
    p7 = (p6 * alpha + o6 * vw);
    o7 = (o6 * alpha - p6 * vw);

    q1 = a[1] + p1 * a[2] + p2 * a[3] + p3 * a[4] +
      p4 * a[5] + p5 * a[6] + p6 * a[7];
    r1 =        o1 * a[2] + o2 * a[3] + o3 * a[4] +
      o4 * a[5] + o5 * a[6] + o6 * a[7];
    q2 = b[1] + p1 * b[2] + p2 * b[3] + p3 * b[4] +
      p4 * b[5] + p5 * b[6] + p6 * b[7] + p7;
    r2 =        o1 * b[2] + o2 * b[3] + o3 * b[4] +
      o4 * b[5] + o5 * b[6] + o6 * b[7] + o7;

    H[w] = (q1 * q2 + r1 * r2) / (q2 * q2 + r2 * r2);
  }

  return H/sqrt(PI);
}

define smaxwell(e, t, ep) {
  variable x, y, c, sx, sy;

  c = 0.282095;
  x = e/t;
  y = ep/t;
  sx = sqrt(x);
  sy = sqrt(y);
  c = c*(exp(-(sx-sy)^2) - exp(-(sx+sy)^2))/(sy*t);
  return c;
}
  
define maxwell(e, t) {
  variable x, c;
  
  c = 1.12837967;
  x = e/t;
  x = c*sqrt(x)*exp(-x)/t;
  
  return x;
}

define powerlaw(e, a, e0, e1) {
  variable x, w, f;

  f = Double_Type[length(e)];
  w = where(e <= e1 and e >= e0);
  if (length(w) == 0) {
    return f;
  }
  if (a == 1.0) {
    x = 1.0/(log(e1) - log(e0));
  } else {
    x = 1.0 - a;
    x = x/(e1^x - e0^x);
  }
  f = x*e^(-a);
  
  return f;
}

define fmaxwell(e, e0, t1, t2) {
  variable p2, p, a, x1, es, e0s, t1s, x2, x3, y;
  
  p2 = 1.0 - t1/t2;
  p = sqrt(p2);
  a = 1.0/(2.0*t2*p);
  x1 = (e0 - e/p2)/t2;
  es = sqrt(e);
  e0s = sqrt(e0)*p2;
  t1s = sqrt(t1)*p;
  x2 = (es + e0s)/t1s;
  x2 = log_erfc(x2)-x1;
  x3 = (es - e0s)/t1s;
  x3 = log_erfc(x3)-x1;
  y = exp(x3)*(1.0 - exp(x2-x3));
  y = y*a;

  return y;
}

define rerate(e, de, g0, b, c, mode) {
  variable r, n, p, p0, f0, f1, str, m, t, x, w, w1, i, j, flev;
  variable RBohr, Hartree, RateAU, c10, Me;

  RBohr = 5.29177E-9;
  Hartree = 27.2114;
  RateAU = 4.13414E16;
  c10 = 2.99792;
  Me = 511.99E3; 

  n = length(e);
  
  p = 1E20*PI^2*RBohr^2*Hartree^2/RateAU;
  if (mode <= 1) p *= c10*sqrt(2.0*c[0,*]/Me);
  str = p*(b[1,*]+1.0)/(2.0*(g0+1.0));
  str *= c[1,*]*c[2,*]/c[0,*];

  f0 = min(b[3,*]);
  f1 = max(b[3,*]);
  m = f1-f0+1;
  r = Double_Type[m,n];
  flev = [f0:f1];

  if (mode > 0) {
    if (length(de) == 1) {
      p0 = 1.0/sqrt(2.0*PI*de*de);
      for (i = 0; i < n; i++) {
	x = (c[0,*] - e[i])/de;
	w = where(x > -10.0 and x < 10.0);
	if (length(w) == 0) continue;
	t = str[w]*p0*exp(-0.5*x[w]*x[w]);
	for (j = f0; j <= f1; j++) {
	  w1 = where(b[3,w] == j);
	  if (length(w1) == 0) continue;      
	  r[j-f0,i] = sum(t[w1]);
	}
      }
    } else {
      p0 = sqrt((de[0]*2.0*c[0,*]+de[1]^2/8.0));
      for (i = 0; i < n; i++) {
	x = (c[0,*] - e[i])/(p0*sqrt(c[0,*]));
	w = where(x > -10.0 and x < 30.0);
	if (length(w) == 0) continue;
	t = str[w]*fmaxwell(e[i], c[0,w], de[0], de[1]);
	for (j = f0; j <= f1; j++) {
	  w1 = where(b[3,w] == j);
	  if (length(w1) == 0) continue;      
	  r[j-f0,i] = sum(t[w1]);
	}
      }
    }
  } else {
    p0 = 1.12838;
    for (i = 0; i < n; i++) {
      x = c[0,*]/e[i];
      t = p0*str*sqrt(x)*exp(-x)/e[i];
      for (j = f0; j <= f1; j++) {
	w1 = where(b[3,*] == j);
	if (length(w1) == 0) continue;
	r[j-f0,i] = sum(t[w1]);
      }
    }
  }
  
  return (flev, r);
  w = where(max(r,1) > 1E-20);
  return (flev[w],r[w,*]);
}

define mexpint1(x) {
  variable r, q2, e;

  if (x > 10) {
    return (1.0-1.0/x+2.0/(x*x)-6.0/(x*x*x))/x;
  }
  if (x < 0.02) {
    return exp(x)*(-log(x) - 0.5773 + x);
  } 
  if (x < 1.5) e = -0.5;
  else e = 0.5;
  r = (1+x)/x;
  q2 = 1.0/((1.0+x)*(1.0+x));
  return log(r) - (0.36+0.03*((x+0.01)^e))*q2;
}
  
define planckpi(y) {
  variable t, x, f, a, d, n, b;

  if (y <= 1E-6) return 1/y;
  if (y >= 10.0) return exp(-y + log(mexpint1(y)));
  
  a = exp(-10.0 + log(mexpint1(10.0)));
  
  b = log10(y);
  n = 1000;
  d = (1.0 - b)/(n-1.0);
  x = 10^[b:1.0+0.1*d:d];
  t = log(x);
  f = 1.0/(exp(x)-1.0);
  a += interp_linear_integ(t, f, t[0], t[-1]);
  return a;
}

define drsupa(n, z0, d0, t0, ic) {
  variable d, t, p, y, m, i, a, e;

  m = length(n);
  p = Double_Type[m];
  if (ic == 0) {
    a = 4.6E-8;
    d = d0/z0^7.0;
    t = t0/z0^2;
    e = const->Ryd_eV/(n*n);
    y = e/t;
    for (i = 0; i < m; i++) {
      p[i] = 1.0/(1.0 + a*n[i]^8*d*sqrt(y[i])*exp(-y[i])*mexpint1(y[i]));
    }
  } else if (ic == 1) {
    a = 0.605;
    e = const->Ryd_eV*z0*z0/(n*n);
    y = e/t0;
    for (i = 0; i < m; i++) {
      p[i] = 1.0/(1.0 + a*d0*planckpi(y[i]));
    }
  } else {
    a = 1.14E-5;
    p = 1.0/(1.0 + a*(d0/(3+t0))*exp((-2*t0-6.0)*log(z0)+(2*t0+6)*log(n)));
  }

  return p;
}

define drsupf(n, n0) {
  variable a, p, x;
  
  if (n0 < 0) {
    return 1.0/(1.0 + n*0.0);
  }
    
  x = n/n0;
  p = [3.9826, 0.7534, 4.2858, 5.5715];
  a = p[0]/(1.0+(x/p[1])^p[2]) + p[3];
  return 1.0/(1.0 + x^a);
}

define drsup(z0, d0, t0) {
  variable p0, p, y, n, n0, s, x0, x1, yn, a, d, t, xd, xt;

  if (d0 <= 0 or t0 <= 0) return -1.0;

  d = exp(log(d0)-7*log(z0));
  t  = exp(log(t0)-2.0*log(z0));
  xd = log(d);
  xt = log(t);  

  p0 = [7.966038e+00, 3.571021e-02, -1.949113e-01];
  p = [8.334987e+00, -4.561742e-02, 9.596166e-02];

  a = (d/sqrt(t))^(1.0/7.0);
  y = (p[0] + p[1]*xd + p[2]*xt)/a;
  n0 = (p0[0] + p0[1]*xd + p0[2]*xt)/a;

  s = 13.6/t;  
  a = s/(n0*n0);
  yn = n0*exp((log(mexpint1(a))-a)/7.0);
  if (yn > y) {
    x1 = n0;
    x0 = 0.5*n0;
    while (1) {
      a = s/(x0*x0);
      yn = x0*exp((log(mexpint1(a))-a)/7.0);
      if (yn < y) break;
      x0 *= 0.5;
    }
  } else {
    x0 = n0;
    x1 = 2.0*n0;
    while (1) {
      a = s/(x1*x1);
      yn = x1*exp((log(mexpint1(a))-a)/7.0);
      if (yn > y) break;
      x1 *= 2.0;
    }
  }

  while (1) {
    n0 = 0.5*(x0 + x1);
    a = s/(n0*n0);
    yn = n0*exp((log(mexpint1(a))-a)/7.0);
    if (abs(1-yn/y) < 1e-2) break;
    if (yn < y) {
      x0 = n0;
    } else {
      x1 = n0;
    }
  }

  return n0;
}

define fbr(x, p) {
  variable x3, np;

  x3 = double(x)*x*x;
  np = length(p);
  if (np == 3) {
    return (p[0]+p[1]/x3)/(p[0]+p[2]/x3);
  } else {
    return (p[0]+p[1]/x3)/(p[2]+p[3]/x3);
  }
}

define drtop(n0, n1, pbr, ns) {
  variable a, b, t, i, x, br;
  
  x = [n0:n1];
  br = Double_Type[length(x)];
  if (length(pbr) > 2) {
    br = fbr(x, pbr);
  } else {
    br[*] = 1.0;
  }
  if (ns > 0) {
    br *= drsupf(x, ns);
  }
  a = 1.0/n0;
  a = br[0]*a*a*a;
  b = 0.0;
  for (i = n0+1; i <= n1; i++) {
    t = 1.0/i;
    t = t*t*t;
    b = b + br[i-n0]*t;
  }
  
  return 1.0+b/a;
}

define branching(g0, b, c, n, channels) {
  variable p, str, wc, w, i;

  p = Double_Type[length(n)];
  wc = where(b[2,*] >= channels[0] and 
	     b[2,*] <= channels[1]);
  if (length(wc) == 0) return p;
  str = (b[1,wc]+1.0)*c[1,wc]/c[0,wc];
  for (i = 0; i < length(n); i++) {
    w = where(b[4,wc] == n[i]);
    if (length(w) == 0) continue;
    p[i] = sum(str[w]*c[2,wc[w]])/sum(str[w]);
  }

  return p;
}

define drrate(e, de, g0, b, c, mode, nmax, channels, pbr, eref) {
  variable r, n, p, p0, str, t, x, w, i, sig, er;
  variable RBohr, Hartree, RateAU, c10, Me, h;
  variable n0, c1, c2, wc, wlo, wg, alpha, j;

  n = length(e);
  r = Double_Type[n];
  if (g0 == NULL or length(b) == 0) return r;

  RBohr = 5.29177E-9;
  Hartree = 27.2114;
  RateAU = 4.13414E16;
  c10 = 2.99792;
  Me = 511.99E3; 
  h = 4.1357E-15/(2.0*PI);
  
  er = c[0,*];
  w = where(er < 0);
  if (length(w) > 0) {
    er[w] -= eref;
  }
  n0 = max(b[4,*]);
  if (length(channels) == 2) {
      wc = where(er > 0 and 
		 b[2,*] >= channels[0] and 
		 b[2,*] <= channels[1] and 
		 b[4,*] <= nmax[1] and
		 b[4,*] >= nmax[0]);      
  } else {
    wc = where(er > 0 and 
	       b[4,*] <= nmax[1] and b[4,*] >= nmax[0]);
  }
  if (length(wc) == 0) return r;

  p = 1E20*PI^2*RBohr^2*Hartree^2/RateAU;
  if (mode <= 1) p *= c10*sqrt(2.0*er[wc]/Me);
  str = p*(b[1,wc]+1.0)/(2.0*(g0+1.0));
  str *= c[1,wc]*c[2,wc]/er[wc];

  if (nmax[1] > n0) {    
    if (pbr != NULL) {
      c1 = drtop(n0, nmax[1], pbr);
      w = where(b[4,wc] == n0);
      str[w] *= c1;
    }
  }
  if (mode > 0) {
    if (length(de) == 1) {
      if (typeof(de) == Array_Type) {
	sig = sqrt(2*de[0]*er[wc]);
      } else {
	sig = er[wc];
	sig[*] = de;
      }
      p0 = 1.0/sqrt(2.0*PI*sig*sig);
      for (i = 0; i < n; i++) {
	x = (er[wc] - e[i])/sig;
	w = where(x > -10.0 and x < 10.0);
	if (length(w) == 0) continue;
	t = str[w]*p0[w]*exp(-0.5*x[w]*x[w]);
	r[i] = sum(t);
      }
    } else {
      p0 = sqrt(de[0]*2.0*er[wc]+de[1]^2/8.0);
      for (i = 0; i < n; i++) {
	x = (er[wc] - e[i])/p0;
	w = where(x > -10.0 and x < 30.0);
	if (length(w) == 0) continue;
	t = str[w]*fmaxwell(e[i], er[wc[w]], de[0], de[1]);
	r[i] = sum(t);
      }
    }
  } else if (mode == 0) {
    p0 = 1.12838;
    for (i = 0; i < n; i++) {
      x = er[wc]/e[i];
      t = p0*str*sqrt(x)*exp(-x)/e[i];
      r[i] = sum(t);
    }
  } else if (mode == -1) {
    for (i = 0; i < n; i++) {
      t = str*smaxwell(er[wc],e[i],de);
      r[i] = sum(t);
    }
  }
  return r;
}
  
define drfitf(x, p) {
  variable n, i, f;

  n = length(p);
  f = 0.0;
  for (i = 0; i < n; i += 2) {
    f += sqrt(p[i+1])*p[i+1]*p[i]*exp(-p[i+1]/x);
  }
  f /= sqrt(x)*x;
  
  return f;
}

define logdrfitf(x, p) {
  return drfitf(x, 10^p);
}

private variable arrayfun;
define marray_fit(lo, hi, p) {
  variable y;

  y = @arrayfun(lo, p);

  return y;
}

define fitarray(x, y, wt, pars, pmin, pmax, fun) {
  variable lo, hi, w, k, np, pname, i, p0, p1, s;

  set_fit_method("marquardt;tol=1e-10;max_loops=10000;delta=1E-6");
  set_fit_statistic("chisqr;sigma=data");

  w = array_sort(x);
  lo = x[w];
  hi = make_hi_grid(lo);
  
  if (wt == NULL) {
    wt = 1/y;
    set_fit_statistic("chisqr;sigma=lsq");
  }
  s = 1.0/wt;  
  k = define_counts(lo, hi, y[w], s[w]);
  arrayfun = fun;
  np = length(pars);
  pname = String_Type[np];
  for (i = 0; i < np; i++) {
    pname[i] = sprintf("P%d", i);
  }
  add_slang_function("marray", pname);
  fit_fun("marray(1)");
  
  for (i = 0; i < np; i++) {
    set_par(i+1, pars[i]);
    if (pmin == NULL) {
      p0 = 1;
    } else {
      p0 = pmin[i];
    }
    if (pmax == NULL) {
      p1 = 0;
    } else {
      p1 = pmax[i];
    }
    if (p0 == p1) {
      freeze(i+1);
    }
    if (p1 > p0) {
      set_par(i+1, pars[i], 0, p0, p1);
    }
  }
  xnotice(k);
  s = fit_counts();
  
  p0 = @pars;
  for (i = 0; i < np; i++) {
    p0[i] = get_par(i+1);
  }  

  delete_data(k);

  return (p0, s);
}

define drfit(temp, rts0, mterms, eps, p0, pmin0, pmax0) {
  variable minr, maxr, maxt, y, p, w, wts, yfit, s, rsd, tol;
  variable w1, w2, e, nterms, logt, rts, pmin, pmax;

  maxt = max(temp);
  maxr = max(rts0);
  w = where(rts0 > 1E-6*maxr);
  minr = min(rts0[w]);
  
  logt = log(temp[w]);
  rts = rts0[w]/minr;

  if (p0 == NULL) {
    e = exp(sum(rts*logt)/sum(rts));
    p0 = [1.0, e];
    y = drfitf(temp[w], p0);
    p0[0] = sum(rts)/sum(y);
    nterms = 1;
  } else {
    nterms = length(p0)/2;
    w1 = [0:2*nterms-1:2];
    w2 = w1 + 1;
    p0[w1] /= minr;
    if (pmin0 != NULL) {
      pmin0[w1] /= minr;
      pmax0[w1] /= minr;
      pmin0 = log10(pmin0);
      pmax0 = log10(pmax0);
    }
  }

  wts = 1/rts;
  p0 = log10(p0);
  while (1) {
    w1 = [0:2*nterms-1:2];
    w2 = w1 + 1;
    if (pmin0 == NULL or nterms > length(p0)/2) {
      pmin = @p0;
      pmax = @p0;
      pmin[*] = -30.0;
      pmax[*] = 30.0;
      pmin[w2] = -5.0;
      pmax[w2] = log10(maxt);
    } else {
      pmin = @pmin0;
      pmax = @pmax0;
    }
    (p,s) = fitarray(temp[w], rts, wts, p0, pmin, pmax, &logdrfitf);
    yfit = logdrfitf(temp[w], p);
    rsd = abs(rts-yfit);
    tol = max(rsd/rts);
    vmessage("NTERMS=%2d Tol=%10.3E", nterms, tol);
    if (tol < eps or nterms == mterms) break;
    s = int(urand()*length(w));
    e = log10(temp[w[s]]);
    if (e < -5.0) e = -5.0*0.99;
    p0 = [p, log10(rts[s]),e];
    nterms++;
  }
    
  w1 = [0:2*nterms-1:2];
  w2 = w1 + 1;
  s = array_sort(p[w2]);
  p[w1] = p[w1[s]];
  p[w2] = p[w2[s]];
  p = 10^p;
  p[w1] *= minr;
  yfit = drfitf(temp, p);
  return (p,yfit,tol);
}
	
define drsupfitf(x, p) {
  variable a, b;

  b = p[2]/x;
  a = p[0] + p[1]*b^p[3]*exp(-(b^p[4]));

  return a;
}

define drsupfit(t, s) {
  variable n, p, p0, w, pmin, pmax, y;

  n = length(t);
  p = Double_Type[5];
  y = -log(s);
  p[0] = y[0];  
  p[1] = max(y);
  w = where(y == max(y));
  w = w[0];
  p[2] = t[w];
  p[3] = 1.0;
  p[4] = 1.0;

  pmin = @p;
  pmax = @p;
  pmin[*] = 0.0;
  pmin[0] = 0.0;
  pmax[0] = 1E10;
  pmax[1] = 1E10;
  pmax[2] = max(t);
  pmin[3] = 0.0;
  pmax[3] = 5.0;
  pmin[4] = 0.0;
  pmax[4] = 5.0;
  p0 = @p;
  set_fit_method("marquardt;tol=1e-10;max_loops=10000;delta=1E-6");
  (p, ) = array_fit(t, y, NULL, p0, pmin, pmax, &drsupfitf);
  return p;
}
   
define population(cas, r) {
  variable p, i, j, n, m, s;

  (s,,) = array_info(r);
  n = s[0];
  m = s[1];
  p = Double_Type[n,m];
  
  for (i = n-1; i >= 1; i--) {
    p[i,*] = r[i,*];
    for (j = i+1; j < n; j++) {
      p[i,*] += p[j,*]*cas[i,j];
    }
    if (cas[i,i] > 0) {
      p[i,*] /= cas[i,i];
    } else {
      p[i,*] = 0.0;
    }
  }
  
  return p;
}

define rdump(fn, m) {
  variable f, n, sz, t;
  variable rs, i, s;
  variable dump0 = struct {
    nele, id, ib, ilev, j, ibase, vnl, energy, n, trate, br
  };
  variable dump5 = struct {
    i0, i1, rd, ri
  };
  variable Hartree=27.2113962;

  f = fopen(fn, "r");
  () = fseek(f, 0, SEEK_END);
  sz = ftell(f);
  () = fseek(f, 0, SEEK_SET);

  if (m == 0) {
    rs = 4*2+3*4+4*8;
    n = sz/rs;
    dump0.nele = Short_Type[n];
    dump0.id = Integer_Type[n];
    dump0.ib = Integer_Type[n];
    dump0.ilev = Integer_Type[n];
    dump0.j = Short_Type[n];
    dump0.ibase = Short_Type[n];
    dump0.vnl = Short_Type[n];
    dump0.energy = Double_Type[n];
    dump0.n = Double_Type[n];
    dump0.trate = Double_Type[n];
    dump0.br = Double_Type[n];
    
    for (i = 0; i < n; i++) {
      s = fread(&t, Short_Type, 1, f);
      dump0.nele[i] = t;
      s = fread(&t, Integer_Type, 1, f);
      dump0.id[i] = t;
      s = fread(&t, Integer_Type, 1, f);
      dump0.ib[i] = t;
      s = fread(&t, Integer_Type, 1, f);
      dump0.ilev[i] = t;
      s = fread(&t, Short_Type, 1, f);
      dump0.j[i] = t;
      s = fread(&t, Short_Type, 1, f);
      dump0.ibase[i] = t;
      s = fread(&t, Short_Type, 1, f);
      dump0.vnl[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump0.energy[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump0.n[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump0.trate[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump0.br[i] = t;
    }
    
    dump0.energy *= Hartree;
    () = fclose(f);
    return dump0;
  } else {
    rs = 2*4+2*8;
    n = sz/rs;
    dump5.i0 = Integer_Type[n];
    dump5.i1 = Integer_Type[n];
    dump5.rd = Double_Type[n];
    dump5.ri = Double_Type[n];
    for (i = 0; i < n; i++) {
      s = fread(&t, Integer_Type, 1, f);
      dump5.i0[i] = t;
      s = fread(&t, Integer_Type, 1, f);
      dump5.i1[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump5.rd[i] = t;
      s = fread(&t, Double_Type, 1, f);
      dump5.ri[i] = t;
    }
    
    () = fclose(f);
    return dump5;
  }
}

define mkbc(dump0, dump5, ilev0) {
  variable w, g0, bnd, b, c, nr;
  
  g0 = dump0.j[ilev0];
  w = where(dump5.i1 == ilev0);
  nr = length(w);
  if (nr == 0) return;
  bnd = dump5.i0[w];
  
  b = Integer_Type[5,nr];
  c = Double_Type[3,nr];
  b[0,*] = bnd;
  b[1,*] = dump0.j[bnd];
  b[2,*] = dump0.ibase[w];
  b[4,*] = dump0.vnl[bnd]/100;
  c[0,*] = dump0.energy[bnd] - dump0.energy[ilev0];
  c[1,*] = dump5.rd[w];
  c[2,*] = dump0.br[bnd];

  return (g0, b, c);
}

define cre(e, de, mode, ilev0, ilev1, dump0, dump5) {
  variable i, k, n1, w0, w, n0, g0, nr, i1, nlev1, rt;
  variable b, c, nr0, nr1, w1, w2, t, q, q1, s, s1;
  
  nlev1 = length(ilev1);
  rt = Double_Type[nlev1,length(e)];

  g0 = dump0.j[ilev0];
  w = where(dump5.i1 == ilev0);
  nr = length(w);
  if (nr == 0) return rt;
  
  s = array_sort(dump5.i0[w]);
  
  for (i1 = 0; i1 < nlev1; i1++) {
    w1 = where(dump5.i1 == ilev1[i1]);
    nr1 = length(w1);
    if (nr1 == 0) continue;
    s1 = array_sort(dump5.i0[w1]);    
    b = Integer_Type[5,nr1];
    c = Double_Type[3,nr1];
    nr0 = 0;
    for (i = 0, k = 0; i < nr and k < nr1;) {
      q = w[s[i]];
      q1 = w1[s1[k]];
      if (dump5.i0[q] == dump5.i0[q1]) {
	t = dump5.i0[q];
	b[0,nr0] = t;
	b[1,nr0] = dump0.j[t];
	b[2,nr0] = dump0.ibase[t];
	b[3,nr0] = ilev1[i1];
	b[4,nr0] = dump0.vnl[t]/100;
	c[0,nr0] = dump0.energy[t] - dump0.energy[ilev0];
	c[1,nr0] = dump5.rd[q];
	c[2,nr0] = dump5.rd[q1]/dump0.trate[t];
	nr0++;
	i++;
	k++;
      } else if (dump5.i0[q] < dump5.i0[q1]) {
	i++;
      } else {
	k++;
      }
    }
    if (nr0 > 0) {
      b = b[*,[0:nr0-1]];
      c = c[*,[0:nr0-1]];
      (,rt[i1,*]) = rerate(e, de, g0, b, c, mode);
    }
  }

  return rt;
}
  
define mkcas(d1, imin, imax) {
  variable w, cas, n, i;

  n = imax-imin+1;
  w = where(d1.i0 <= imax and d1.i0 >= imin);
  cas = Double_Type[n,n];
  for (i = 0; i < length(w); i++) {
    cas[d1.i1[w[i]]-imin,d1.i0[w[i]]-imin] = d1.rd[w[i]];
  }
  
  for (i = 0; i < n; i++) {
    cas[i,i] = sum(cas[*,i]);
  }

  return cas;
}

define casmat(c) {
  variable i, j, n, k, r;

  n = array_shape(c)[0];
  r = Double_Type[n,n];
  for (i = 1; i < n; i++) {
    r[i,i] = 1.0;
    for (j = i-1; j >= 0; j--) {
      r[j,i] = c[j,i]/c[i,i];
      for (k = j+1; k < i; k++) {
	r[j,i] += r[k,i]*c[j,k]/c[k,k];
      }
    }
  }
  return r;
}

define rrmx(fn, i0, i1) {
  variable f, buf, e, s, np, n0, nk, n, i;
  variable k0, k1, j0, j1, et, et0, ist, x0, x1;

  f = fopen(fn, "r");
  
  () = fgets(&buf, f);
  () = sscanf(buf, "### %d %d %d", &np, &n0, &nk);
  e = Double_Type[np];
  s = @e;

  n = 0;
  ist = 0;
  while (1) {
    i = fgets(&buf, f);
    if (i <= 0) break;
    if (strlen(buf) < 2) continue;
    if (buf[[0:1]] == "# ") {
      () = sscanf(buf, "# %d %d %d %d %f", &k0, &j0, &k1, &j1, &et);
      if (k0 != i0 or k1 != i1) {
	ist = 0;
      } else {
	ist = 1;
	et0 = et;
      }
    } else {
      if (buf[0] != '#') {
	if (ist) {
	  () = sscanf(buf, "%f %f", &x0, &x1);
	  e[n] = x0;
	  s[n] = x1;
	  n++;
	}
      }
    }
  }

  () = fclose(f);
  n = [0:n-1];
  e = e[n];
  s = s[n];
  return (e, s, et0);
}

define effcs(t, e, s) {
  variable nt, g, i, x, w, f;

  nt = length(t);
  g = Double_Type[nt];
  w = [0:length(e)-2];
  for (i = 0; i < nt; i++) {
    x = e/t[i];
    f = s[w]*(x[w+1]-x[w])*exp(-x[w]);
    g[i] = sum(f);
  }

  return g;
}

define osplot(x, y, c, x0, x1, y0, y1, ymin) {
  variable i;
  
  for (i = 0; i < length(x); i++) {
    color(c);
    if (x[i] >= x1 or x[i] < x0) continue;
    oplot([x[i], x[i]], [ymin, y[i]]);
  }
}

define splot(x, y, c, x0, x1, y0, y1, ymin) {
  variable i, d;
  
  if (x1 <= x0) {
    x0 = min(x);
    x1 = max(x);
    d = 0.1*(x1-x0);
    x0 -= d;
    x1 += d;
  } 
  if (y1 <= y0) {
    y0 = 0.0;
    y1 = max(y);
    d = 0.1*(y1-y0);
    y1 += d;
  }
  
  xrange(x0, x1);
  yrange(y0, y1);

  for (i = 0; i < length(x); i++) {
    color(c);
    if (x[i] >= x1 or x[i] < x0) continue;
    if (i == 0) plot([x[i], x[i]], [ymin, y[i]]);
    else oplot([x[i], x[i]], [ymin, y[i]]);
  }
}
 
