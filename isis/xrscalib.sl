() = evalfile("ebitutil");
() = evalfile("rcols");
() = evalfile("minimize");
() = evalfile("conv");
  
define convsp(s, mc) {
  variable sc, a, w;

  sc = collapse(s, 0);  
  sc.err = (0.25 + sqrt(0.75+sc.value));
  if (mc > 0) {
    sc.value = conv(0.5*(sc.bin_lo+sc.bin_hi), sc.value, 
		    mc*(sc.bin_hi[0]-sc.bin_lo[0]));
    sc.value = abs(sc.value);
    %a = sqrt(mc*mc+3.5^2)/3.5;
    %sc.err = (0.25 + sqrt(0.75+a*sc.value))/a;
  }
  return sc;
}

define invcalib(x, c) {
  variable q, n, y, w;

  if (abs(c[2]) < 1E-10) {
    return (x-c[0])/c[1];
  }
  
  q = c[1]^2 - 4.0*(c[0]-x)*c[2];
  n = length(x);
  y = Double_Type[n];
  w = where(q < 0);
  if (length(w) > 0) {
    y[w] = 1E100;
  }
  w = where(q == 0);
  if (length(w) > 0) {
    y[w] = -c[1]/(2.0*c[2]);
  }
  w = where(q > 0);
  if (length(w) > 0) {
    q = sqrt(q[w]);
    y[w] = (-c[1] + q)/(2.0*c[2]);
  }

  if (typeof(x) != Array_Type) y = y[0];
  return y;
}

define apply1calib(pac) {
  variable c0, c1, c3, c4, c5, c6;
  variable i, np, p, w, nr, s0, x, xi, xj, j;
  variable k, v, nx, x0, x1, ck, i0, i1, dx;

  p = pac[0];
  x = pac[2];
  xi = pac[3];
  c0 = pac[4];
  c1 = pac[5];
  c3 = pac[6];
  c4 = pac[7];
  c5 = pac[8];
  c6 = pac[9];
  w = where(c3 == p);
  nr = length(w);
  s0 = rhista(pac[1]);
  nx = length(s0.bin_hi);
  for (j = 0; j < nr; j++) {
    k = w[j];
    ck = [c4[k], c5[k], c6[k]];
    x0 = c0[k];
    x1 = c1[k];
    v = where(x >= x0 and x < x1);
    if (length(v) > 0) {
      xi[v] = invcalib(x[v], ck);
    }
  }
  xj = x/xi;
  x0 = interpol(s0.bin_lo, xi, xj)*s0.bin_lo;
  x1 = interpol(s0.bin_hi, xi, xj)*s0.bin_hi;
  dx = x1-x0;
  w = where(dx > 0);
  i0 = w[0];
  if (i0 > 0) {
    w = [0:i0-1];
    dx[w] = 0;
    x0[w] = s0.bin_lo*x0[i0]/s0.bin_lo[i0];
    x1[w] = s0.bin_hi*x1[i0]/s0.bin_hi[i1];
  }
  w = where(dx < 0);
  if (length(w) > 0) {
    i1 = w[0];
    w = [i1:nx-1];
    x0[w] = s0.bin_lo*x0[i1-1]/s0.bin_lo[i1-1];
    x1[w] = s0.bin_hi*x1[i1-1]/s0.bin_hi[i1-1];
  } else {
    i1 = nx;
  }
  s0.value = agrebin(s0.bin_lo, s0.bin_hi, x0, x1,
		     s0.value, NULL);
  s0.err = sqrt(s0.value);
  whista(pac[1]+".csp", s0);
}

define applycalib(ipx, sfn, cfn) {
  variable c, np, x, i, xi, pac;
  
  c = rcols(cfn, [0,1,2,3,4,5,6], ["F","F","I","I","F","F","F"], 
	    '#', 0, -1, 1);
  np = length(ipx);
  x = [min(c[0]):max(c[1]):0.005];

  pac = List_Type[np];
  for (i = 0; i < np; i++) {
    xi = @x;
    pac[i] = {ipx[i], sfn[i], x, xi, c[0], c[1], c[3], c[4], c[5], c[6]};
  }
  parallel_map(Void_Type, &apply1calib, List_Type, pac);
}

define chpar(a0, s) {
  variable a0p, a;

  a = Double_Type[3];
  if (s.xr2 > 0) {
    a0p = a0[0]*s.xm0 - a0[0]*s.xr0 - a0[1]*(s.xm0/a0[0] - s.xr2)*s.xm0;
    a0p /= s.xm0 - a0[0]*s.xr2;
    a[0] = s.xr0 - a0p*s.xr2;
    a[1] = a0p - a0[1]*s.xr2;
    a[2] = a0[1];
  } else {
    if (abs(a0[1]) > 1E-16) {
      a0p = a0[0] - a0[1]*s.xm0/a0[0];
    } else {
      a0p = a0[0];
    }
    a[0] = 0.0;
    a[1] = a0p;
    a[2] = a0[1];
  }

  return (a, a0p);
}
  
define minimfun(p, s) {
  variable a, b, y0, z0, y1, x, y, w, w0, w1, n, x0, x1, lo, hi, v1;
  
  (a, ) = chpar(p, s);
  v1 = @(s.s1.value);
  w0 = where(s.s0.bin_lo >= s.xr0 and s.s0.bin_hi <= s.xr1);
  x0 = [s.s0.bin_lo[w0[0]], s.s0.bin_hi[w0[-1]]];
  if (s.xr2 > 0) {
    x1 = invcalib(x0, a);
  } else {
    x1 = invcalib(x0, a);
  }
  w1 = where(s.s1.bin_lo >= x1[0] and s.s1.bin_hi <= x1[1]);  
  n = length(w1);
  if (n < 3) return 1E30;
  x0 = [s.s1.bin_lo[w1], s.s1.bin_hi[w1[-1]]];
  x0 = dpoly(x0, a);
  lo = x0[[0:n-1]];
  hi = x0[[1:n]];
  w = where(hi <= lo);
  if (length(w) > 0) return 2E30;

  y0 = s.s0.value[w0];
  z0 = s.s0.err[w0];
  w = where(s.s0.bin_lo[w0] >= lo[0] and s.s0.bin_hi[w0] <= hi[-1]);
  y1 = Double_Type[length(y0)];
  y1[w] = rebin(s.s0.bin_lo[w0[w]], s.s0.bin_hi[w0[w]], lo, hi, v1[w1]);

  w1 = where(s.wn[w0]);  
  y0 = y0[w1];
  y1 = y1[w1];

  x = sum(y0);
  y = sum(y1);
  b = abs(x-y)/sqrt(x+y);
  if (b > 20.0) return 3E30;
  
  x = ((y0-y1)/z0[w1])^2;
  b = sum(x)/length(x);  
  b = b*sum(y0*y0)/sum(y0*y1);

  return b;
}

define align_linear(s) {
  variable a0, a, amin, amax, na, da, c, p, w, i;

  if (s.xr2 > 0) {
    a0 = s.xr0/s.xr2;
    amin = a0*0.9;
    amax = a0*1.1;
    c = 0.999*s.xm0/s.xr2;
    if (amax > c) amax = c;
    na = 250;
  } else {
    a0 = 1.0;
    amin = 0.6;
    amax = 1.4;
    na = 500;
  }
  da = (amax - amin)/(na-1.0);
  c = (da*s.xm0)/(s.s0.bin_hi[0]-s.s0.bin_lo[0]);
  vmessage("smoothing in align_linear: %f", c);
  if (c < 1.0) c = 1.0;
  s.s0 = convsp(s.sr0, c);
  s.s1 = convsp(s.sr1, c);
  a = [amin:amax+da*0.1:da];
  p = @a;
  for (i = 0; i < na; i++) {
    p[i] = minimfun([a[i],0.0], s);
  }
  c = min(p);
  w = where(p == c)[0];
  a0 = a[w];
  return (a0, a, p);
}

define find_minimum(a, b, p, s) {
  variable w, na, nb, i, j, id, a0;

  na = length(a);
  nb = length(b);
  w = where(p == min(p));
  w = w[0];
  i = w/nb;
  j = w mod nb;
  id = set_minim([a[i],b[j]], &minimfun, s);
  vmessage("Init Val: %2d %2d %10.3E %10.3E", i, j, a[i], b[j]);
  set_par(1, a[i], 0, a[i]-(a[1]-a[0]), a[i]+(a[1]-a[0]));
  set_par(2, b[j], 0);
  set_fit_method("subplex;tol=1E-6;maxnfe=10000");
  (a0, w) = get_minim();
  delete_data(id);
  set_fit_method("marquardt");
  set_fit_statistic("chisqr");  
  (a0, ) = chpar(a0, s);
  return (a0, w);
}

define align_minim(sr0, sr1, xr0, xr1, xr2, xri, nt, sig, sa, sb, na, nb, xd) {
  variable s, a0, b0, a0p, a, b, w, i, j, da, db, d, p, s0, s1, sc0;

  s0 = collapse(sr0, 0);
  s1 = collapse(sr1, 0);
  s = struct{sr0, sr1, s0, s1, xr0, xr1, xr2, xm0, xm1, wn};  
  s.wn = Integer_Type[length(sr0.bin_lo)];
  if (length(xri) > 1) {
    for (i = 0; i < length(xri); i+= 2) {
      s.wn = s.wn or (sr0.bin_lo >= xri[i] and sr0.bin_lo <= xri[i+1]);
    }
  }
  s.wn = not s.wn;
  s.sr0 = s0;
  s.sr1 = s1;
  s.xr0 = xr0;
  s.xr1 = xr1;
  s.xr2 = xr2;
  if (na <= 0) {
    na = 25;  
  }
  if (nb <= 0) {
    if (s.xr2 < 0) {
      nb = 41;
    } else {
      nb = 25;
    }
  }
  if (sa <= 0) {
    sa = 1.0;
  }
  if (sb <= 0) {
    sb = 1.0;
  }
  if (sig <= 3.5) {
    sig = 3.5;
  }
  sig = sqrt(sig*sig+nt*nt);
  
  sc0 = convsp(s0, nt);
  w = where(sc0.bin_lo >= xr0 and sc0.bin_hi <= xr1);
  b = max(sc0.value[w]);
  i = where(sc0.value[w] == b);
  i = w[i[0]];
  s.xm0 = 0.5*(sc0.bin_lo[i]+sc0.bin_hi[i]);  

  s.s0 = s0;
  s.s1 = s1;
  if (xd == 0.0) {
    (a0,,) = align_linear(s);
  } else {
    a0 = xd;
  }
  s.s0 = sc0;
  s.s1 = convsp(s1, nt);
  (b, a0p) = chpar([a0,0.0], s);
  if (s.xr2 > 0) {
    s.xm1 = invcalib(s.xm0, b);
  } else {
    s.xm1 = invcalib(s.xm0, b);
  }
  b0 = s.xm0/s.xm1;
  vmessage("a0: %10.3E %10.3E", a0, b0);
  vmessage("xr: %12.5E %12.5E %12.5E %10.3E %10.3E", 
	   s.xr0, s.xr1, s.xr2, s.xm0, s.xm1);
  da = sa*6.0*sig*(s0.bin_hi[0]-s0.bin_lo[0])*a0/s.xr1;
  d = s.xr0/s.xr1;
  if (d < 0.05) d = 0.05;
  db = sb*0.4*b0*d/s.xr1;
  p = db;
  if (s.xr2 > 0) {
    db *= s.xm0/(s.xr1-s.xr0);
  }
  d = da*2.0/(na-1.0);
  a = [a0-da:a0+da+0.1*d:d];
  if (nb > 1) {
    d = 2*db/(nb-1.0);
    b = [-db:db+0.1*d:d];
  } else {
    b = [0.0];
  }

  p = Double_Type[na,nb];
  for (i = 0; i < na; i++) {
    for (j = 0; j < nb; j++) {
      p[i,j] = minimfun([a[i],b[j]], s);
    }
  }

  (a0, w) = find_minimum(a, b, p, s);
  vmessage("Best Fit: %10.3E %10.3E %10.3E", w, get_par(1), get_par(2));

  return (a0, a, b, p, s);
}

% aling_pix
% s0, s1 are the two histograms.
% xr is an array of boundaries to break up s0 into segments.
% xri is an array of data ranges to be ignored in the fitting.
% returns
% c[nregion, 3], the polynomial coeff. to align s1 with s0 in all regions.
% x1, x0, the arrays map s1 scale to s0 scale.
define align_pix(sr0, sr1, xr, xri, nt0, ctr, xd) {
  variable np, c, r, n, t, i, nt, sig, sa, sb, na, nb;
  variable x0, x1, w, j, s0, s1, wn;

  sig = ctr[0];
  sa = 1.0;
  sb = 1.0;
  na = 0;
  nb = 0;
  if (length(ctr) == 2) {
    sa = ctr[1];
  }
  if (length(ctr) == 3) {
    na = ctr[2];
  }
  if (length(ctr) == 4) {
    sb = int(ctr[3]);
  }
  if (length(ctr) == 5) {
    nb = int(ctr[4]);
  }
  
  s0 = collapse(sr0, 0);
  s1 = collapse(sr1, 0);
  if (length(xri) > 1 and xri[0] <= 0.0) {
    w = where(s0.bin_lo >= xri[1]);
    c = sum(s0.value[w])/sum(s1.value[w]);
  } else {
    c = sum(s0.value)/sum(s1.value);
  }
  s1.value *= c;
  vmessage("spec 1 scaled: %10.3E", c);

  np = 3;
  c = Double_Type[length(xr)-1, np];
  r = Double_Type[np];
  n = length(xr);
  for (i = 1; i < n; i++) {
    if (i > length(nt0)) nt = nt0[-1];
    else nt = nt0[i-1];
    if (i > 1) {
      if (t > 0) {
	t = invcalib(xr[i-1], r);
      } else {
	t = invcalib(xr[i-1], r);
      }
    } else {
      t = -1.0;      
    }    
    (r,,,,) = align_minim(s0, s1, 
			  xr[i-1], xr[i], t, xri, nt, sig, sa, sb, na, nb, xd);
    c[i-1,*] = r;
    vmessage("Reg:%2d %10.3E %10.3E %12.5E, %12.5E, %12.5E", 
	     i-1, xr[i-1], xr[i], c[i-1,0], c[i-1,1], c[i-1,2]);
    () = fflush(stdout);
  }

  x0 = [sr0.bin_lo, sr0.bin_hi[-1]];
  x1 = @x0;
  for (i = 0; i < n-1; i++) {
    w = where(x0 >= xr[i] and x0 < xr[i+1]);
    if (length(w) > 0) {
      r = c[i,*];
      x1[w] = x0[w]/invcalib(x0[w], r);
    }
  }
  w = where(x0 < xr[0]);
  if (length(w) > 0) {
    x1[w] = xr[0]/invcalib(xr[0], c[0,*]);
  }
  w = where(x0 >= xr[-1]);
  if (length(w) > 0) {
    x1[w] = xr[-1]/invcalib(xr[-1], c[-1,*]);
  }
  t = x0/x1;
  nt = length(t);
  r = Integer_Type[nt];
  r[0] = 0;
  j = 1;
  for (i = 1; i < nt; i++) {
    if (t[i] > t[r[j-1]]) {
      r[j] = i;
      j++;
    }
  }
  r = r[[0:j-1]];
  x0 = x0[r];
  x1 = x1[r];

  return (c, x0, x1);
}

define combpoly(c0, c1, x0, x1, nx, np) {
  variable x, y, c, p, s;

  x = [x0:x1:(x1-x0)/(nx-1)];
  y = ipoly(x, [0.0,c1]);
  y = ipoly(y, [0.0,c0]);
  
  c = Double_Type[np];
  c[0] = c1[0]*c0[0];
  
  set_fit_method("marquardt;max_loops=1000;tol=1e-8;delta=1E-10");
  set_fit_statistic("chisqr;sigma=data");
  (p,s) = array_fit(y, x/y, NULL, c, NULL, NULL, &dpoly);
  
  return p;
}
  
% this function takes the coeff. of piece-wise polynomials returned by xrscalib
% in file cfn, and a list of voltage values for the reference pixel in file rfn,
% create a list of voltage for each pixel by inverting the piece=wise polynomial.
% it then fit the two set of voltages with a single polynomial with order np.
% return the coeff in file ofn.coef. the two list of valtages are output in 
% ofn%d.txt for each pixel. If rfn is NULL, then the list of reference voltages 
% is generated automatically from the limits given in cfn.
define xrspoly(ofn, cfn, rfn, np) {
  variable c, r, w, s, k, y, d, p, i, j, f;

  if (stat_file(cfn) == NULL) return -1;
  c = rcols(cfn, [0,1,2,3,4,5,6], 0, ' ', 0, -1, 1);
  c[2] = int(c[2]);
  c[3] = int(c[3]);
  if (rfn == NULL) {
    w = where(c[2] == c[2][0]);
    r = Double_Type[0];
    for (i = 0; i < length(w); i++) {
      y = c[0][w[i]] + 0.05;
      d = c[1][w[i]] - 0.05;
      if (d > y) {
	r = [r, [y:d:0.01]];
      }
    }
  } else {
    r = readcol(rfn, 5);
  }
  f = fopen(ofn+".coef", "w");
  if (f == NULL) {
    return -1;
  }
  y = @r;
  set_fit_method("marquardt;max_loops=1000;tol=1e-8;delta=1E-10");
  set_fit_statistic("chisqr;sigma=data");
  for (i = min(c[2]); i <= max(c[2]); i++) {
    w = where(c[2] == i);
    if (length(w) == 0) continue;
    for (j = 0; j < length(w); j++) {
      k = w[j];
      if (j == 0) {
	s = where(r < c[1][k]);
      } else {	
	s = where(r >= c[0][k] and r < c[1][k]);
      }
      if (length(s) > 0) {
	d = [c[4][k],c[5][k],c[6][k]];
	y[s] = invcalib(r[s], d);
      }      
    }
    d = Double_Type[np];
    k = w[0];
    d[0] = c[4][k];
    d[1] = c[5][k];
    d[2] = c[6][k];
    writecol(ofn+sprintf("%d.txt", c[3][k]), y, r);
    (p,s) = array_fit(y, r/y, NULL, d, NULL, NULL, &dpoly);
    for (j = 0; j < np; j++) {
      () = fprintf(f, "%2d %2d %2d %12.5E\n", i, c[3][w[0]], j, p[j]);
    }
  }

  () = fclose(f);
  return 0;
}

define align1pix(p) {
  variable a, x0, x1, s0, s1, lo, hi;

  if (p[0] == "") {     
    return;
  }
  vmessage(p[0]);
  s0 = p[1];
  s1 = p[2];
  () = fflush(stdout);
  (a,x0,x1) = align_pix(s0, s1, p[3], p[4], p[5], p[6], p[7]);
  () = fflush(stdout);
  lo = interpol(s1.bin_lo, x0/x1, x1)*s1.bin_lo;
  hi = interpol(s1.bin_hi, x0/x1, x1)*s1.bin_hi;
  s1.value = agrebin(s0.bin_lo, s0.bin_hi, lo, hi, s1.value, NULL);
  s1.bin_lo = @(s0.bin_lo);
  s1.bin_hi = @(s0.bin_hi);
  return {a,s1};
}

% xrscalib
% sr is the histogram of the reference pix. may be an integer which
%    picks out the sr-th spectrum in fspecs as the reference.
% ipx is the pixel ids of the specs. same length as fspecs.
% fspecs, file contains a list of file names for the spectra. may be
%    a list of file names.
% pref, the file name prefix used to construct output file names.
% fxr, file contains a list of region boundaries. may be a list.
% pr, the plotting parameters, pr[0] and pr[1] and the x-axis range,
%     pr[2] is the number of pixels to be plotted on one page,
%     pr[3] is the x-value where pixel names are written.
% rdsp_fun, is a function that reads the spectrum file and returns
%   a histogram. For two column data with mac line ending, use
%   machist in ebitutil.sl. for histograms written by whista, use
%   rhista.
% fxri, file contains the data ranges to be ignored in the fit, may be a list.
% ss, is the smoothing factor for the histograms
% ctr, additional control params
%
% xrscalib writes 3 files, "pref.coef" contains all polynomial coeff.
% "pref.ps" is the PS file of the plot with aligned spectra.
% "pref.tsp" is the histogram of the summed spectrum.
define xrscalib(sr, ipx, fspecs, pref, fxr, pr, rdsp_fun, 
		fxri, ss, ctr, xd) {
  variable f1, fn, id, ym, yd, yt, np, q, ip, a, m, s1, i;
  variable ip0, nsp, s0, label, specs, xr, j, npo, xri, xr0;

  npo = 3;
  if (typeof(fspecs) == String_Type) {
    specs = rcols(fspecs, [0], ["A"], ' ', 0, -1, 1);
    specs = specs[0];
  } else {
    specs = @fspecs;
  }
  if (typeof(fxr) == String_Type) {
    xr0 = readcol(fxr, 1);
  } else {
    xr0 = @fxr;
  }
  if (typeof(fxri) == String_Type) {
    xri = readcol(fxri, 1);
  } else {
    if (length(fxri) > 1) {
      xri = @fxri;
    } else {
      xri = fxri[0];
    }
  }

  xr0 = double(xr0);
  pr = double(pr);
  xri =  double(xri);
  
  ip0 = -1;
  if (typeof(sr) == Integer_Type) {
    ip0 = where(ipx == sr)[0];
    sr = (@rdsp_fun)(specs[ip0]);    
  }
  s0 = struct{bin_lo, bin_hi, value, err};
  s0.bin_lo = @sr.bin_lo;
  s0.bin_hi = @sr.bin_hi;
  s0.value = @sr.value;
  s0.err = @sr.err;

  f1 = fopen(sprintf("%s.coef", pref), "w");
  id = open_plot(sprintf("%s.ps/cps", pref));
  m = where(s0.bin_lo > pr[0] and s0.bin_hi < pr[1]);
  ym = 1.1*max(s0.value[m]);
  yd = 0.3*ym;
  xrange(pr[0], pr[1]);
  yrange(-0.02*ym, 0.95*ym+yd*(pr[2]));  
  np = int(pr[2]);
  color(2);
  set_line_width(1);
  hplot(sr);
  set_line_width(2);
  charsize(0.8);
  color(1);
  label = "reference: ";
  if (ip0 >= 0) label += specs[ip0];
  xylabel(pr[3], yd*0.3, label);
  q = 1;  
  nsp = length(specs);
  variable apps = List_Type[nsp];
  for (ip = 0; ip < nsp; ip++) {
    xr = @xr0;
    fn = specs[ip];
    m = stat_file(fn);
    if (m == NULL) {
      fn = "";
    } else if (ip == ip0) {
      whista(fn+".csp", sr);
      xr = @xr0;
      for (i = 0; i < length(xr)-1; i++) {
	() = fprintf(f1, "%12.5E %12.5E %2d %2d", xr[i], xr[i+1], ip, ipx[ip]);
	() = fprintf(f1, " %12.5E %12.5E %12.5E", 0.0, 1.0, 0.0);
	() = fprintf(f1, " %s\n", fn);
      }
      () = fprintf(f1, "\n");
      fn = "";
    }
    if (fn != "") {
      s1 = (@rdsp_fun)(fn);
    } else {
      s1 = NULL;
    }
    apps[ip] = {fn, sr, s1, xr, xri, ss, ctr, xd};
  }
  variable aprs = parallel_map(List_Type, &align1pix, apps);
  for (ip = 0; ip < nsp; ip++) {
    if (apps[ip][0] == "") continue;
    xr = apps[ip][3];
    fn = apps[ip][0];
    a = aprs[ip][0];
    s1 = aprs[ip][1];
    s0.value += s1.value;    
    whista(fn+".csp", s1);
    color(1);
    yt = yd+yd*((q-1) mod np);
    set_line_width(1);
    ohplot(s0.bin_lo, s0.bin_hi, s1.value+yt);
    color(2);
    set_line_width(2);
    charsize(0.8);
    xylabel(pr[3], yt+yd*0.3, fn);    
    for (i = 0; i < length(xr)-1; i++) {
      () = fprintf(f1, "%12.5E %12.5E %2d %2d", xr[i], xr[i+1], ip, ipx[ip]);
      for (j = 0; j < npo; j++) {
	() = fprintf(f1, " %12.5E", a[i,j]);
      }
      () = fprintf(f1, " %s\n", fn); 
    }
    () = fprintf(f1, "\n");
    if ((q mod np) == 0 or ip == nsp-1) {
      color(3);
      set_line_width(1);
      ohplot(s0.bin_lo, s0.bin_hi, s0.value/(q+1.0)+yt+yd);
      color(1);
      set_line_width(3);
      charsize(1.0);
      xylabel(pr[3], yt+yd+0.3*yd, "total");
      if (ip < nsp-1) {
	_pgpage();
	color(2);
	set_line_width(1);
	hplot(sr);
	color(1);
	set_line_width(2);
	charsize(0.8);
	xylabel(pr[3], yd*0.3, label);
      }
    }
    q++;
  }
  () = fclose(f1);
  close_plot(id);
  fn = sprintf("%s.tsp", pref);
  s0.err = sqrt(s0.value);
  whista(fn, s0);
}

define calib1d0(sr, ipx, ds, dr, xr0, pr0, xri0, ss, ctr) {
  variable fs, np, i, s0, s1, c, x0, x1, xr, pr, xri, nd, d, fs1, isr;

  isr = where(ipx == sr)[0];
  xr0 = double(xr0);
  xri0 = double(xri0);
  ctr = double(ctr);
  
  np = length(ipx);
  fs = String_Type[np];
  for (i = 0; i < np; i++) {
    fs[i] = sprintf("/vh%02d.txt", ipx[i]);
  }
  nd = length(ds);  

  for (i = 0; i < nd; i++) {
    d = ds[i];
    vmessage("processing dir %s", d);
    () = fflush(stdout);
    if (dr != NULL and dr != d) {
      s0 = rhista(dr + fs[isr]);
      s1 = rhista(d + fs[isr]);
      (c, x0, x1) = align_pix(s0, s1, [xr0[0],xr0[-1]], 0, ss, ctr, 1.0);
      xr = xr0/interpol(xr0, x0, x1);
      pr = pr0/interpol(pr0, x0, x1);
      pr[2] = pr0[2];
      if (length(xri0) > 1) {
	xri = xri0/interpol(xri0, x0, x1);
      } else {
	xri = xri0[0];
      }
    } else {
      xr = xr0;
      pr = pr0;
      xri = xri0;
    }
    xrscalib(sr, ipx, d+fs, d+"/calib1d", xr, pr, &rhista, xri, ss, ctr, 0);
    if (dr != NULL and dr != d) {
      fs1 = [dr, d] + "/calib1d.tsp";
      xrscalib(0, [0,1], fs1, d+"/calibxd", [xr0[0],xr0[-1]], pr0, &rhista, xri0, ss, ctr, 1.0);
      s0 = rhista(fs1[0]);
      s1 = rhista(d+"/calibxd.tsp");
      s1.value -= s0.value;
      s1.err = sqrt(s1.value);
      whista(d+"/calibxd.tsp", s1);
    }
  }
}

define calib1d(sr, ipx, ds, dr, xr, pr, xri, ss, ctr) {
  variable fs, np, i, s0, s1, c, x0, x1, cx;
  variable nd, d, fs1, isr, fsc, fsx, idr, ip, p;

  nd = length(ds);  
  if (nd <= 1 or dr == NULL) {
    calib1d0(sr, ipx, ds, dr, xr, pr, xri, ss, ctr);
    return;
  }
  isr = where(ipx == sr)[0];
  xr = double(xr);
  xri = double(xri);
  ctr = double(ctr);
  pr = double(pr);
  
  np = length(ipx);
  fs = String_Type[np];
  fsc = @fs;
  fsx = @fs;
  for (i = 0; i < np; i++) {
    fs[i] = sprintf("/vh%02d.txt", ipx[i]);
    fsc[i] = fs[i]+".csp";
  }
  idr = where(ds == dr)[0];
  vmessage("process refd %s", dr);
  xrscalib(sr, ipx, dr+fs, dr+"/calibrd", [xr[0],xr[-1]], pr, &rhista, xri, ss, ctr, 0);
  for (i = 0; i < nd; i++) {
    d = ds[i];
    vmessage("apply calib %s", d);
    applycalib(ipx, d+fs, dr+"/calibrd.coef");
  }
  for (p = 0; p < np; p++) {
    cx = sprintf("/calibxd%02d", ipx[p]);
    vmessage("process xd pix: %d", ipx[p]);
    fsx[p] = cx+".tsp";
    xrscalib(idr, [0:(nd-1)], ds+fsc[p], dr+cx, [xr[0],xr[-1]], pr, &rhista, xri,ss,ctr,0.0);
  }
  vmessage("process final refd %s", dr);
  xrscalib(sr, ipx, dr+fsx, dr+"/calib1d", xr, pr, &rhista, xri, ss,ctr,0.0);
}

define combsp(d0, d1, fs, md) {
  variable s0, s1, w, i, x0;
  if (md == 0) {
    s0 = machist(d0+"/"+fs);
    s1 = machist(d1+"/"+fs);
  } else {
    s0 = rhista(d0+"/"+fs);
    s1 = rhista(d1+"/"+fs);
  }
  w = where(s1.bin_lo > 5.0);
  i = where(s1.value[w] == max(s1.value[w]));
  x0 = (s1.bin_lo[w])[i[0]];
  w = where(s1.bin_lo > x0-0.25);
  s0.value[w] += s1.value[w];
  return s0;
}

define combnife(dlist, ipx, md) {
  variable r, i, p, s;
  r = rcols(dlist, [8,9], ["A", "A"], '#', 0, -1, 1);
  for (i = 0; i < length(r[0]); i++) {
    for (p = 0; p < length(ipx); p++) {
      if (md == 0) {
	s = combsp(r[0][i], r[1][i], sprintf("Volts_Pixel%d_Hist.txt", ipx[p]), md);
      } else {
	s = combsp(r[0][i], r[1][i], sprintf("vh%02d.txt", ipx[p]), md);
      }
      whista(sprintf(r[0][i]+"/vh%02d.txt", ipx[p]), s);
    }
  }
}

define addspec(pref, ds, ipx, cdir, cref, dr) {
  variable fn, i, j, nd, np, nc, s0, s, sfn, idc, rdir, idr, ddc;
  variable nr;

  idr = where(ds == dr)[0];
  if (cdir[0] == NULL or cdir[0] == "NULL" or cdir[0]=="null") {
    s = rhista(ds[idr]+"/calib1d.tsp");
    whista(pref+".tsp", s);
    return;
  }
  nd = length(ds);
  np = length(ipx);
  i = is_substr(cdir[0],"%");
  if (i > 0) {
    s = NULL;
    for (i = 0; i < nd; i++) {
      for (j = 0; j < np; j++) {
	sfn = sprintf(ds[i]+"/"+cdir[0], ipx[j]);
	vmessage(sfn);
	if (s == NULL) {
	  s = rhista(sfn);
	} else {
	  s0 = rhista(sfn);
	  s.value += s0.value;
	}
      }
    }
    s.err = sqrt(s.value);
    whista(pref+".tsp", s);
    return;
  }
  nc = length(cdir);
  nr = length(cref);
  idc = Array_Type[nc];
  ddc = String_Type[nc];
  for (i = 0; i < nc; i++) {
    if (stat_file(cdir[i]) != NULL) {
      idc[i] = (rcols(cdir[i], [0], ["I"], ' ', 0, -1, 1))[0];
      j = is_substr(cdir[i], ".");
      if (j > 1) {
	ddc[i] = cdir[i][[:j-2]];
      }
    } else {
      idc[i] = [0:(nd-1)];
      ddc[i] = cdir[i];
    }
    if (i < nr) {
      ddc[i] = ddc[i]+cref[i];
    } else {
      ddc[i] = ddc[i]+cref[nr-1];
    }
  }
  rdir = ddc[0];
  if (stat_file(rdir+"/calibrd.coef") != NULL) {
    for (i = 0; i < nd; i++) {
      vmessage("apply calibrd %s: %s", rdir, ds[i]);
      sfn = String_Type[np];
      for (j = 0; j < np; j++) {
	sfn[j] = sprintf("%s/vh%02d.txt", ds[i], ipx[j]);
      }
      applycalib(ipx, sfn, rdir+"/calibrd.coef");
    }
    for (j = 0; j < np; j++) {
      vmessage("apply calibxd %s: %d", rdir, ipx[j]);
      sfn = String_Type[nd];
      for (i = 0; i < nd; i++) {
	sfn[i] = sprintf("%s/vh%02d.txt.csp", ds[i], ipx[j]);
      }
      applycalib(idc[0], sfn, rdir+sprintf("/calibxd%02d.coef", ipx[j]));
      s = NULL;
      for (i = 0; i < nd; i++) {
	if (i == 0) s = rhista(sfn[i]+".csp");
	else {
	  s0 = rhista(sfn[i]+".csp");
	  s.value += s0.value;
	}
      }
      s.err = sqrt(s.value);
      whista(sprintf("%s/calibxd%02d.tsp", ds[idr], ipx[j]), s);
    }
    sfn = String_Type[np];
    for (j = 0; j < np; j++) {
      sfn[j] = sprintf("%s/calibxd%02d.tsp", ds[idr], ipx[j]);
    }
  } else {    
    sfn = String_Type[np];
    for (j = 0; j < np; j++) {
      sfn[j] = sprintf("%s/vh%02d.txt", ds[i], ipx[j]);
    }
  }
  vmessage("apply calib1d %s: %s", rdir, ds[idr]);
  applycalib(ipx, sfn, rdir+"/calib1d.coef");
  
  s = NULL;  
  for (j = 0; j < np; j++) {
    if (j == 0) {
      s = rhista(sfn[j]+".csp");
    } else {
      s0 = rhista(sfn[j]+".csp");
      s.value += s0.value;
    }
  }
  s.err = sqrt(s.value);
  whista(pref+".tsp", s);
}

define tabvolts(ofn, fn, dlist, ipx, cdir, cref, dr) {
  variable r, s, i, j, k, m, vfn, vf, d0, dt, d, idr, tdr, wr, u, q, t;
  variable c0, c1, c2, cx, w, nr, cfn, v0, v1, v2, v3, nw, wd, nwd, ip, np, nd;
  variable nc, nf, idc, ddc;

  t = is_substr(fn, ":");
  if (t > 0) {
    t = eval("["+fn+"]");
    m = length(t);
    r = Array_Type[5];
    r[0] = String_Type[m];
    r[1] = String_Type[m];
    r[2] = t;
    r[3] = @t;
    r[4] = @t;
    for (i = 0; i < m; i++) {
      r[0][i] = cdir[0];
      r[1][i] = sprintf("%d", i);
      r[4][i] = 1.0;
    }
  } else {    
    r = rcols(fn, [0,1,2,3,4], ["A","A","F","F","F"], '#', 0, -1, 1);
  }
  s = array_sort(r[0]);
  for (i = 0; i < length(r); i++) {
    r[i] = r[i][s];
  }
  nr = length(r[0]);
  nd = length(dlist);  
  np = length(ipx);
  nc = length(cdir);
  nf = length(cref);
  idc = Array_Type[nc];
  ddc = String_Type[nc];
  for (i = 0; i < nc; i++) {
    if (stat_file(cdir[i]) != NULL) {
      idc[i] = (rcols(cdir[i], [0], ["I"], ' ', 0, -1, 1))[0];
      j = is_substr(cdir[i], ".");
      if (j > 1) {
	ddc[i] = cdir[i][[:j-2]];
      }
    } else {
      idc[i] = [0:(nd-1)];
      ddc[i] = cdir[i];
    }
  }
  vfn = ofn + ".vts";
  vf = fopen(vfn, "w");
  idr = where(dlist == dr)[0];  
  d0 = "";
  
  for (j = 0; j < nr; j++) {
    if (r[0][j] != d0) {
      d0 = r[0][j];
      u = where(ddc == d0);
      if (length(u) == 1) {
	u = u[0];
	if (u < nf) tdr = cref[u];
	else tdr = cref[nf-1];
      } else {
	u = 0;
	tdr = cref[0];
      }
      d = d0+tdr;
      cfn = d+"/calib1d.coef";
      c0 = rcols(cfn, [0,1,2,3,4,5,6], ["F","F","I","I","F","F","F"], 
		 '#', 0, -1, 1);
      c1 = Array_Type[np];
      for (m = 0; m < np; m++) {
	cfn = d+sprintf("/calibxd%02d.coef", ipx[m]);
	if (stat_file(cfn) != NULL) {
	  c1[m] = rcols(cfn, [0,1,2,3,4,5,6], ["F","F","I","I","F","F","F"], 
			'#', 0, -1, 1);
	} else {
	  c1[m] = NULL;
	}
      }
      cfn = d+"/calibrd.coef";
      if (stat_file(cfn) != NULL) {
	c2 = rcols(cfn, [0,1,2,3,4,5,6], ["F","F","I","I","F","F","F"], 
		   '#', 0, -1, 1);
      } else {
	c2 = NULL;
      }
    }    
    v0 = r[3][j];
    for (m = 0; m < np; m++) {
      ip = ipx[m];	
      w = where(c0[0] <= v0 and v0 < c0[1] and c0[3] == ip);
      nw = length(w);
      if (nw == 0) continue;
      k = w[0];
      v1 = invcalib(v0, [c0[4][k], c0[5][k], c0[6][k]]);
      if (c2 == NULL) {
	dt = dlist[idr][[-5:]];
	() = fprintf(vf, "%2d %5s %2d %6s %6s %13.7E %13.7E %12.5E %13.7E\n",
		     0, dt, ip, r[0][j], r[1][j], r[2][j], v1, r[4][j], v0);
	continue;
      }
      cx = c1[m];
      if (cx == NULL) continue;
      for (i = 0; i < nd; i++) {
	dt = dlist[i][[-5:]];	
	wd = where(cx[0] <= v1 and v1 < cx[1] and cx[3]==idc[u][i]);
	nwd = length(wd);
	if (nwd != 1) continue;
	q = wd[0];
	v2 = invcalib(v1, [cx[4][q], cx[5][q], cx[6][q]]);
	wr = where(c2[0] <= v2 and v2 < c2[1] and c2[3] == ip);
	if (length(wr) != 1) continue;
	t = wr[0];
	v3 = invcalib(v2, [c2[4][t], c2[5][t], c2[6][t]]);
	() = fprintf(vf, "%2d %5s %2d %6s %6s %13.7E %13.7E %12.5E %13.7E\n",
		     i, dt, ip, r[0][j], r[1][j], r[2][j], v3, r[4][j], v0);
      }
    }
  }
  () = fclose(vf);
}

define strfmtime(tm) {
  return sprintf("%d/%d/%d-%d/%d/%d",
		 tm.tm_year+100, tm.tm_mon+1, tm.tm_mday, 
		 tm.tm_hour, tm.tm_min, tm.tm_sec);
}

define calibpoly(ofn, npo, dlist, md) {
  variable r, s, i, j, k, m, pfn, vfn, pf, d0, dt, d;
  variable c0, c1, w, nr, cfn, v0, v1, v2, np, id, x, y, mf;
  variable wf, wname, wpo, tm, npp;
  
  wpo = 6;
  pfn = ofn + ".cpc";
  vfn = ofn + ".vts";

  if (npo <= 1) {
    vmessage("Polynormial order <= 1, fitting disabled");
    () = fflush(stdout);
    return;
  }
  if (npo > wpo) {
    vmessage("Maximum polynomial order is %d", wpo);
    () = fflush(stdout);
    return;
  }
  if (npo > 1) {
    set_fit_method("marquardt;max_loops=10000;tol=1e-8;delta=1E-10");
    set_fit_statistic("chisqr;sigma=data");
    pf = fopen(pfn, "w");
    r = rcols(vfn, 1+[0,1,2,3,4,5,6,7], ["A","I","A","A","F","F","F","F"], 
	      '#', 0, -1, 1);
    
    id = open_plot(sprintf("%s.ps/cps", ofn));
    j = 0;
    np = 6;
    w = Integer_Type[np];
    w[*] = 1;
    multiplot(w);
    if (md == 0) {
      v0 = min(r[4]);
      v1 = max(r[4]);
    } else {
      v0 = min(r[7]);
      v1 = max(r[7]);
    }
    v2 = v1-v0;
    v0 -= 0.1*v2;
    v1 += 0.1*v2;
    xrange(v0, v1);
    for (i = 0; i < length(dlist); i++) {
      d0 = dlist[i];
      dt = d0[[-5:]];
      y = readcol(d0[[:-4]], 1);
      y = int(y[0] - 2082844800.0);
      tm = localtime(y);
      w = tm.tm_year;
      if (w > 100) w -= 100;
      wname = sprintf("Gain%02d%02d%02dT", w, tm.tm_mon+1, tm.tm_mday);
      wname += d0[[-2:]];
      wname += "_gain";
      wf = fopen(wname+".itx", "w");
      () = fprintf(wf, "IGOR\n");
      () = fprintf(wf, "WAVES/D/N=(%d,%d) %s\n", 32, wpo, wname);
      () = fprintf(wf, "BEGIN\n");      
      vmessage("fit polynomials for %s", dt);
      () = fflush(stdout);
      for (m = 0; m < 32; m++) {
	w = where(r[0] == dt and r[1] == m);
	if (length(w) == 0) {
	  () = fprintf(pf, "%2d %5s %2d %10.3E %10.3E", i, dt, m, 0.0, 0.0);
	  c1 = Double_Type[npo];
	  for (k = 0; k < npo; k++) {
	    () = fprintf(pf, " %13.6E", c1[k]);
	    () = fprintf(wf, " %13.6E", c1[k]);
	  }
	  () = fprintf(pf, "\n");	
	  for (; k < wpo; k++) {
	    () = fprintf(wf, " %3.1f", 0.0);
	  }
	  () = fprintf(wf, "\n");
	  continue;
	}
	
	if (length(w) < npo) {
	  npp = length(w);
	} else {
	  npp = npo;
	}
	if ((j mod np) == 0 and j > 0) {
	  _pgpage();
	}
	c0 = Double_Type[npp];
	v0 = r[5][w];
	if (md == 0) {
	  v1 = r[4][w]/v0;
	} else {
	  v1 = r[7][w]/v0;
	}
	v2 = r[6][w]/max(r[6][w]);
	w = array_sort(v0);
	v0 = v0[w];
	v1 = v1[w];
	v2 = v2[w];
	x = v0;
	y = v1;
	(c1, s) = array_fit(x, y, v2, c0, NULL, NULL, &dpoly);
	() = fflush(stdout);
	() = fprintf(pf, "%2d %5s %2d %10.3E %10.3E", i, dt, m, min(v0), max(v0));
	for (k = 0; k < npp; k++) {
	  () = fprintf(pf, " %13.6E", c1[k]);
	  () = fprintf(wf, " %13.6E", c1[k]);
	}
	for (; k < npo; k++) {
	  () = fprintf(pf, " %13.6E", 0.0);
	  () = fprintf(wf, " %13.6E", 0.0);
	}
	() = fprintf(pf, "\n");
	for (; k < wpo; k++) {
	  () = fprintf(wf, " %3.1f", 0.0);
	}
	() = fprintf(wf, "\n");
	v1 = v1*v0;
	v2 = v1-dpoly(v0, [0.0,c1]);
	mf = log10(max(abs(v2)));
	if (mf < 0) {
	  mf = 1 - int(mf);
	} else {
	  mf = -int(mf);
	}
	v2 *= 10^mf;
	pointstyle(-8);
	plot(v1, v2);
	pointstyle(-1);
	s = max(v2) - min(v2);
	v0 = max(v1)-0.5*(max(v1)-min(v1));
	v1 = max(v2)-0.2*s;
	xylabel(v0, v1, latex2pg(sprintf("%s pix %02d x 10^{%d}", dt, m, -mf)));
	v2 = moment(abs(v2));
	xylabel(v0, v1-0.15*s, sprintf("ave/max/min: %3.1f/%3.1f/%3.1f", 
				       v2.ave, v2.max, v2.min));
	j++;
      }
      () = fprintf(wf, "END\n");
      () = fprintf(wf, "X SetScale/P x 0,1,\"\", %s;", wname);
      () = fprintf(wf, "  SetScale/P y 0,1,\"\", %s;", wname);
      () = fprintf(wf, "  SetScale d 0,0,\"\", %s\n", wname);
      () = fprintf(wf, "X Note %s, \"GAINCAL:POLY;\\rVERSION:1;\\r", wname);
      () = fprintf(wf, "CREATETIME:%s;\\r", strfmtime("%Y%m%d-%H%M%S", tm));
      tm = localtime(_time());
      () = fprintf(wf, "MODTIME:%s\\r\"\n", strfmtime("%Y%m%d-%H%M%S", tm));
      () = fclose(wf);
    }
    close_plot(id);
    () = fclose(pf);
  }
}  

define mergepoly(fns, dlist) {
  variable wpo, npo, pfn, nfn, pf, wf, r, a, np, f, b, d0, dt, y;
  variable tm, w, i, j, k, m, wname, t, c;

  wpo = 6;
  pfn = fns[0] + ".cpc";
  nfn = length(fns)-1;
  pf = fopen(pfn, "w");

  r = Array_Type[nfn];
  np = Integer_Type[nfn];
  for (i = 0; i < nfn; i++) {
    f = fopen(fns[i+1]+".cpc", "r");
    if (-1 == fgets(&b, f)) continue;
    () = fclose(f);
    b = strchop(strcompress(b, " \t"), ' ', 0);
    np[i] = length(b)-5;
    f = String_Type[np[i]+5];
    f[0] = "I";
    f[1] = "A";
    f[2] = "I";
    f[[3:]] = "F";  
    a = rcols(fns[i+1]+".cpc", [0:length(b)-1], f, ' ', 0, -1, 1);
    r[i] = a;
  }
  npo = max(np);
  for (i = 0; i < length(dlist); i++) {
    d0 = dlist[i];
    dt = d0[[-5:]];
    y = readcol(d0[[:-4]], 1);
    y = int(y[0] - 2082844800.0);
    tm = localtime(y);
    w = tm.tm_year;
    if (w > 100) w -= 100;
    wname = sprintf("Gain%02d%02d%02dC", w, tm.tm_mon+1, tm.tm_mday);
    wname += d0[[-2:]];
    wname += "_gain";
    wf = fopen(wname+".itx", "w");
    for (m = 0; m < 32; m++) {
      c = Double_Type[npo+2];
      for (j = 0; j < nfn; j++) {
	w = where(r[j][1] == dt and r[j][2] == m);
	if (length(w) > 0) {
	  t = w[0];
	  if (r[j][3][t] != 0 or r[j][4][t] != 0) {
	    for (k = 0; k < np[j]+2; k++) {
	      c[k] = r[j][3+k][t];
	    }
	  }
	}
      }
      () = fprintf(pf, "%2d %5s %2d %10.3E %10.3E", i, dt, m, c[0], c[1]);
      for (k = 0; k < npo; k++) {
	() = fprintf(pf, " %13.6E", c[2+k]);
	() = fprintf(wf, " %13.6E", c[2+k]);
      }
      () = fprintf(pf, "\n");
      for (; k < wpo; k++) {
	() = fprintf(wf, " %3.1f", 0.0);
      }
      () = fprintf(wf, "\n");
    }
    () = fprintf(wf, "END\n");
    () = fprintf(wf, "X SetScale/P x 0,1,\"\", %s;", wname);
    () = fprintf(wf, "  SetScale/P y 0,1,\"\", %s;", wname);
    () = fprintf(wf, "  SetScale d 0,0,\"\", %s\n", wname);
    () = fprintf(wf, "X Note %s, \"GAINCAL:POLY;\\rVERSION:1;\\r", wname);
    () = fprintf(wf, "CREATETIME:%s;\\r", strfmtime("%Y%m%d-%H%M%S", tm));
    tm = localtime(_time());
    () = fprintf(wf, "MODTIME:%s\\r\"\n", strfmtime("%Y%m%d-%H%M%S", tm));
    () = fclose(wf);
  }
  () = fclose(pf);
}

define size_file(fn) {
  variable f, n;

  f = fopen(fn, "r");
  () = fseek(f, 0, SEEK_END);
  n = ftell(f);
  () = fclose(f);
  return n;
}

define convevts(ifn, ofn) {
  variable f, buf, n, w;
  
  f = fopen(ifn, "r");  
  if (f == NULL) {
    vmessage("cannot open file %s", ifn);
    () = fflush(stdout);
    return;
  }
  () = fseek(f, 0, SEEK_END);
  n = ftell(f);
  () = fseek(f, 0, SEEK_SET);
  () = fread(&buf, Char_Type, n, f);
  if (typeof(buf) == BString_Type) buf = bstring_to_array(buf);
  () = fclose(f);
  
  f = fopen(ofn, "w");  
  if (buf[0] != '#') {
    () = fprintf(f, "#  ");
  }
  w = where(buf == '\r');
  if (length(w) > 0) {
    buf[w] = '\n';
  }
  () = fwrite(buf, f);
  () = fclose(f);
}

define idcolumn(fn) {
  variable f0, h;

  f0 = fopen(fn, "r");
  () = fgets(&h, f0);
  h = strtrim(h);
  h = strchop(strcompress(h, " \t"), ' ', 0);
  if (length(h) < 2) {
    return String_Type[0];
  }
  h = h[[1:]];
  () = fclose(f0);

  return h;
}

define evtsfmt(fn) {
  variable f, b;
  
  f = fopen(fn, "r");
  () = fread(&b, Char_Type, 1, f);
  () = fclose(f);
  if (b == '#') return 1;
  else return 0;
}

variable filter_args = struct{r, c, w, w0, efn, volts, region, 
			      x0, y0, xs, ys, md, use_ph};
filter_args.use_ph = 0;
define filter_fun(f) {
  variable i, a, na, i0, i1, x0, y0, x, y;
  variable xc0, xc1, yc0, yc1, ch, w0, xr, yr, xs, ys;

  a = strchop(strcompress(f, " \t"), '?', 0);
  na = length(a);  
  if (na > 2) {
    vmessage("Invalid filter format. No filter applied.");
  } else {
    if (na == 2) {
      a[0] = strcompress(a[0], " \t");
      a[1] = strcompress(a[1], " \t");
      if (filter_args.use_ph == 1) {
	if (a[0] == "PH" or a[0] == "Ph") {
	  a[0] = "Volts";
	}
	if (a[1] == "PH" or a[1] == "Ph") {
	  a[1] = "Volts";
	}
      }
      i0 = where(filter_args.c == a[0]);
      if (length(i0) != 1) {
	vmessage("Invalid filter format. No filter applied.");
	return;
      }
      i0 = i0[0];
      i1 = where(filter_args.c == a[1]);
      if (length(i1) != 1)  {
	vmessage("Invalid filter format. No filter applied.");
	return;
      }
      i1 = i1[0];
      x = filter_args.r[i0];
      y = filter_args.r[i1];      
      if (a[0] == "Volts") x = filter_args.volts;
      if (a[1] == "Volts") y = filter_args.volts;

      if (filter_args.w0 != NULL) {
	w0 = where(filter_args.w0);
      } else {
	w0 = Integer_Type[length(x)];
      }
      if (filter_args.md < 2) {
	() = open_plot();      
	resize(20, 0.75);
	charsize(1.5);      
	_pgsvp(0.075, 0.925, 0.075, 0.9);
      }
      x0 = min(x[w0]);
      y0 = min(y[w0]);
      if (x0 < 0) x0 = 0;
      if (y0 < 0) y0 = 0;
      filter_args.x0 = x0;
      filter_args.y0 = y0;
      x -= x0;
      y -= y0;
      connect_points(0);
      pointstyle(-1);
      color(1);
      limits;
      if (a[0] == "Volts") {
	w0 = w0[where(x[w0] > 0)];
      } else if (a[1] == "Volts") {
	w0 = w0[where(y[w0] > 0)];
      }
      plot(x[w0], y[w0]);
      xr = [min(x[w0]), max(x[w0])];
      yr = [min(y[w0]), max(y[w0])];
      if (filter_args.md == 0) {
	filter_args.w = Integer_Type[length(x)];
      }
      xs = NULL;
      ys = NULL;      
      if (filter_args.md > 0) {
	xs = filter_args.xs;
	ys = filter_args.ys;
	connect_points(1);
	if (xs != NULL) {
	  for (i = 0; i < length(xs); i++) {
	    color(3);
	    if (length(xs[i]) > 0) {
	      oplot(xs[i], ys[i]);
	    }
	  }
	}
      }
      while (1) {
	_pgsci(2);
	() = _pgband(7, 0, 0, 0, &xc0, &yc0, &ch);
	i = zoom_plot(x[w0], y[w0], xr, yr, ch, xs, ys);
	if (i >= 0) {
	  _pgsci(2);
	  () = _pgband(7, 0, 0, 0, &xc0, &yc0, &ch);
	}
	if (ch == 'q' or ch == 'Q') break;
	if (filter_args.md > 0) continue;
	if (ch == 'n' or ch == 'N') {
	  color(1);
	  connect_points(0);
	  pointstyle(-1);
	  plot(x[w0], y[w0]);
	  filter_args.w[*] = 0;
	  xs = NULL;
	  ys = NULL;
	  vmessage("restart filter selection");
	  () = fflush(stdout);
	  continue;
	}
	_pgsci(2);
	() = _pgband(2, 0, xc0, yc0, &xc1, &yc1, &ch);
	if (ch == 'q' or ch == 'Q') break;
	if (ch == 'n' or ch == 'N') {
	  color(1);
	  connect_points(0);
	  pointstyle(-1);
	  plot(x[w0], y[w0]);
	  filter_args.w[*] = 0;
	  xs = NULL;
	  ys = NULL;
	  vmessage("restart filter selection");
	  () = fflush(stdout);
	  continue;
	}
	_pgsci(1);
	connect_points(1);
	if (xc0 > xc1) (xc0, xc1) = (xc1, xc0);
	if (yc0 > yc1) (yc0, yc1) = (yc1, yc0);
	if (xc0 < xr[0]) xc0 = xr[0];
	if (xc1 > xr[1]) xc1 = xr[1];
	if (yc0 < yr[0]) yc0 = yr[0];
	if (yc1 > yr[1]) yc1 = yr[1];
	if (filter_args.md == 0) {
	  filter_args.w = (filter_args.w or (x > xc0 and x <= xc1 and 
					     y > yc0 and y <= yc1));	
	  vmessage("Applying filter: \"%f<%s and %s<=%f and %f<%s and %s<=%f\"",
	  xc0+x0, a[0], a[0], xc1+x0, yc0+y0, a[1], a[1], yc1+y0);
	}
	() = fflush(stdout);
	color(3);
	oplot([xc0,xc1,xc1,xc0,xc0], [yc0,yc0,yc1,yc1,yc0]);
	if (xs == NULL) {
	  xs = [xc0, xc1, xc1, xc0, xc0];
	  ys = [yc0, yc0, yc1, yc1, yc0];
	} else {
	  xs = [xs, xc0, xc1, xc1, xc0, xc0];
	  ys = [ys, yc0, yc0, yc1, yc1, yc0];
	}
	connect_points(0);
	color(1);
      }
      if (filter_args.md <= 0) {
	if (xs == NULL) {
	  i0 = 0;
	} else {
	  i0 = length(xs);
	  xs += x0;
	  ys += y0;
	}      
	filter_args.region = Array_Type[i0/5];
	for (i = 0; i < i0; i += 5) {
	  filter_args.region[i/5] = [xs[i], xs[i+1], ys[i], ys[i+2]];
	}
	if (filter_args.md == 0) {
	  if (i0 == 0) {
	    filter_args.w[*] = 1;
	  }
	}
      }
    } else if (na == 1) {
      for (i = 0; i < length(filter_args.c); i++) {
	(f, ) = strreplace(f, filter_args.c[i], 
			   sprintf("filter_args.r[%d]", i), strlen(f));
      }
      vmessage("Applying filter: \"%s\"", f);
      () = fflush(stdout);
      eval(sprintf("filter_args.w = (%s);", f));
    }
  }
}

define readevts(ie, fn0, pr, ifg, filter, pc) {
  variable pix, flag, volts, ut, ep, ic0, w, s, q, a, nc, nw, fn, r, rc, i, j;
  
  if (ifg == NULL) ifg = 0x3;

  if (filter_args.efn != fn0) {
    filter_args.efn = fn0;  
    if (evtsfmt(fn0) == 0) {
      fn = fn0 + ".m2u";
      s = stat_file(fn);
      w = stat_file(fn0);
      if (s == NULL) {
	convevts(fn0, fn);
      } else {
	if (s.st_mtime <= w.st_mtime) {
	  convevts(fn0, fn);
	}
      }
      if (size_file(fn) <= size_file(fn0)) {
	convevts(fn0, fn);
      }
    } else {
      fn = fn0;
    }
    ic0 = idcolumn(fn);
    nc = length(ic0);
    if (nc < 4) {
      vmessage("%s does not have enough columns");
      return NULL;
    }
    r = Array_Type[nc];
    rc = Struct_Type[1+nc];
    for (i = 0; i <= nc; i++) {    
      rc[i] = struct {value};
      if (i == 0) {
	rc[i].value = fn;
      } else {
	rc[i].value = i;
      }
    }
    readcol(__push_args(rc));
    for (i = nc-1; i >= 0; i--) {
      r[i] = ();
    }
    filter_args.r = r;
    if (filter_args.use_ph) {
      w = where(ic0 == "Ph" or ic0 == "PH");
      i = where(ic0 == "Volts");
      if (length(w) == 1 and length(i) == 1) {
	w = w[0];
	i = i[0];
	(ic0[w],ic0[i]) = (ic0[i], ic0[w]);
      }
    }
    filter_args.c = ic0;
  } else {
    r = filter_args.r;
    ic0 = filter_args.c;
  }
  w = where(ic0 == "Pixel");
  if (length(w) == 0) {
    vmessage("no Pixel column in %s", fn0);
    () = fflush(stdout);
    return NULL;
  }
  pix = r[w[0]];
  w = where(ic0 == "Flags");
  if (length(w) == 0) {
    vmessage("no Flags column in %s", fn0);
    () = fflush(stdout);
    return NULL;
  }
  flag = r[w[0]];
  w = where(ic0 == "UTime");
  if (length(w) == 0) {
    vmessage("no UTime column in %s", fn0);
    () = fflush(stdout);
    return NULL;
  }
  ut = r[w[0]];
  if (pr != NULL) {
    w = where(ic0 == "EbitPhase");
    if (length(w) == 0) {
      vmessage("no EbitPhase column in %s", fn0);
      () = fflush(stdout);
      ep = @ut;
      ep[*] = 0.0;
    } else {
      ep = r[w[0]];
    }
  }

  w = where(ic0 == "Volts");
  if (length(w) == 0) {
    vmessage("no Volts column in %s", fn0);
    () = fflush(stdout);
    return NULL;
  }
  vmessage("%d", w[0]);
  volts = @(r[w[0]]);
  filter_args.w = NULL;
  filter_args.xs = NULL;
  filter_args.ys = NULL;
  if (pc != NULL) {
    (,w,) = array_info(pc);
    if (w == 1) {
      for (i = 0; i < 32; i++) {
	s = where(pix == i);
	if (length(s) == 0) continue;
	q = where(pc[0] == ie and pc[1] == i);
	if (length(q) == 0) continue;
	a = interpol(ut[s], 0.5*(pc[5][q]+pc[6][q]), pc[7][q]);
	w = where(a+1.0 != 1.0);
	volts[s[w]] = volts[s[w]]/a[w];
	w = where(a+1.0 == 1.0);
	volts[s[w]] = 0.0;
      }
    } else {
      for (i = 0; i < 32; i++) {
	s = where(pix == i and (volts <= pc[i,0]*0.95 or volts >= pc[i,1]*1.05));
	volts[s] = 0.0;
	s = where(pix == i and volts > pc[i,0]*0.95 and volts < pc[i,1]*1.05);
	if (length(s) == 0) continue;
	volts[s] = dpoly(volts[s], pc[i,[2:]])*volts[s];
      }
    }
  }
  filter_args.volts = volts;

  w = where(flag >= 64);
  if (length(w)>0) {
    flag[w] -= 64;
  }
  flag = int(flag) & 0x3F;
  pix = int(pix);
  w = not (flag & (~ifg));
  if (length(where(w)) == 0) {
    vmessage("no good evts after flag filtering");
    () = fflush(stdout);
  }
  if (pr != NULL) {
    s = 1;
    if (length(pr) == 1) {
      s = ep >= pr[0];
    } else {
      if (pr[1] > pr[0]+1E-10) {
	vmessage("%f %f", pr[0], pr[1]);
	s = ep >= pr[0] and ep <= pr[1];
      }
    }
    w = w and s;
  }
  filter_args.w0 = w;
  if (length(where(w)) == 0) {
    vmessage("no good evts after flag and ebitphase filtering");
    () = fflush(stdout);
  } else {
    if (filter != NULL) {
      if (filter != "") {
	if (filter[[0:0]] != "#") {
	  filter_fun(filter);
	}
      }
    }
    if (filter_args.w != NULL) {
      w = w and filter_args.w;
    }
  }
  w = where(w);
  nc = length(volts);
  nw = length(w);
  vmessage("event file %s: total events %d, events after filtering %d", 
	   fn0, nc, nw);
  () = fflush(stdout);
  pix = pix[w];
  flag = flag[w];
  volts = volts[w];
  ut = ut[w];
  return (pix, flag, volts, ut);
}

% reads event files from a list of file names in evtf. make histograms for
% different time intervals.
% fn, the output file name giving a list a directories and time intervals for
%     the histograms.
% ief, a list of file indexes for the event files.
% dv, the volts bin size.
% vmin, minimum of the volts for the histogram.
% vmax, maximum of the volts for the histogram.
% ip, the pixels to include in the histograms.
define splitevts(ief, listf, ofn, dv, vmin, vmax, ip, mc, fpc) {
  variable f, vlo, vhi, m, i, j, nw, w, nip, dir, i0, nc, ep, cts, pc;
  variable pix, flag, volts, ut, tlo, thi, nt, s, h, sfn, ie, p, evtf;

  dv = double(dv);
  vmin = double(vmin);
  vmax = double(vmax);

  p = rcols(listf, [0:5], ["A","I","F","F","I","R"], '#', 0, -1, 1);
  if (ip == NULL) ip = [0:31];
  vlo = [vmin:vmax:dv];
  vhi = vlo + dv;
  evtf = String_Type[length(ief)];
  for (m = 0; m < length(evtf); m++) {
    ie = ief[m];
    evtf[m] = p[0][ie];
    vmessage("processing event file %s", evtf[m]);
    () = fflush(stdout);
    pc = fpc;
    if (pc != NULL) {
      pc += sprintf(".%02d", ie);
      if (stat_file(pc) != NULL) {
	nc = 8;
	f = String_Type[nc];
	f[*] = "F";
	f[0] = "I";
	f[1] = "I";
	f[4] = "I";
	pc = rcols(pc, [0:nc-1], f, ' ', 0, -1, 1);
      } else {
	pc = NULL;
      }
    }
    f = fopen(sprintf("%s%02d", ofn, ie), "w");
    filter_args.md = 0;
    (pix, flag, volts, ut) = readevts(m, evtf[m], 
				      [p[2][ie],p[3][ie]], p[4][ie],
				      p[5][ie], pc);
    nip = length(ip);
    nt = 0;
    cts = p[1][ie];
    for (i = 0; i < nip; i++) {
      w = where(pix == ip[i]);
      nw = length(w);
      vmessage("pixel %d has %d events", ip[i], nw);
      () = fflush(stdout);
      if (nw == 0) continue;
      if (nt == 0) {
	if (cts < 0) {
	  nt = nw/abs(cts);
	  if (nt == 0) nt = 1;
	} else {
	  nt = cts;
	}
	i0 = i;
	j = [0:nw-1:(nw/nt)];
	tlo = ut[w[j]];
	tlo[0] = min(ut);
	if (nt > 1) {
	  thi = make_hi_grid(tlo);
	  thi[-1] = max(ut);
	} else {
	  thi = [max(ut)];
	}
      }
      for (j = 0; j < nt; j++) {
	dir = sprintf("%s%02dt%02d", ofn, ie, j);
	if (i == i0) {
	  if (stat_file(dir) == NULL) {
	    () = system(sprintf("mkdir %s", dir));
	  } else {
	    () = system(sprintf("rm -rf %s/vh*.txt", dir));
	  }
	  () = fprintf(f, "%16.5f %16.5f %s\n", tlo[j], thi[j], dir);
	}
	s = where(ut[w] >= tlo[j] and ut[w] < thi[j]);
	if (length(s) > 0) {
	  h = histogram(volts[w[s]], vlo, vhi);
	} else {
	  h = histogram([vlo,vhi], vlo, vhi);
	}
	sfn = sprintf("%s/vh%02d.txt", dir, ip[i]);
	writecol(sfn, vlo, vhi, h, sqrt(h));      
      }
    }
    () = fclose(f);
    for (j = 0; j < nt; j++) {
      sfn = sprintf("%s%02dt%02d/vh.ps", ofn, ie, j);
      h = open_plot(sprintf("%s/cps", sfn));
      charsize(1.3);
      limits;
      multiplot([1,1,1,1]);
      cts = 0;
      for (i = 0; i < nip; i++) {
	sfn = sprintf("%s%02dt%02d/vh%02d.txt", ofn, ie, j, ip[i]);
	if (stat_file(sfn) == NULL) continue;
	if ((cts mod 4) == 0 and cts > 0) {
	  _pgpage();
	}
	s = collapse(rhista(sfn), mc);
	color(1);
	hplot(s);
	color(2);
	xylabel(vmax-0.2*(vmax-vmin), 0.8*max(s.value), 
		sprintf("pix %02d", ip[i]));
	cts++;
      }
      close_plot(h);
    }
  }
}

define wrtevts(evtf, m) {
  variable f, fmt, args, w, i, n, c;

  c = @(filter_args.c);
  if (filter_args.use_ph) {
    w = where(c == "Volts");
    i = where(c == "PH" or c == "Ph");
    w = w[0];
    i = i[0];
    (c[w], c[i]) = (c[i], c[w]);
  }
  vmessage("writting event file %s", evtf);
  () = fflush(stdout);
  n = length(filter_args.c);
  if (m == 0) {
    f = fopen(evtf, "w");
    () = fprintf(f, "#");
    for (i = 0; i < n; i++) {
      () = fprintf(f, " %s", c[i]);
    }
  } else {
    f = fopen(evtf, "a+");
  }
  fmt = "%d %9.6f %13.5f %13.10f %d %16.5f";
  if (n == 7) {
    fmt += " %12.5f\n";
  } else {
    fmt += "\n";
  }
  args = Struct_Type[n];
  w = filter_args.w0;
  if (filter_args.w != NULL) w = w and filter_args.w;
  w = where(w);
  for (i = 0; i < n; i++) {
    args[i] = struct {value};
    if (filter_args.c[i] == "Pixel" or 
	filter_args.c[i] == "Flags") {
      args[i].value = int(filter_args.r[i][w]);
    } else if (filter_args.c[i] == "Volts") {
      args[i].value = filter_args.volts[w];
    } else {
      args[i].value = filter_args.r[i][w];
    }
  }
  () = fprintf(f, "\n");
  () = array_map(Int_Type, &fprintf, f, fmt, __push_args(args));
  () = fclose(f);
}

define mevent(ief, listf, ofn) {
  variable p, evtf, m, ie, pix, flag, volts, ut, f;

  p = rcols(listf, [0:5], ["A","I","F","F","I","R"], '#', 0, -1, 1);
  evtf = String_Type[length(ief)];
  for (m = 0; m < length(evtf); m++) {
    ie = ief[m];
    evtf[m] = p[0][ie];
    (pix, flag, volts, ut) = readevts(m, evtf[m], 
				      [p[2][ie], p[3][ie]], p[4][ie], 
				      p[5][ie], NULL);
    wrtevts(ofn, m);
  }
}
    
define applyd(ief, listf, fpc, idc, ift, iwt, edc) {
  variable p, m, ie, evtf, pc, f, pix, flag, volts, ut, nc, i, w, s;

  p = rcols(listf, [0:5], ["A","I","F","F","I","R"], '#', 0, -1, 1);
  for (m = 0; m < length(ief); m++) {
    ie = ief[m];
    evtf = p[0][ie];
    vmessage("processing event file %s", evtf);
    () = fflush(stdout);
    if (idc == 0) {
      pc = NULL;
    } else {
      pc = fpc;
      pc += sprintf(".%02d", ie);
      if (stat_file(pc) == NULL) {
	vmessage("no dcc file for evt file %d: %s", ie, evtf);
	() = fflush(stdout);
	continue;
      } else {
	nc = 8;
	f = String_Type[nc];
	f[*] = "F";
	f[0] = "I";
	f[1] = "I";
	f[4] = "I";
	pc = rcols(pc, [0:nc-1], f, ' ', 0, -1, 1);
      }
    }
    filter_args.md = 0;
    if (ift) {
      (pix, flag, volts, ut) = readevts(m, evtf, [p[2][ie],p[3][ie]],
					p[4][ie], p[5][ie], pc);
    } else {
      (pix, flag, volts, ut) = readevts(m, evtf, NULL, p[4][ie], NULL, pc);
    }
    if (iwt) {
      wrtevts(evtf+edc, 0);
    }
    filter_args.md == 1;
    filter_args.xs = NULL;
    filter_args.ys = NULL;
    filter_fun("UTime ? Volts");
  }
}

define dcorr_region(fdc0, ie) {
  variable f, nr, v0, v1, bt0, bt1, xy, i, v, fdc;

  fdc = sprintf("%s.%02d", fdc0, ie);
  nr = length(filter_args.region);
  if (nr == 0) {
    if (stat_file(fdc) == NULL) {
      if (stat_file(fdc0) == NULL) {
	v0 = Double_Type[0];
	v1 = v0;
	bt0 = v0;
	bt1 = v0;
      } else {
	(xy, v0, v1) = readcol(fdc0, 1, 3, 4);
	i = where(xy == 0);
	bt0 = v0[i];
	bt1 = v1[i];
	i = where(xy == 1);
	v0 = v0[i];
	v1 = v1[i];
      }
    } else {
      (xy, v0, v1) = readcol(fdc, 1, 3, 4);
      i = where(xy == 0);
      bt0 = v0[i];
      bt1 = v1[i];
      i = where(xy == 1);
      v0 = v0[i];
      v1 = v1[i];
    }
  } else {
    v = -1.0;
    v0 = Double_Type[0];
    v1 = Double_Type[0];
    for (i = 0; i < nr; i++) {
      xy = filter_args.region[i];
      if (xy[2] <= v) break;
      v0 = [v0, xy[2]];
      v1 = [v1, xy[3]];
      v = v0[-1];
    }
    if (i == nr) {
      xy = filter_args.region[0];
      bt0 = [xy[0]];
      bt1 = [xy[1]];
    } else {
      bt0 = Double_Type[0];
      bt1 = Double_Type[0];
      for (; i < nr; i++) {
	xy = filter_args.region[i];
	bt0 = [bt0, xy[0]];
	bt1 = [bt1, xy[1]];
      }
    }
    f = fopen(fdc, "w");
    for (i = 0; i < length(bt0); i++) {
      () = fprintf(f, "0 %2d %15.6f %15.6f\n", i, bt0[i], bt1[i]);
    }
    for (i = 0; i < length(v0); i++) {
      () = fprintf(f, "1 %2d %8.5f %8.5f\n", i, v0[i], v1[i]);
    }
    () = fclose(f);
  }
  return (v0, v1, bt0, bt1);
}

define dcorrv(nd, na, ndp, v0, v1, bt0, bt1, pc) {
  variable nt, t0, t1, w, ut, volts, volts0, pix, nip, dt, nit, c, t, x;
  variable vdc, vr, nw, sv, n, i, j, q, nc, w0, it0, it1, nv, niw, vip, y;
 
  t0 = Double_Type[0];
  t1 = Double_Type[0];
  it0 = @t0;
  it1 = @t1;
  w0 = where(filter_args.w0);
  w = where(filter_args.c == "UTime");
  w = w[0];
  ut = filter_args.r[w][w0];
  volts = filter_args.volts[w0];
  w = where(filter_args.c == "Pixel");
  w = w[0];
  pix = filter_args.r[w][w0];
  w = where(filter_args.c == "Volts");
  w = w[0];
  volts0 = filter_args.r[w][w0];
  w0 = (volts >= v0[0] and volts <= v1[0]);
  nv = length(v0);
  for (i = 1; i < nv; i++) {
    w0 = w0 or (volts >= v0[i] and volts <= v1[i]);
  }
  w0 = where(w0);
  ut = ut[w0];
  volts = volts[w0];
  pix = pix[w0];
  volts0 = volts0[w0];
  w0 = (ut >= bt0[0] and ut <= bt1[0]);
  for (i = 1; i < length(bt0); i++) {
    w0 = w0 or (ut >= bt0[i] and ut <= bt1[i]);
  }
  w0 = where(w0);
  ut = ut[w0];
  volts = volts[w0];
  pix = pix[w0];
  volts0 = volts0[w0];
  w0 = Integer_Type[length(pix)];
  for (i = 0; i < length(pc[*,0]); i++) {
    if (pc[i,1] <= pc[i,0]) {
      w0 = w0 or (pix == i);
    }
  }
  w0 = where(not w0);
  ut = ut[w0];
  volts = volts[w0];
  pix = pix[w0];
  volts0 = volts0[w0];
  
  nip = length(pc[*,0]);
  nw = Integer_Type[nip];
  for (i = 0; i < length(bt0); i++) {
    for (j = 0; j < nip; j++) {
      if (pc[j,1] <= pc[j,0]) nw[j] = 0;
      else {
	w = where(pix == j and ut >= bt0[i] and ut <= bt1[i]);
	nw[j] = length(w);
      }
    }
    nc = sum(nw);
    n = int(nc/nd);
    if (n == 0) n++;
    dt = (bt1[i] - bt0[i])/n;
    for (q = 0; q < n; q++) {
      t0 = [t0, bt0[i] + q*dt];
      t1 = [t1, t0[-1]+dt];
    }
    nc = min(nw[where(nw>0)]);
    if (na > 0) {
      n = int(nc/na);
      if (n == 0) n++;
      dt = (bt1[i] - bt0[i])/n;
      for (q = 0; q < n; q++) {
	it0 = [it0, bt0[i] + q*dt];
	it1 = [it1, it0[-1]+dt];
      }
    }
  }
  if (na == 0) {
    it0 = [bt0[0]];
    it1 = [bt1[-1]];
    n = 1;
  }
  x = 0.5*(t0+t1);
  x = [t0[0], x, t1[-1]];
  nt = length(t0);
  nit = length(it0);
  vdc = Double_Type[nv, nip, nt];
  vip = Double_Type[nv, nip, nit];
  nw = Integer_Type[nv, nt];
  niw = Integer_Type[nv, nip, nit];
  for (i = 0; i < nv; i++) {    
    vr = Double_Type[nt];
    for (q = 0; q < nt; q++) {
      vmessage("AllPix: %2d %7.4f %7.4f %04d/%04d %16.5f %16.5f", 
	       i, v0[i], v1[i], q, nt, t0[q], t1[q]);
      () = fflush(stdout);
      w = where(ut >= t0[q] and ut <= t1[q] and
		volts >= v0[i] and volts <= v1[i]);      
      nw[i,q] = length(w);
      if (nw[i,q] > 0) {
	vr = sum(volts[w])/nw[i,q];
      }
      for (j = 0; j < nip; j++) {
	if (pc[j,1] <= pc[j,0]) vdc[i,j,q] = 0.0;
	else {
	  vdc[i,j,q] = ipoly(vr, [0.0,pc[j,[2:]]]);
	}
      }      
    }
    for (q = 0; q < nit; q++) {
      vmessage("OnePix: %2d %7.4f %7.4f %04d/%04d %16.5f %16.5f", 
	       i, v0[i], v1[i], q, nit, it0[q], it1[q]);
      () = fflush(stdout);
      for (j = 0; j < nip; j++) {
	w = where(pix == j and ut >= it0[q] and ut <= it1[q] and
		  volts >= v0[i] and volts <= v1[i]);
	niw[i,j,q] = length(w);
	if (niw[i,j,q] == 0) {
	  vip[i,j,q] = 0.0;
	} else {
	  y = vdc[i,j,*];
	  y = [vdc[i,j,0], y, vdc[i,j,-1]];
	  t = volts0[w]/interpol_points(ut[w], x, y);
	  vip[i,j,q] = sum(t);
	}
      }
    }
    for (j = 0; j < nip; j++) {
      sv = vip[i,j,*];
      nc = niw[i,j,*];
      w = sum(nc);
      if (w > 0) {
	vip[i,j,*] /= sum(sv)/w;
      }
    }
  }
  vr = Double_Type[nip, nit];
  for (j = 0; j < nip; j++) {
    for (q = 0; q < nit; q++) {
      sv = vip[*,j,q];
      nc = niw[*,j,q];
      w = sum(nc);
      if (w > 0) {
	vr[j,q] = sum(sv)/w;
      } else {
	vr[j,q] = 0.0;
      }
    }
  }
  for (q = 0; q < nt; q++) {
    t = 0.5*(t0[q]+t1[q]);
    w = where(it0 <= t and it1 >= t);
    w = w[0];
    for (j = 0; j < nip; j++) {
      vdc[*,j,q] *= vr[j,w];
    }
  }

  vip = @vdc;
  for (i = 0; i < nv; i++) {
    for (j = 0; j < nip; j++) {
      vr = mean(vdc[i,j,*]);
      if (1.0+vr != 1.0) {
	vip[i,j,*] = vdc[i,j,*]/vr;
      } else {
	vip[i,j,*] = 0.0;
      }
    }
  }
  niw = Double_Type[nip, nt, ndp];
  for (j = 0; j < nip; j++) {
    for (q = 0; q < nt; q++) {
      sv = _reshape(vip[*,j,q], [nv]);
      i = where(vdc[*,j,q] > 0);
      if (length(i) == nv) {	
	if (ndp > 1) {
	  nc = _reshape(vdc[*,j,q], [nv]);
	  (c, ) = array_fit(nc, sv, sqrt(nw[*,q]), Double_Type[ndp], 
			    NULL, NULL, &dpoly);
	  niw[j,q,*] = c;
	} else {
	  nc = _reshape(nw[*,q], [nv]);
	  niw[j,q,0] = sum(sv*nc)/sum(nc);
	}
      }
    }
  }
  return (t0, t1, vdc, niw);
}

define dcorr(fdc, ief, listf, pfn, iad) {
  variable f, b, np, pc, r, p, m, ie, evtf, w, i, j, v, v0, v1, c, nd, na;
  variable nw, vdc, t0, t1, bt0, bt1, x, y, pix, flag, volts, ut;
  variable nv, vr, q, im, ndp, nip, nt, xs, ys;

  if (length(iad) != 3) {
    vmessage("last argument must contain 3 ints");
    () = fflush(stdout);
    return;
  }
  nd = iad[0];
  na = iad[1];
  ndp = 1;
  im = iad[2];
  if (nd <= 0 or na < 0) {
    vmessage("counts must be > 0 in the last argument");
    () = fflush(stdout);
    return;
  }

  f = fopen(pfn, "r");
  if (-1 == fgets(&b, f)) return;
  () = fclose(f);
  b = strchop(strcompress(b, " \t"), ' ', 0);
  np = length(b)-5;  
  f = String_Type[np+5];
  f[0] = "I";
  f[1] = "A";
  f[2] = "I";
  f[[3:]] = "F";
  nip = 32;
  pc = Double_Type[nip,np+2];
  r = rcols(pfn, [0:length(b)-1], f, ' ', 0, -1, 1);  
  p = rcols(listf, [0:5], ["A","I","F","F","I","R"], '#', 0, -1, 1);
  for (m = 0; m < length(ief); m++) {
    ie = ief[m];
    evtf = p[0][ie];
    vmessage("processing evtfile %s", evtf);
    () = fflush(stdout);
    w = where(r[1] == sprintf("%02dt00", ie));
    if (length(w) != 32) {
      vmessage("no cpc for %s, skipping.", evtf);
      () = fflush(stdout);
      continue;
    }
    for (i = 0; i < length(w); i++) {
      for (j = 0; j < np+2; j++) {
	pc[i,j] = r[j+3][w[i]];
      }
    }
    filter_args.md = -1;
    (pix, flag, volts, ut) = readevts(m, evtf, 
				      [p[2][ie],p[3][ie]], p[4][ie],
				      "UTime ? Volts", pc);
    (v0, v1, bt0, bt1) = dcorr_region(fdc+sprintf(".dcr%d", im), ie);
    nv = length(v0);
    if (nv == 0) {
      vmessage("no regions to dcorr");
      () = fflush(stdout);
      continue;
    }
    f = fopen(fdc+sprintf(".dcc%d.%02d", im, ie), "w");
    c = ndp;
    if (c > nv) c = nv;
    (t0, t1, vdc, nw) = dcorrv(nd, na, ndp, v0, v1, bt0, bt1, pc);
    nt = length(t0);
    i = where(filter_args.c == "UTime");
    i = i[0];    
    color(1);
    connect_points(0);
    (q,j,x,y) = _pgqwin();
    xrange(q, j);
    yrange(x, y);
    w = where(filter_args.w0);
    plot(filter_args.r[i][w]-filter_args.x0, filter_args.volts[w]-filter_args.y0);
    filter_args.xs = Array_Type[nip*nv];
    filter_args.ys = Array_Type[nip*nv];
    for (i = 0; i < nip; i++) {
      for (j = 0; j < nt; j++) {
	if (pc[i,1] <= pc[i,0]) {
	  x = 0.0;
	  y = 0.0;
	} else {
	  x = ipoly(v0[0], [0.0,pc[i,[2:]]]);
	  y = ipoly(v1[-1], [0.0,pc[i,[2:]]]);
	}
	() = fprintf(f, "%2d %2d %7.4f %7.4f %5d %16.5f %16.5f %9.7f",
		     ie, i, x, y, j, t0[j], t1[j], nw[i,j,0]);
	for (q = 1; q < c; q++) {	  
	  () = fprintf(f, " %12.5f",  nw[i,j,q]);
	}	
	() = fprintf(f, " %-s\n", evtf);
      }
      () = fflush(f);
      for (q = 0; q < nv; q++) {
	vr = vdc[q,i,0];
	if (1+vr != 1) {
	  j = (i mod 14)+2;
	  connect_points(1);
	  color(j);	
	  x = 0.5*(t0+t1) - filter_args.x0;
	  y = dpoly(vdc[q,i,*], pc[i,[2:]])*vdc[q,i,*] - filter_args.y0;
	  oplot(x, y);
	  filter_args.xs[q*nip+i] = x;
	  filter_args.ys[q*nip+i] = y;
	} else {
	  filter_args.xs[q*nip+i] = Double_Type[0];
	  filter_args.ys[q*nip+i] = Double_Type[0];
	}
      }      
    }
    () = fclose(f);
    filter_args.md = 2;
    filter_fun("UTime ? Volts");
  }
}

define merged(ief, ofn, mfn) {
  variable n, m, i, j, q, a, pc, w, s, f, ie, fn, nc, ip0, ip1, x, b;

  n = length(mfn);
  for (m = 0; m < length(ief); m++) {
    ie = ief[m];
    vmessage("merging evtfile %2d", ie);
    () = fflush(stdout);
    pc = NULL;
    for (i = 0; i < n; i++) {
      fn = mfn[i] + sprintf(".%02d", ie);
      if (stat_file(fn) != NULL) {
	nc = 9;
	f = String_Type[nc];
	f[*] = "F";
	f[0] = "I";
	f[1] = "I";
	f[4] = "I";
	f[-1] = "A";
	if (pc == NULL) {
	  pc = rcols(fn, [0:nc-1], f, ' ', 0, -1, 1);
	  pc[4][*] = i;
	} else {
	  a = rcols(fn, [0:nc-1], f, ' ', 0, -1, 1);
	  a[4][*] = i;
	  for (j = 0; j < nc; j++) {
	    pc[j] = [pc[j], a[j]];
	  }
	}
      }
    }
    if (pc == NULL) {
      vmessage("no files to merge for evtfile %2d", ie);
      () = fflush(stdout);
      continue;
    }
    f = fopen(sprintf("%s.%02d", ofn, ie), "w");
    ip0 = min(pc[1]);
    ip1 = max(pc[1]);
    for (i = ip0; i <= ip1; i++) {
      w = where(pc[1] == i);
      if (length(w) == 0) continue;
      a = Array_Type[nc];
      for (j = 0; j < nc; j++) a[j] = pc[j][w];
      x = 0.5*(a[5]+a[6]);
      s = array_sort(x);
      for (j = 0; j < nc; j++) a[j] = a[j][s];      
      for (q = 1; q < length(a[0]); q++) {
	if (a[4][q] != a[4][q-1]) {
	  w = where(a[4] == a[4][q]);
	  s = where(a[4] == a[4][q-1]);
	  if (q < s[-1]) {
	    b = interpol_points(x[q], x[s], a[7][s]);
	  } else {
	    b = a[7][q-1];
	  }
	  a[7][w] *= b/a[7][q];
	  a[4][w] = a[4][q-1];
	} 
      }
      for (q = 0; q < length(a[0]); q++) {
	() = fprintf(f, "%2d %2d %7.4f %7.4f %5d %16.5f %16.5f %9.7f",
		     a[0][q], a[1][q], a[2][q], a[3][q], q, 
		     a[5][q], a[6][q], a[7][q]);
	() = fprintf(f, " %-s\n", a[8][q]);
      }
    }
    () = fclose(f);
  }
}

define prepln(hdir, lnf, x0, x1) {
  variable s, xr, yr, w, i, fn;

  x0 = double(x0);
  x1 = double(x1);
  s = rhista(hdir+".tsp");
  xr = [x0, x1];
  w = where(s.bin_lo > x0 and s.bin_hi < x1);
  yr = [min(s.value[w]), max(s.value[w])];
  yr = [yr[0]-0.1*(yr[1]-yr[0]), yr[1]+0.1*(yr[1]-yr[0])];
  xrange(x0, x1);  
  hplot(s);
  resize(25, 0.5);
  _pgsvp(0.075, 0.925, 0.075, 0.9);
  charsize(1.5);
  for (i = 0; i < length(lnf); i++) {
    fn = hdir+lnf[i]+".txt";
    title("prepare lines for " + fn);
    hplot(s);
    vmessage("Prepare lines for " + fn);
    () = fflush(stdout);
    prepare_lines(fn, s, xr, yr, 1, 0);
  }
}

define fitsp(hdir, lnf, a, e, x0, x1, w, p, alpha, nb) {
  variable s, x, y, c, i, j, y0, y1, fn1, fn2, f, r;

  if (p != NULL) {
    f = fopen(p, "r");
    if (-1 == fgets(&y, f)) return;
    () = fclose(f);
    y = strchop(strcompress(y, " \t"), ' ', 0);
    s = length(y)-5;
    f = String_Type[s+5];
    f[0] = "I";
    f[1] = "A";
    f[2] = "I";
    f[[3:]] = "F";
    x = rcols(p, [0:s+4], f, ' ', 0, -1, 1);
    i = where(x[2] == w and x[1] == x1);
    if (length(i) == 0) return;
    i = i[0];
    p = Double_Type[s];
    for (j = 0; j < s; j++) {
      p[j] = x[5+j][i];
    }
    f = fopen(a+".lno", "w");
    () = fprintf(f, "Abund  1  1.0E5\n");
    c = 0;
    for (i = 0; i < length(lnf); i++) {
      fn1 = hdir + lnf[i] + ".txt";
      r = rcols(fn1, [0,1,2,3,4,5,6], ["I","I","I","I","F","F","F"], 
		'#', 0, -1, 1);
      r[6] = dpoly(r[4], [0.0, p]);
      for (j = 0; j < length(r[0]); j++) {
	() = fprintf(f, "%3d %2d %4d %4d %13.7E %11.4E %13.7E 1.0 0.0 0 0\n",
		     c, 1, 0, c, r[4][j], 1.0, r[6][j]);
	c++;
      }
    }
    () = fclose(f);
    return;
  }
  if (strlow(a) == "volts") {
    a = hdir + lnf[0] + ".txt";
    fn1 = NULL;
    fn2 = NULL;
  } else {
    i = where(const->AtomicSymbol == a);
    if (length(i) == 0) {
      vmessage("Error: Atomic symbol %s is invalid");
      () = fflush(stdout);
      exit(1);
    }
    if (length(lnf) != 2) {
      vmessage("Error: FPrep must have two line files for H- and He-like series");
      () = fflush(stdout);
      exit(1);
    }
    fn1 = hdir+lnf[0]+".txt";
    fn2 = hdir+lnf[1]+".txt";
  }
  x0 = double(x0);
  x1 = double(x1);
  s = rhista(hdir+".tsp");
  i = where(s.bin_lo > 0);
  s.bin_lo = s.bin_lo[i];
  s.bin_hi = s.bin_hi[i];
  s.value = s.value[i];
  s.err = s.err[i];
  xrange(x0, x1);
  color(1);
  hplot(s);
  resize(25, 0.5);
  charsize(1.5);
  _pgsvp(0.075, 0.925, 0.075, 0.9);
  hplot(s);
  if (1+w == 1) {
    _pgsci(2);    
    vmessage("Choose Gaussian FWHM");
    () = fflush(stdout);
    color(2);
    (x0, x1, y0, y1) = _pgqwin();
    xylabel(x0+0.1*(x1-x0), y1+0.05*(y1-y0), "Choose Gaussian FWHM");
    color(1);
    _pgsci(2);
    () = _pgband(7, 0, 0, 0, &x, &y, &c);
    i = zoom_plot(s, NULL, [x0,x1], [y0,y1], c, NULL, NULL);
    if (i >= 0) {
      color(2);
      (x0, x1, y0, y1) = _pgqwin();
      xylabel(x0+0.1*(x1-x0), y1+0.05*(y1-y0), "Choose Gaussian FWHM");
      _pgsci(2);
      () = _pgband(7, 0, 0, 0, &x, &y, &c);
    }
    w = x;
    _pgsci(2);
    () = _pgband(4, 0, x, y, &x, &y, &c);
    if (x > w) {
      w = (x-w)/2.35;
    } else {
      w = (w-x)/2.35;
    }
  } else {
    w = w/2.35;
  }
  hdir = hdir + e;
  xrsrsp(hdir+".rsp", w, alpha, nb);
  fit_kspec(hdir, s, a, fn1, fn2, hdir+".rsp", 2);
}

define pltsp(fn) {
  variable s, xr, yr, xc0, yc0, ch, i, w, h, s0, xr1, yr1, xs, ys, xt, yt;

  if (stat_file(fn) == NULL) {
    vmessage("no spec file %s", fn);
    () = fflush(stdout);
    return;
  }
  s0 = rhista(fn);
  s = collapse(s0, 1);
  () = fflush(stdout);
  color(1);
  () = open_plot();
  resize(25, 0.5);
  charsize(1.5);
  _pgsvp(0.075, 0.935, 0.075, 0.9);
  title(fn);
  limits;
  hplot(s);
  xr = [0.0, 0.0];
  yr = [0.0, 0.0];
  xr1 = @xr;
  yr1 = @yr;
  (xr[0], xr[1], yr[0], yr[1]) = _pgqwin();
  xs = NULL;
  while (1) {
    _pgsci(2);
    () = _pgband(0, 0, 0, 0, &xc0, &yc0, &ch);
    if (ch == 'q' or ch == 'Q') break;    
    if (ch == 's' or ch == 'S') {
      resize(0.0, 100.0);
      (, , ,h) = _pgqvsz(2);
      resize(0.0, 0.01);
      (, w, ,) = _pgqvsz(2);
      resize(w*0.1, h/w);
      hplot(s);
      xs = NULL;
    }
    if (ch >= '1' and ch <= '9') {
      i = ch - '0';
      vmessage("collapsing the spectrum by %d", i);
      () = fflush(stdout);
      s = collapse(s0, i);
      limits;      
      hplot(s);
      (xr[0], xr[1], yr[0], yr[1]) = _pgqwin();
      xs = NULL;
    }
    i = zoom_plot(s, NULL, xr, yr, ch, NULL, NULL);
    if (i < 0) {
      vmessage("x = %12.5E y = %12.5E", xc0, yc0);
      () = fflush(stdout);
      (xr1[0], xr1[1], yr1[0], yr1[1]) = _pgqwin();
      xt = xr1[1] - 0.3*(xr1[1]-xr1[0]);
      yt = yr1[1] + 0.035*(yr1[1]-yr1[0]);
      if (xs != NULL) {
	color(0);
	xylabel(xt, yt, sprintf("(%g, %g)", xs, ys));
      }
      color(1);
      xylabel(xt, yt, sprintf("(%g, %g)", xc0, yc0));
      xs = xc0;
      ys = yc0;
    } else {
      xs = NULL;
    }
  }
}

define mergeln(ofn, mfn, lab, cnts) {  
  comblines(ofn, mfn+".lno", lab, cnts);
}

define apply1gain(pag) {
  variable i, m, ds, sd, st, r, dt, w, t, c, j, s, sfn, lo, hi;
  variable xlo, xhi, b, np;
  
  i = pag[0];
  m = pag[1];
  ds = pag[2];
  r = pag[3];
  xlo = pag[4];
  xhi = pag[5];
  np = pag[6];
  vmessage("processing %s pix %d %d", ds[i], m, length(xlo));
  () = fflush(stdout);
  dt = ds[i][[-5:]];
  w = where(r[1] == dt and r[2] == m);
  if (length(w) == 0) {
    w = where(r[1] == "00t00" and r[2] == m);
  }
  if (length(w) == 0) {
    return;
  }
  t = w[0];
  c = Double_Type[np];
  for (j = 0; j < np; j++) c[j] = r[5+j][t];
  w = where(c != 0);
  if (length(w) == 0) {
    return;
  }
  sfn = sprintf("%s/vh%02d.txt", ds[i], m);
  if (stat_file(sfn) == NULL) {
    return;
  }
  s = rhista(sfn);
  w = where(s.bin_lo >= r[3][t] and s.bin_hi <= r[4][t]);      
  lo = dpoly(s.bin_lo, c)*s.bin_lo;
  hi = dpoly(s.bin_hi, c)*s.bin_hi;

  if (xlo == NULL) {
    for (j = w[0]; j >= 0; j--) {
      if (hi[j] <= lo[j] or lo[j] <= 0) break;
    }
    t = j+1;
    for (j = w[-1]; j < length(lo); j++) {
      if (hi[j] <= lo[j]) break;
    }
    j--;
    b = mean(hi[w]-lo[w]);
    xlo = [lo[t]:hi[j]:b];
    xhi = xlo + b;
  }
  s.value = agrebin(xlo, xhi, lo, hi, s.value, w);
  s.err = sqrt(s.value);
  s.bin_lo = @xlo;
  s.bin_hi = @xhi;
  return {s};
}

define applygain(pfn, ds, ofn, de, emin, emax, ipx) {
  variable f, b, r, rc, np, nd, npx, ndp, i, j, k, pag;
  variable sfn, xlo, xhi, st, sd, s, s1;
  
  f = fopen(pfn, "r");
  if (-1 == fgets(&b, f)) return;
  () = fclose(f);
  b = strchop(strcompress(b, " \t"), ' ', 0);
  np = length(b)-5;  
  f = String_Type[np+5];
  f[0] = "I";
  f[1] = "A";
  f[2] = "I";
  f[[3:]] = "F";
  
  r = rcols(pfn, [0:length(b)-1], f, ' ', 0, -1, 1);
  nd = length(ds);
  npx = length(ipx);
  st = NULL;
  if (de <= 0.0) {
    xlo = NULL;
    xhi = NULL;  
  } else {
    xlo = [emin:emax:de];
    xhi = xlo + de;
  }
  ndp = nd*npx;
  k = 0;
  pag = List_Type[ndp];
  rc = {};
  for (i = 0; i < length(r); i++) {
    list_append(rc, r[i]);
  }
  for (i = 0; i < nd; i++) {
    for (j = 0; j < npx; j++) {
      pag[k] = {i,ipx[j],ds,rc,xlo,xhi,np};
      k++;
    }
  }
  s = List_Type[ndp];
  for (k = 0; k < ndp; k++) {
    s[k] = apply1gain(pag[k]);
    if (length(s[k]) > 0) {
      k++;
      break;
    }
  }
  s1 = s[k-1][0];
  if (xlo == NULL) {
    xlo = @s1.bin_lo;
    xhi = @s1.bin_hi;
    for (i = k; i < ndp; i++) {
      pag[i][4] = xlo;
      pag[i][5] = xhi;
    }
  }
  if (k < ndp) {
    s1 = parallel_map(List_Type, &apply1gain, pag[[k:ndp-1]]);
    for (i = k; i < ndp; i++) {
      s[i] = s1[i-k];
    }
  }
  s1 = s[k-1][0];
  st = collapse(s1, 1);
  sd = collapse(s1, 1);
  st.value = 0.0;
  k = 0;
  for (i = 0; i < nd; i++) {
    sd.value = 0.0;
    for (j = 0; j < npx; j++) {
      if (length(s[k]) == 0) continue;
      s1 = s[k][0];
      sfn = sprintf("%s/ch%02d.txt", ds[i], ipx[j]);      
      whista(sfn, s1);
      sd.value = sd.value + s1.value;
      k++;
    }
    st.value = st.value + sd.value;
    sd.err = sqrt(sd.value);
    sfn = sprintf("%s/ch.tsp", ds[i]);
    whista(sfn, sd);
  }
  st.err = sqrt(st.value);
  whista(ofn, st);
}

define cevent(ief, listf, pfn, ef) {
  variable p, m, i, j, w, t, c, ie, evtf, f, b, np, r;
  variable pix, flag, volts, ut, ip;

  p = rcols(listf, [0:5], ["A","I","F","F","I","R"], '#', 0, -1, 1);
  f = fopen(pfn, "r");
  if (-1 == fgets(&b, f)) return;
  () = fclose(f);
  b = strchop(strcompress(b, " \t"), ' ', 0);
  np = length(b)-5;  
    f = String_Type[np+5];
  f[0] = "I";
  f[1] = "A";
  f[2] = "I";
  f[[3:]] = "F";
  r = rcols(pfn, [0:length(b)-1], f, ' ', 0, -1, 1);
    
  for (m = 0; m < length(ief); m++) {
    ie = ief[m];
    evtf = p[0][ie];
    vmessage("processing event file %s", evtf);
    () = fflush(stdout);
    filter_args.md = 0;
    (pix, flag, volts, ut) = readevts(m, evtf, [p[2][ie],p[3][ie]],
				      p[4][ie], p[5][ie], NULL);
    ip = where(filter_args.c == "Pixel");
    ip = ip[0];
    for (i = 0; i < 32; i++) {
      w = where(r[1] == sprintf("%02dt00", ie) and r[2] == i);
      if (length(w) == 0) continue;
      t = w[0];
      c = Double_Type[np];
      for (j = 0; j < np; j++) c[j] = r[5+j][t];
      w = where(c != 0);
      if (length(w) == 0) continue;
      vmessage("Pixel %2d", i);
      () = fflush(stdout);
      w = where(filter_args.r[ip] == i);
      if (length(w) == 0) continue;
      filter_args.volts[w] = dpoly(filter_args.volts[w], c)*filter_args.volts[w];
    }    
    wrtevts(evtf+ef, 0);
  }
}
