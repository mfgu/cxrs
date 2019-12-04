define gderiv(x, n) {
  variable y, y2, p;

  y = @x;
  y2 = y^2;
  if (n == 1) {
    p = -y;
  } else if (n == 2) {
    p = y2 - 1.0;
  } else if (n == 3) {
    p = y*(3 - y2);
  } else if (n == 4) {
    p = 3 + y2*(-6 + y2);
  } else if (n == 5) {
    p = y*(-15 + y2*(10 - y2));
  } else if (n == 6) {
    p = -15 + y2*(45 + y2*(-15 + y2));
  } else if (n == 7) {
    p = -y*(-105 + y2*(105 + y2*(-21 + y2)));
  } else if (n == 8) {
    p = 105 + y2*(-420 + y2*(210 + y2*(-28 + y2)));
  } else if ( n == 9) {
    p = -y*(945 + y2*(-1260 + y2*(378 + y2*(-36 + y2))));
  } else if (n == 10) {
    p = -945 + y2*(4725 + y2*(-3150 + y2*(630 + y2*(-45 + y2))));
  }
  
  return p*exp(-y2*0.5);
}

define gausswv(n, s, nd) {
  variable g, v, e, f, p, w;

  g = 1.77245385; % gamma(0.5);
  for (v = 1; v <= n; v++) {
    g *= (v-0.5); % gamma(n+0.5);
  }
  
  v = double([0:nd-1]);
  v[[nd/2+1:]] = v[[nd/2+1:]] - nd;
  v = 2.0*PI/nd * v;
  e = sqrt(2.0*PI*s);
  f = -((1j)^n)/g;
  p = s * v;
  w = e * f * p^n * exp(-p*p*0.5);
  return w;
}

define cwt(s, sc, x0, x1, dx) {
  variable x, y, yf, ns, nd, c, i, wf, w, nw, scw;

  if (x0 == NULL) {
    x0 = s.bin_lo[0];
  } 
  if (x1 == NULL) {
    x1 = s.bin_hi[-1];
  }

  scw = struct {bin_lo, bin_hi, value, err, thresh, scale,
		wt, ns, nd, imax, iridge, peaks};
  
  scw.scale = @sc;
  
  w = where(s.bin_lo >= x0 and s.bin_hi <= x1);
  if (dx == NULL) {
    scw.bin_lo = s.bin_lo[w];
    scw.bin_hi = s.bin_hi[w];
    scw.value = s.value[w];
    scw.err = s.err[w];
    y = s.value[w];
  } else {
    x = [double(w[0]):w[-1]+dx*0.1:dx];
    y = interpol(x, double(w), s.value[w]);
    scw.value = @y;
    scw.err = interpol(x, double(w), s.err[w]);
    scw.bin_lo = interpol(x, double(w), s.bin_lo[w]);
    scw.bin_hi = interpol(x, double(w), s.bin_hi[w]);
  }
  
  ns = length(sc);
  nd = length(w);
  scw.ns = ns;
  scw.nd = nd;

  if (nd mod 2) {
    nd++;
    y = [y, 0.0];
  }
  
  yf = fft(y, 1);
  
  c = Double_Type[ns,nd];
  for (i = 0; i < ns; i++) {
    wf = gausswv(2, sc[i], nd);
    c[i,*] = Real(fft(yf*Conj(wf), -1));
  }
  
  if (nd > scw.nd) {
    c = c[*,[0:scw.nd-1]];
  }

  scw.wt = c;

  return scw;
}

define localmax1d(y, ws, t) {
  variable n, ws2, m, i0, i1, w;

  if (ws < 3) ws = 3;
  ws2 = ws/2;
  n = length(y);
  m = Integer_Type[n];
  for (i0 = 0; i0 < n; i0 += ws2) {
    i1 = i0 + ws-1;
    if (i1 >= n) i1 = n-1;
    w = where(y[[i0:i1]] == max(y[[i0:i1]]));
    w = i0 + w[0];
    if (w > i0 and w < i1) {
      m[w] = 1;
    }
  }
  i0 = where(m > 0 and y < t);
  m[i0] = 0;

  w = where(m > 0);
  for (i0 = 1; i0 < length(w); i0++) {
    if (w[i0]-w[i0-1] < ws) {
      if (y[w[i0-1]] < y[w[i0]]) m[w[i0-1]] = 0;
      else {
	m[w[i0]] = 0;
	w[i0] = w[i0-1];
      }
    }
  }

  return m;
}
  
define localmax(s) {
  variable ns, nd, m, i;

  ns = s.ns;
  nd = s.nd;
  
  s.imax = Integer_Type[ns, nd];
  
  for (i = 0; i < ns; i++) {
    s.imax[i,*] = localmax1d(s.wt[i,*], int(s.scale[i]*2+2), s.thresh);
  }
}

define idridge(s, g) {
  variable r, ns, nd, nr, w0, w, nw0, nw, i, j, k, d, q, p, v, sc, m;

  sc = s.scale;
  m = s.imax;
  ns = s.ns;
  nd = s.nd;
  r = Integer_Type[ns,nd];

  nr = 0;
  for (i = ns-1; i >= 0; i--) {
    w = where(m[i,*] > 0);
    nw = length(w);
    if (nw > 0) {
      for (q = 0; q <= g; q++) {
	p = i+q+1;
	if (p < ns) {
	  w0 = where(m[p,*] > 0);
	  nw0 = length(w0);
	} else {
	  nw0 = 0;
	}
	v = where(r[i,w] == 0);
	if (length(v) == 0) break;
	for (j = 0; j < nw0; j++) {
	  k = where(r[i,w] == r[p,w0[j]]);
	  if (length(k) > 0) continue;
	  d = abs(w0[j] - w[v]);
	  k = where(d == min(d));
	  k = k[0];
	  if (d[k] < 2*(sc[i]+1)) {
	    r[i,w[v[k]]] = r[p,w0[j]];
	  }
	}
      }
      for (j = 0; j < nw; j++) {
	if (r[i,w[j]] == 0) {
	  nr++;
	  r[i,w[j]] = nr;
	}
      }
    }
  }
  
  s.iridge = r;
}
  
define idpeak(scw, snr, sr) {
  variable ir, r, x, sc, n, ns, nd, s, p, i, w, wi, wj, sw, xw, rw, k, j;

  ns = scw.ns;
  nd = scw.nd;
  sc = scw.scale;
  ir = scw.iridge;
  r = scw.wt;
  x = 0.5*(scw.bin_lo + scw.bin_hi);
  n = max(ir);
  s = Double_Type[5,n];  
  p = 0;
  for (i = 1; i <= n; i++) {
    w = where(ir == i);
    wi = w/nd;
    wj = w mod nd;    
    sw = sc[wi];
    xw = x[wj];
    rw = r[w];
    k = where(rw == max(rw));
    k = k[0];
    j = where(sw == min(sw));
    j = j[0];
    if (rw[k]/scw.err[wj[k]] < snr) continue;
    if (sr != NULL) {
      if (min(sw) > sr) continue;
    }
    s[0,p] = xw[j]; %line center
    s[1,p] = sw[int(k*0.9)]; %line width
    s[2,p] = rw[k]; %line peak count
    s[3,p] = wj[j]; %line center in channel
    s[4,p] = i; %ridge id.
    p++;
  }
  
  i = array_sort(s[0,[0:p-1]]);
  scw.peaks = s[*,i];
}

define wvpick(w, rg, snr, sr) {
  variable i, k;

  w.thresh = w.err;
  i = where(w.err <= 0);
  k = where(w.err > 0);
  w.thresh[k] = w.err[k]*snr*0.75;
  w.thresh[i] = median(w.err[k]^2/w.value[k])*snr^2;

  localmax(w);
  idridge(w, rg);
  idpeak(w, snr, sr);
}

% detect peaks with wavelet transformation
% s is the spectrum
% rg the maximum gap in the ridge identification, suggest 0.
% snr signal-noise-ratio of the peak, suggest 3.
% sr minimum scale the peak must be detected, suggest typical line width in pixel.
% x0, x1 the spectral range within which to search, if NULL, entire coerage.
% dx, the rebin factor of the sepctrum before search. < 1 oversample. NULL=1.
define wvpeak(s, rg, snr, sr, x0, x1, dx) {
  variable w, sc;

  sc = [1.0:20.1];
  if (dx != NULL) {
    sc /=dx;
    if (sr != NULL) sr /= dx;
  }

  w = cwt(s, sc, x0, x1, dx);
  wvpick(w, rg, snr, sr);

  return w;
}

% write the peaks to a file for mspec fit. 
% xw is the the line gap between groups (in pixels).
define writewvpeak(fn, s, xw) {
  variable f, i, w, j, k, p, q, x, m, sig;

  f = fopen(fn, "w");
  if (f == NULL) {
    vmessage("cannot open file %s", fn);
    return NULL;
  }
  
  w = median(s.peaks[1,*]);
  k = array_sort(s.peaks[0,*]);
  x = s.peaks[0,k];
  q = shift(x, 1) - x;
  sig = interpol(x, 0.5*(s.bin_hi+s.bin_lo), s.bin_hi-s.bin_lo);
  p = where(q > xw*sig or q < 0);
  for (i = 0; i < length(k); i++) {
    j = k[i];
    m = where(p >= i);
    m = m[0];
    () = fprintf(f, "%5d %3d %5d %5d %11.4E %11.4E %11.4E 1.0 0.0 0 0\n",
		 i, 1, m, int(s.peaks[4,j]), s.peaks[0,j], 
		 s.peaks[2,j]*s.peaks[1,j]*2, s.peaks[0,j]);
  }
  
  () = fclose(f);
}

% fit the peaks identified by wvpeak using mspec model.    
define fitwvpeak(s, lfn, rfn, pid) {
  variable id, j, w, x0, x1, m, xw, a, c, mj, mt, f, k, nw;

  set_fit_statistic("chisqr;sigma=lsq");

  delete_data(all_data());
  () = define_counts(s);
  
  mt = get_data_counts(1);
  mt.value[*] = 0.0;
  mt.err[*] = 0.0;

  read_response(rfn);
  read_lines(lfn, 0, 0);
  
  color(1);
  hplot(s);
  mj = max(_spec.ilines[0][2,*]);
  f = fopen(sprintf("flines%s.txt", pid), "w");
  id = 0;
  for (j = 0; j <= mj; j++) {
    w = where(_spec.ilines[0][2,*] == j);
    nw = length(w);
    _spec.ilines[0][[0:1],*] = 2;
    _spec.ilines[0][[0:1],w] = 0;
    x0 = min(_spec.lines[0][0,w]);
    x1 = max(_spec.lines[0][0,w]);
    if (j > 0) {
      m = where(_spec.ilines[0][2,*] == j-1);
      a = _spec.lines[0][0,m[-1]];
      a = 0.7*a + 0.3*x0;
      x0 -= 25.0*_spec.width[0];
      if (x0 < a) x0 = a;
    } else {
      x0 -= 25.0*_spec.width[0];
    }
    if (j < mj) {
      m = where(_spec.ilines[0][2,*] == j+1);
      a = _spec.lines[0][0,m[0]];
      a = 0.7*a + 0.3*x1;
      x1 += 25.0*_spec.width[0];
      if (x1 > a) x1 = a;
    } else {
      x1 += 25.0*_spec.width[0];
    }
    
    _spec.chan0 = 0.5*(x0+x1);
    xnotice(1, x0, x1);
    xw = interpol(0.5*(x0+x1), 0.5*(s.bin_lo+s.bin_hi), s.bin_hi-s.bin_lo)[0];
    set_abund(1, x0, x1);

    if (_spec.abund[0] > 0) {      
      _spec.lines[0][1,*] *= _spec.abund[0];
      _spec.abund[0] = 1.0;
    } else {
      _spec.abund[0] = 1.0;
    }
    set_mspec(xw, 100.0, 0, 0);

    set_fit_method("marquardt;tol=1E-8;max_loops=10000");    
    set_par(1, 0.01, 0, 0.005, 0.02);
    () = fit_counts();
    list_par;
    xw = dpoly(_spec.lines[0][0,w]-_spec.chan0, _spec.width)*2.35;
    a = dpoly(_spec.lines[0][0,w]-_spec.chan0, _spec.bkgd)*xw;
    for (k = 0; k < nw; k++) {
      c = _spec.lines[0][1,w[k]];
      if (c > 5.0) {	
	() = fprintf(f, "%5d %5d %5d %7.3f %10.3E %5.3f %10.3E\n",
		     id, j, _spec.ilines[0][3,w[k]],
		     _spec.lines[0][0,w[k]], c,
		     xw[k], a[k]);
	id++;
      }
    }
    m = get_model_counts(1);
    color(2);
    w = where(m.bin_lo >= x0 and m.bin_hi <= x1);
    ohplot(m.bin_lo[w], m.bin_hi[w], m.value[w]);
    mt.value[w] = m.value[w];
    color(1);    
  }
  () = fclose(f);
  whista(sprintf("msp%s.txt", pid), mt);
}
  
define picklines(s0, xr, a) {
  variable x0, y0, m, i;

  x0 = [s0.bin_lo[0],s0.bin_hi[-1]];
  y0 = [0.0,0.0];
  
  for (i = 0; i < length(xr)-1; i++) {
    m = wvpeak(s0, 0, a, 4.0, xr[i], xr[i+1], NULL);
    if (length(m.peaks) > 0) {
      x0 = [x0, m.peaks[0,*]];
      y0 = [y0, m.peaks[2,*]];
    } else {
      while (length(m.peaks) < 10) {
	a *= 0.75;
	m = wvpeak(s0, 0, a, 4.0, xr[i], xr[i+1], NULL);
      }
      x0 = [x0, m.peaks[0,*]];
      y0 = [y0, m.peaks[2,*]];
    }
  }
  i = array_sort(x0);
  x0 = x0[i];
  y0 = y0[i];
  return (x0, y0);
}
