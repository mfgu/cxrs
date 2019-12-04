() = evalfile("ebitutil");

define conv(x, y, sig) {
  variable a, yk, yc, n, xd, yd;
  
  if (_NARGS == 4) a = ();
  else a = 0.0;

  n = length(x);
  xd = Double_Type[n*2];
  yd = @xd;
  xd[[n:]] = x-x[0];
  xd[[0:n-1]] = -reverse(x-x[0])-(x[1]-x[0]);
  yd[[n:]] = y;
  
  if (a == 0.0) {
    yk = dgauss(xd, [1.0, 0.0, sig]);
  } else {
    yk = dvoigt(xd, [1.0, 0.0, sig, a]);
  }
  yc = Real(fft(fft(yd, 1)*Conj(fft(yk, 1)), -1));
  yc = yc[[0:n-1]];
  yc *= sum(y)/sum(yc);
  return yc;
}

define convd(e, s, sig) {
  variable i, x, y, x0, x1, dx, sm, w, r, n, a;

  if (_NARGS == 4) a = ();
  else a = 0.0;
  if (typeof(sig) == Array_Type) {
    sm = 5.0*max(sig);
    x0 = min(e)-sm;
    x1 = max(e)+sm;
    sm = min(sig);
    dx = 0.1*sm;
    x = [x0:x1:dx];
    n = length(x);
    y = Double_Type[n];
    for (i = 0; i < length(e); i++) {    
      if (a == 0.0) {
	r = (e[i]-x)/sig[i];
	r *= r;
	w = where(r < 15.0);
	y[w] += (1.0/(sqrt(2.0*PI)*sig[i]))*s[i]*exp(-0.5*((e[i]-x[w])/sig[i])^2); 
      } else {
	r = (e[i]-x)/sig[i];
	r = abs(r);
	w = where(r < sqrt(25.0+a*1e4));
	y += (1.0/(1.41421*sig[i]))*s[i]*voigt(a, (x-e[i])/(1.41421*sig[i]));
      }
    }
  } else {
    sm = 5.0*sig;
    x0 = min(e)-sm;
    x1 = max(e)+sm;
    sm = min(sig);
    dx = 0.1*sm;
    x = [x0:x1:dx];
    n = length(x);
    y = Double_Type[n];
    x0 = x - 0.5*dx;
    x1 = x + 0.5*dx;
    r = histogram(e, x0, x1, &w);
    r = @x;
    for (i = 0; i < length(w); i++) {
      r[i] = sum(s[w[i]]);
    }

    y = conv(a, x, r, sig);
  }
  return (x,y);
}

define convfn(fn, x0, x1, sig) {
  variable e, s, lam, w, i, x, y, a;

  if (_NARGS == 5) a = ();
  else a = 0.0;

  (e, s) = readcol(fn, 5, 7);
  lam = const->hc/e;
  w = where(lam >= x0 and lam <= x1);
  i = where(s[w] > 1e-4*max(s[w]));
  (x, y) = convd(a, lam[w[i]], s[w[i]], sig);
  return (x, y);
}
