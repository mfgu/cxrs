() = evalfile("rcols");


% evaluate a polynomial.
define dpoly(x, c) {
  variable n, i, y;

  n = length(c);

  if (typeof(x) == Array_Type) {
    y = Double_Type[length(x)];
    y[*] = c[-1];
  } else {
    y = c[-1];
  }
  for (i = n-2; i >= 0; i--) {
    y = y*x + c[i];
  }
  
  return y;
}

% evaluate an inverse polynomial
define ipoly(y, c) {
  variable n, i, j, x, x1, y1, c1, m;
 
  n = length(c);
  m = length(y); 
  x = Double_Type[m];
  if (n == 2) {
    x = (y - c[0])/c[1];
  } else if (n == 3) {
    if (abs(c[2]) < 1E-10) {
      x = (y - c[0])/c[1];
    } else {
      x[*] = 1E100;
      c1 = c[1]^2 - 4.0*(c[0]-y)*c[2];    
      i = where(c1 >= 0.0);
      if (length(i) > 0) {
	if (c[1] > 0) {
	  x[i] = (-c[1] + sqrt(c1[i]))/(2.0*c[2]);
	} else {
	  x[i] = (-c[1] - sqrt(c1[i]))/(2.0*c[2]);
	}
      }
    }
  } else {
    c1 = Double_Type[n-1];
    for (i = 1; i < n; i++) {
      c1[i-1] = c[i]*i;
    }
    x = [(y-c[0])/c[1]];
    x = (y-c[0]-dpoly(x,c[[2:]])*x^2)/c[1];    
    for (j = 0; j < m; j++) {
      while (1) {
	i = dpoly(x[j], c1);
	if (i*c[1] <= 0) {
	  x[j] = 1E100;
	  break;
	}
	y1 = dpoly(x[j], c);
	x1 = x[j] + (y[j] - y1)/i;
	if (abs(x1) > 1E-8) {
	  if (abs((x[j]-x1)/x1) < 1E-6) break;
	} else {
	  if (abs(x[j]-x1) < 1E-10) break;
	}
	x[j] = x1;
      }    
    }
  }
  if (typeof(y) != Array_Type) x = x[0];

  return x;
}

define newspec() {
  variable s;

  s = struct {bin_lo, bin_hi, value, err};
  
  return s;
}

% rhist reads the kmax histogram file
% sp = rhist(fn)
define rhist(fn) {
  variable sp, sz, n, x, i, f, nt, y, fmt;
  variable endian = ">";

  f = fopen(fn, "r");
  () = fseek(f, 0, SEEK_END);
  sz = ftell(f);
  () = fseek(f, 512, SEEK_SET);
  
  n = (sz - 112)/4;
  sp = struct {bin_lo, bin_hi, value, err};
  sp.bin_lo = [1.0:n+1.0:1.0];
  sp.bin_hi = sp.bin_lo + 1.0;
  sp.value = Double_Type[n];
  fmt = sprintf("%sk", endian);
  for (i = 0; i < n; i++) {
    nt = fread(&x, Char_Type, 4, f);
    if (typeof(x) != BString_Type) x = array_to_bstring(x);
    sp.value[i] = double(unpack(fmt, x));
  }

  () = fclose(f);

  sp.err = sqrt(sp.value);  

  return sp;
}

% rkmax reads the kmax event file
% evts = rkmax([endian,] fn)
% evts is a 2-dim 4-byte integer array [n,np];
% where np is the number of parameters per event,
% and n is the number of total events
% fn is the file name
% endian is optional parameter for byte-order of the file,
% it may be ">", "<", or "=" for big, little, or native endian,
% default is for big endian.
define rkmax(fn) {
  variable f, sz;
  variable h, nb, nt;
  variable i, k, x, fmt, fmt1;
  variable np, n, evts;
  variable endian = ">";

  if (_NARGS == 2) endian = ();

  f = fopen(fn, "r");
  () = fseek(f, 0, SEEK_END);
  sz = ftell(f);
  () = fseek(f, 0, SEEK_SET);
  nt = fread(&h, 362, f);
  if (typeof(h) != BString_Type) h = array_to_bstring(h);
  nt = fread(&h, 2, f);
  if (typeof(h) != BString_Type) h = array_to_bstring(h);
  fmt = sprintf("%1sJ", endian);
  np = unpack(fmt, h);
  nt = fread(&h, 370, f);
  if (typeof(h) != BString_Type) h = array_to_bstring(h);
  
  n = (sz-734)/(np*4);
  evts = Integer_Type[n,np];
  k = 0;
  fmt = sprintf("x4%sK", endian);
  fmt1 = sprintf("%sK%d", endian, np);
  while (1) {
    nt = fread(&h, 8, f);
    if (nt == -1) break;
    if (typeof(h) != BString_Type) h = array_to_bstring(h);
    nb = unpack(fmt, h);
    for (i = 0; i < nb; i++) {
      nt = fread(&h, np*4, f);
      if (nt == -1) break;
      if (typeof(h) != BString_Type) h = array_to_bstring(h);
      x = unpack(fmt1, h);
      evts[k,*] = x;
      k++;
    }
    if (nt == -1) break;
  }

  () = fclose(f);

  evts = evts[[0:k-1],*];

  % for all ADC channels, cutoff bits higher than 12. 
  evts[*,[1:]] = evts[*,[1:]]&0xfff;

  return evts;
}

% make a histogram of the paramter i of evts
% h = mkhist([w,] evts, i, xmin, xmax, n)
% h is the histogram structure returned
% evts is the evts array
% i is the parameter index.
% xmin is the minimum of the channel
% xmax is the maximum of the channel
% n is the number of bins
% w is an optional index array to filter evts.
define mkhist(evts, i, xmin, xmax, n) {
  variable w;

  if (_NARGS == 6) w = ();
  else w = [0:length(evts[*,i])-1];

  variable lo, hi, h, err, dx;
  
  dx = (xmax - xmin)/double(n);
  lo = [xmin:xmax:dx];
  hi = lo + dx;
  h = histogram(evts[w,i], lo, hi);
  
  variable x = struct {bin_lo, bin_hi, value, err};
  x.bin_lo = lo;
  x.bin_hi = hi;
  x.value = h;
  x.err = 1+sqrt(x.value + 0.75);

  return x;
}

define agrebin(lo, hi, lo0, hi0, y0, w0) {
  variable y, i, w;

  if (w0 != NULL) {
    for (i = w0[0]; i >= 0; i--) {
      if (hi0[i] <= lo0[i]) break;
    }
    w = i+1;
    for (i = w0[-1]; i < length(lo0); i++) {
      if (hi0[i] <= lo0[i]) break;
    }
    w = [w:i-1];    
    i = where(hi0[w] > lo0[w]);
    w = w[i];
  } else {
    w = where(hi0 > lo0);
  }
  y = rebin(lo, hi, lo0[w], hi0[w], y0[w]);
  return y;
}

define collapse(sp, c) {
  variable lo, hi, y, sp1;

  sp1 = struct {bin_lo, bin_hi, value, err};
  if (c <= 1) {
    sp1.bin_lo = @(sp.bin_lo);
    sp1.bin_hi = @(sp.bin_hi);
    sp1.value = @(sp.value);
    sp1.err = @(sp.err);
  } else {
    lo = [sp.bin_lo[0]:sp.bin_hi[-1]:c*(sp.bin_hi[0]-sp.bin_lo[0])];
    hi = make_hi_grid(lo);
    y = rebin(lo, hi, sp.bin_lo, sp.bin_hi, sp.value);
    sp1.bin_lo = lo;
    sp1.bin_hi = hi;
    sp1.value = y;    
    sp1.err = sqrt(y);
  }
  return sp1;
}

define reverse_spec(sp) {
  variable a, lo, hi, sp1;

  sp1 = struct {bin_lo, bin_hi, value, err};

  a = sp.bin_hi[-1]+1.0;
  lo = reverse(a-sp.bin_hi);
  hi = reverse(a-sp.bin_lo);
  sp1.bin_lo = lo;
  sp1.bin_hi = hi;
  sp1.value = reverse(sp.value);
  sp1.err = reverse(sp.err);

  return sp1;
}
  
% read the ccd file from the gold spectromenter
% img = rccd(fn, nx, ny, mode)
% img is the 2-dim img retured.
% fn is the file name
% nx is the number of rows
% ny is the number of columns.
% mode indicates if the files was from windows or mac machine. 
% 0 for windows IPLAB, 
% 1 for mac IPLAB, 
% 2 for tiff in mac, 
% 3 for tiff in windows.
% note that the file store the image in the row-major order.
%%
% in the new version, mode is omitted. the correct mode is 
% determined from the file size.
define rccd(fn, nx, ny) {
  variable endian, mode;
  variable f, n, x, fmt, img;
  variable i, nhdr, sz, isz;

  f = fopen(fn, "r");
  if (f == NULL) return NULL;
  () = fseek(f, 0, SEEK_END);
  sz = ftell(f);
  isz = 2*nx*ny;
  nhdr = sz - isz;
  
  if (nhdr == 4100) {
    mode = 0;
    endian = "<";
  } else if (nhdr == 2122) {
    mode = 1;
    endian = ">";
  } else if (nhdr == 198) {
    mode = 2;
    endian = ">";
  } else if (nhdr == 146) {
    mode = 3;
    endian = "<";
  } else {
    vmessage("Invalid image file");
    return NULL;
  }

  () = fseek(f, nhdr, SEEK_SET);
  fmt = sprintf("%sJ%d", endian, ny);
  img = Double_Type[nx,ny];
  for (i = 0; i < nx; i++) {
    n = fread(&x, ny*2, f);
    if (typeof(x) != BString_Type) x = array_to_bstring(x);
    img[i,*] = unpack(fmt, x);
  }
  
  () = fclose(f);

  return img;
}

define rccdm(fn) {  
  return rccd(fn, 1300, 1340);
}

% sum the ccd image in the cross dispersion direction
% h = ccdsp(img)
% h is the returned histogram structure.
% h.err is undefined here
% img is the input ccd image.
define ccdsp(img, y0, y1) {
  variable x, nx;

  if (y1 <= 0) y1 = length(img[*,0])-1;
  if (y0 <= 0) y0 = 0;
  x = struct {bin_lo, bin_hi, value, err};
  nx = length(img[0,*]);
  x.bin_lo = [1.0:nx+0.1:1.0];
  x.bin_hi = x.bin_lo + 1.0;
  x.value = sum(img[[y0:y1],*], 0);
  x.err = @x.value;
  return x;
}

% shift every other row of the image to get around the 
% alternating pattern in the readout noise.
define shiftccd(img) {
  variable sz, ny, img1, i;

  (sz,,) = array_info(img);
  ny = sz[0];
  img1 = @img;
  for (i = 0; i < ny; i += 2) {
    img1[i,*] = shift(img[i,*], 1);
  }
  return img1;
}

define rotatecoeff(img, n, xr) {
  variable sz, nx, ny, s1, s2, w, x, y1, y2, i;
  variable w0, nw, nb, nr, p, a, b;

  (sz,,) = array_info(img);
  nx = sz[0];
  ny = sz[1];
  s1 = ccdsp(img, 0, n-1);
  s2 = ccdsp(img, nx-n, nx-1);
  nb = length(xr);
  nr = nb/2;
  x = Double_Type[nr];
  b = @x;
  for (i = 0; i < nb; i += 2) {
    w = where(s1.bin_lo >= xr[i] and s1.bin_hi <= xr[i+1]);
    nw = length(w);
    if (nw mod 2) nw--;
    w = w[[1:]];
    w0 = where(s1.value[w] == max(s1.value[w]));
    w0 = w[w0[0]];
    x[i/2] = s1.bin_lo[w0]+0.5;
    y1 = fft(s1.value[w]-min(s1.value[w]), 0);
    y2 = fft(s2.value[w]-min(s2.value[w]), 0);
    p = abs(fft(y1*Conj(y2), 1));
    nw = length(w);
    p = shift(p, nw/2);
    p -= min(p);
    p /= max(p);
    w0 = where(p == 1.0);
    w0 = [w0[0]-3:w0[0]+3];  
    b[i/2] = sum(w0*p[w0])/sum(p[w0]) - nw/2;
    vmessage("shift: %2d %6.3f %6.2f", i, x[i/2], b[i/2]);    
  }

  (a, sz) = array_fit(x, b, NULL, [0.0, 0.0], NULL, NULL, &dpoly);
  a = (a*nx)/(nx-n);

  return a;
}
  
define rotateccd(img, d) {
  variable sz, nx, ny, img1, s, i;
  variable lo0, hi0, lo, hi, x;

  (sz,,) = array_info(img);
  nx = sz[0];
  ny = sz[1];
  img1 = @img;
  (lo0, hi0) = linear_grid(0.0, ny, ny);
  if (typeof(d) != Array_Type) d = [d];
  x = 0.5*(lo0+hi0);
  x = dpoly(x, d);
  for (i = 0; i < nx; i++) {
    s = (x*i)/(nx-1.0);    
    lo = lo0 + s;
    hi = hi0 + s;
    img1[i,*] = rebin(lo0, hi0, lo, hi, img[i,*]);
  }
  return img1;
}

define rotateccd0(img, d) {
  variable sz, nx, ny, img1, s, i;
  (sz,,) = array_info(img);
  nx = sz[0];
  ny = sz[1];
  img1 = @img;
  for (i = 0; i < nx; i++) {
    s = int((d*i)/(nx-1.0));
    img1[i,*] = shift(img[i,*], -s);
  }
  
  return img1;
}
  
% remove the cosmic rays in the grating image
% img1 = badpix(img, nc0, nb, nn, sg, p0, pmax)
% img1 is the cleaned image.
% img is the original image.
% nc0 is number of columns to be processed at a time. suggested: 10
% nb is the threshold number of pixels. suggested: 5
% nn is the threshold gap in the pixel values. suggested: 10
%% the cutoff of pixel values are based on the criteria that number of 
%% pixels within pixel values of gap nn must be greater than nb. 
% sg, p0 are the parameter controlling the spread of the affected pixels. 
%        suggested: (0.0, 0.5)
% pmax is the maximum pixel value allowed. suggested: 1e4
define badpix(img, nc0, nb, nn, sg, p0, pmax) {
  variable sz, nx, ny, nc, nt, sg2;
  variable img1, md, md0, s, k1, k2, xw, yw;
  variable r, w, h, i, i1, n0, j, p, p1;
  
  if (pmax <= 0) pmax = max(img);
  nc = nc0;
  (sz,md,h) = array_info(img);
  nx = sz[0];
  ny = sz[1];
  nt = nx*nc;
  p = Double_Type[nx, ny];
  p1 = Double_Type[nt];
  md0 = 1e30;  
  for (i = 0; i < ny; i += nc0) {
    i1 = i+nc;
    if (i1 > ny) {
      i1 = ny; 
      nc = i1-i;
      nt = nx*nc;
    }
    img1 = _reshape(img[*,[i:i1-1]], nt);
    s = array_sort(img1);
    sz = img1[s];
    w = [int(0.5*nt):nt-1];
    md = sz[w[0]];
    if (md < md0) md0 = md;
    p1[*] = 0.0;
    w = where(sz <= pmax);
    for (j = w[-1]; j > 0; j--) {
      n0 = where(img1 > sz[j]-nn and img1 <= sz[j]);
      if (length(n0) >= nb) {
	break;
      }
    }
    if (j < nt) {
      p1[s[[j:]]] = 1.0;
      p[*,[i:i1-1]] = _reshape(p1[[0:nt-1]], [nx, nc]);
    }    
  }
  img1 = @img;
  if (sg > 0) {
    sg2 = sg*sg;
    w = where(p == 1.0);
    if (length(w) > 0) {
      xw = w/ny;
      yw = w mod ny;
      for (i = 0; i < length(w); i++) {
	for (k1 = xw[i]-2; k1 <= xw[i]+2; k1++) {
	  if (k1 < 0 or k1 >= nx) continue;
	  for (k2 = yw[i]-2; k2 <= yw[i]+2; k2++) {
	    if (k2 < 0 or k2 >= ny) continue;
	    r = exp(-0.5*((k1-xw[i])^2+(k2-yw[i])^2)/sg2);
	    p[k1,k2] += r;
	  }
	}
      }
    }
  }
  w = where(p >= p0);
  if (length(w) > 0) {
    img1[w] = md0;
  }

  return img1;
}

define badpixm(img) {
  return badpix(img, 10, 10, 10, 0.0, 0.5, -1);
}

define dbgauss(x, p) {

  return p[0]*exp(-0.5*((x-p[1])/p[2])^2)+p[3];
}
  
define dgauss(x, p) {
  
  return p[0]*exp(-0.5*((x-p[1])/p[2])^2)/(sqrt(2*PI)*p[2]);
}

% fit function, constant background + multi-gaussian.
define bgauss(x, p) {
  variable y, y1, sig0, ksig, sig, x0, x1, x2, xb, b, kb, s;
  variable n, k, i, a, dx, j;
  
  n = (length(p)-2)/2;

  kb = p[-1];
  b = p[-2];
  ksig = p[-3];
  sig0 = p[-4];

  k = 10;
  if (k > 1) {
    dx = 1.0/k;
    x1 = Double_Type[length(x)*k];
    for (i = 0, j = 0; i < length(x); i++, j += k) {
      x1[[j:j+k-1]] = [x[i]:x[i]+1.0-0.1*dx:dx];
    }
    x2 = x1 + dx;
  } else {
    dx = 1.0;
    x1 = x;
    x2 = x + dx;
  }
  xb = 0.5*(x1+x2);
  
  k = 0;
  y = b + kb*(x-x[0]);
  for (i = 0; i < n; i++) {
    s = p[k];
    k++;
    x0 = p[k];
    k++;  
    sig = sig0 + ksig*(x0 - x[0]);
    a = 1.0/sqrt(2.0*PI*sig*sig);
    y1 = s*a*exp(-0.5*((xb-x0)/sig)^2)*dx;
    y1 = rebin(x, x+1, x1, x2, y1);
    y = y + y1;
  }

  return y;
}

define fitgspec(x, x1, y, x0, sig0, em, q) {
  variable err,p,s,y1, w;
  variable n, k, i, a, pmin, pmax;

  a = sqrt(2.0*PI*sig0*sig0);
  n = length(x0);
  p = Double_Type[n*2+4];
  pmin = Double_Type[n*2+4];
  pmax = Double_Type[n*2+4];
  k = 0;
  for (i = 0; i < length(x0); i++) {
    w = where(x < x0[i] and x1 >= x0[i]);
    p[k] = y[w[0]]*a;
    pmin[k] = 0.0;
    pmax[k] = p[k]*1E3;
    k++;
    p[k] = x0[i];
    pmin[k] = p[k]-10*sig0;
    pmax[k] = p[k]+10*sig0;
    k++;
  }
  p[k] = sig0;
  pmin[k] = 0.0;
  pmax[k] = 10*sig0;
  k++;
  p[k] = 0.0;
  pmin[k] = -1.0;
  pmax[k] = 1.0;
  k++;
  p[k] = mean(y[[0:10]]);
  pmin[k] = 1E-10;
  pmax[k] = 10.0*p[k];
  k++;
  p[k] = 0.0;
  pmin[k] = -0.1*p[k-1]/(x[-1]-x[0]);
  pmax[k] = 0.1*p[k-1]/(x[-1]-x[0]);
  
  err = sqrt(y);
  (p,s) = array_fit(x,y,1.0/err^2,p,pmin,pmax,&bgauss);
  y1 = bgauss(x,p);

  err = em*sqrt(y1*(1.0+(1.0-(p[-2]+p[-1]*(x-x[0]))/y1)*q));
  (p,s) = array_fit(x,y,1.0/err^2,p,pmin,pmax,&bgauss);  
  y1 = bgauss(x,p);

  return p,pmin,pmax,s,y1;
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
      H[w] = exp(-v2)*(1. + fac2*alpha^2 * (1. - 2.*v2)) + fac1 * (alpha/v2);
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

define dvoigt(x, p) {
  variable s2, t, y;

  s2 = sqrt(2.0)*p[2];
  t = (x-p[1])/s2;
  y = (1.0/s2)*p[0]*voigt(p[3], t);
  return y;
}

define dlorentz(x, p) {
  
  return (p[2]/(2.0*PI))/((x-p[1])^2 + 0.25*p[1]^2);
}

% fit function with constant background + multi-voigts.
define bvoigt(x, p) {
  variable y, y1, sig, sig2, a, x0, x1, x2, xb, b, kb, s;
  variable n, k, i, u, dx, j;
  variable ksig, sig0;

  n = (length(p)-3)/2;
  
  b = p[-1];
  a = p[-2];
  ksig = p[-3];
  sig0 = p[-4];

  k = 10;
  if (k > 1) {
    dx = 1.0/k;
    x1 = Double_Type[length(x)*k];
    for (i = 0, j = 0; i < length(x); i++, j += k) {
      x1[[j:j+k-1]] = [x[i]:x[i]+1.0-0.1*dx:dx];
    }
    x2 = x1 + dx;
  } else {
    dx = 1.0;
    x1 = x;
    x2 = x + dx;
  }
  xb = 0.5*(x1+x2);

  k = 0;
  y = b;
  for (i = 0; i < n; i++) {
    s = p[k];
    k++;
    x0 = p[k];
    k++;
    sig = sig0 + ksig*(x0 - x[0]);
    sig2 = sqrt(2.0)*sig;
    u = (xb-x0)/sig2;
    y1 = s*voigt(a, u)*dx/sig2;
    y1 = rebin(x, x+1.0, x1, x2, y1);
    y = y + y1;
  }

  return y;
}
  
define fitcspec(x, x1, y, x0, sig0, a0) {
  variable err, p, s, y1, w;
  variable n, k, i, a, pmin, pmax;

  a = sqrt(2.0*PI*sig0*sig0);
  n = length(x0);
  p = Double_Type[n*2+4];
  pmin = Double_Type[n*2+4];
  pmax = Double_Type[n*2+4];
  k = 0;
  for (i = 0; i < length(x0); i++) {
    w = where(x < x0[i] and x1 >= x0[i]);
    p[k] = y[w[0]]*a;
    pmin[k] = 0.0;
    pmax[k] = 1E3*p[k];
    k++;
    p[k] = x0[i];
    pmin[k] = p[k]-10*sig0;
    pmax[k] = p[k]+10*sig0;
    k++;
  }
  p[k] = sig0;
  pmin[k] = 0.0;
  pmax[k] = sig0*100.0;
  k++;
  p[k] = 0.0;
  pmin[k] = -1.0;
  pmax[k] = 1.0;
  k++;
  p[k] = a0;
  pmin[k] = 0.0;
  pmax[k] = 10.0;
  k++;
  p[k] = mean(y[[0:10]]);
  pmin[k] = 1E-10;
  pmax[k] = 10.0*p[k];
  
  err = 1 + sqrt(y+0.75);
  (p,s) = array_fit(x,y,1.0/err^2,p,pmin,pmax,&bvoigt);
  y1 = bvoigt(x,p);
  err = 1+sqrt(y1+0.75);
  (p,s) = array_fit(x,y,1.0/err^2,p,pmin,pmax,&bvoigt);
  y1 = bvoigt(x,p);

  return p,pmin,pmax,s,y1;
}

% read crystal response data from a file
define rxtal(fn) {
  variable f, d, s, n, i;
  variable x1,x2,x3,x4,x5,x6;

  d = struct {d2, e, rm, rp, rprs, wp, ws};

  f = fopen(fn, "r");
  s = fgetslines(f);
  () = sscanf(s[0], "2d=%f\n", &x1);
  d.d2 = x1;
  n = length(s)-2;
  d.e = Double_Type[n];
  d.rm = @(d.e);
  d.rp = @(d.e);
  d.rprs = @(d.e);
  d.wp = @(d.e);
  d.ws = @(d.e);
  for (i = 0; i < n; i++) {
    () = sscanf(s[i+2], "%f %f %f %f %f %f\n", 
		&x1, &x2, &x3, &x4, &x5, &x6);
    d.e[i] = x1;
    d.rm[i] = x2;
    d.rp[i] = x3;
    d.rprs[i] = x4;
    d.wp[i] = x5;
    d.ws[i] = x6;
  }
 
  d.e = 1E3*(d.e);
  () = fclose(f);

  return d;
}

% read xrs event file in ascii format. 
% evts = rxrs(fn)
% the file must contain 5 params per event:
% Ph, Volts, Flags, UTime, Pixel
define rxrs(fn) {
  variable evts, s, f, n, i, e;
  variable a0, a1, a2, a3, a4, t0;
  
  f = fopen(fn, "r");
  s = fgetslines(f);
  
  n = length(s);
  evts = Double_Type[n, 5];
  for (i = 0; i < 10; i++) {
    e = sscanf(s[i], "%lf%lf%lf%lf%lf", &a0, &a1, &a2, &a3, &a4);
    if (e != 5) break;
    evts[i,0] = a0;
    evts[i,1] = a1;
    evts[i,2] = a2;
    evts[i,3] = a3;
    evts[i,4] = a4;
  }

  t0 = evts[0,3];
%  evts[*,3] = evts[*,3] - t0;
  return evts[[0:i-1],*];
}

% make a xrs spectrum with either Ph or Volts
% sp = xrsp(evts, k, xmin, xmax, npts, gts, det)
% k = 0 for Ph spectrum, 1 for Volts spectrum.
% xmin, xmax, npts give the lower, upper limits and num. points
% gts is an array of even number of elements giving the good time intervals.
% in the form [start1, stop1, start2, stop2, ...]
% det = 0 include all pixels except 18, 28.
% if det > 0, then only pixel=det is include.
define xrssp(evts, k, xmin, xmax, npts, gts, det) {
  variable w, w1, i, sp;

  if (det > 0) {
    w = (evts[*,4] == det);
  } else {
    w = evts[*,4] != 18 and evts[*,4] != 28;
  }
  if (length(gts) > 0){
    w1 = 0;
    for (i = 0; i < length(gts); i += 2) {
      w1 = w1 or (evts[*,3] >= gts[i] and evts[*,3] < gts[i+1]);
    }
  } else {
    w1 = 1;
  }
  w = where(w and w1);

  sp = mkhist(w, evts, k, xmin, xmax, npts);
  
  return sp;
}

% write a histogram to a file in ascii format with writecol.
define whista(fn, sp) {
  if (sp.err != NULL) {
    writecol(fn, sp.bin_lo, sp.bin_hi, sp.value, sp.err);
  } else {
    sp.err = sqrt(sp.value);
    writecol(fn, sp.bin_lo, sp.bin_hi, sp.value);
  }
}

define xy2hist(x, y) {
  variable sp, i;
  
  sp = struct {bin_lo, bin_hi, value, err};
  i = array_sort(x);
  
  sp.bin_lo = x[i] - 0.5*(x[i[1]]-x[i[0]]);
  sp.bin_hi = make_hi_grid(sp.bin_lo);
  sp.value = y[i];
  sp.err = sqrt(abs(y[i]));
  return sp;
}

% read a ascii spectra from file written with writecol.
define rhista(fn) {
  variable sp;

  sp = struct {bin_lo, bin_hi, value, err};
  
  (sp.bin_lo, sp.bin_hi, sp.value) = readcol(fn, 1, 2, 3);
  sp.err = readcol(fn, 4);

  return sp;
}

define rwava(fn) {
  variable sp, f, buf, n, x, i, k;

  sp = struct {bin_lo, bin_hi, value, err};
  
  f = fopen(fn, "r");
  buf = fgetslines(f);
  () = fclose(f);
  if (length(buf) > 1) {
    if (_slang_guess_type(strtrim(buf[0][[:-2]])) == String_Type) buf = buf[[1:]];
    n = length(buf);
    sp.bin_lo = [1.0:n+0.1:1.0];
    sp.bin_hi = sp.bin_lo + 1.0;
    sp.value = array_map(Double_Type, &atof, buf);
    sp.err = 1.0 + sqrt(sp.value + 0.75);
  } else {
    buf = strchop(strcompress(buf[0], " \t"), '\r', 0);
    buf = array_map(String_Type, &strtrim, buf);
    k = array_map(Integer_Type, &isdigit, buf);
    k = where(k);
    n = length(k);
    sp.bin_lo = [1.0:n+0.1];
    sp.bin_hi = sp.bin_lo + 1.0;
    sp.value = Double_Type[n];
    sp.err = Double_Type[n];
    for (i = 0; i < n; i++) {
      sp.value[i] = atof(buf[k[i]]);
      sp.err[i] = 1.0+sqrt(sp.value[i]+0.75);
    }
  }
  return sp;
}
  
% read a histogram in Mac line ending
define machist(fn) {
  variable sp, f, buf, a, n, i, x, y, w;
  
  sp = struct {bin_lo, bin_hi, value, err};
  f = fopen(fn, "r");
  buf = fgetslines(f);
  if (length(buf) == 1) {    
    a = strchop(strcompress(buf[0][[:-2]], " \t"), '\r', 0);
  }
  if (_slang_guess_type(strcompress(a[0], " \t")) == String_Type) a = a[[1:]];
  n = length(a);
  sp.bin_lo = Double_Type[n];
  sp.bin_hi = Double_Type[n];
  sp.value = Double_Type[n];
  sp.err = Double_Type[n];
  for (i = 0; i < n; i++) {
    () = sscanf(a[i], "%f %f", &x, &y);
    sp.bin_lo[i] = x;
    sp.value[i] = y;
  }
  w = where(sp.bin_lo > 0);
  sp.bin_lo = sp.bin_lo[w];
  sp.value = sp.value[w];
  sp.bin_hi = make_hi_grid(sp.bin_lo);
  sp.err = sqrt(sp.value);
  return sp;
}

define igauss(lo, hi, x, sig) {
  variable a, xt, dx, w, y, i, sig0, s, k, blo, bhi, alo, ahi, x0, x1;
  
  sig0 = 0.2*sig;
  a = 1.0/sqrt(2.0*PI);
  y = Double_Type[length(lo)];

  dx = hi - lo;
  s = where(dx <= sig0);
  if (length(s) > 0) {
    xt = 0.5*(lo[s] + hi[s]);
    xt = (x - xt)/sig;
    w = where(xt > -10.0 and xt < 10.0);
    y[s[w]] = a*exp(-0.5*xt[w]*xt[w])*dx[s[w]]/sig;
  }  
  s = where(dx > sig0);
  if (length(s) > 0) {
    i = where(lo[s] <= x and hi[s] > x);
    if (length(i) == 1) {
      i = i[0];
      for (k = i-1; k >= 0; k--) {
	x0 = (hi[s[k]]-x)/sig;
	if (x0 < -10.0) break;
      }      
      w = k+1;
      for (k = i+1; k < length(s); k++) {
	x0 = (lo[s[k]]-x)/sig;
	if (x0 > 10.0) break;
      }
      w = [w:k-1];
      for (i = 0; i < length(w); i++) {
	k = s[w[i]];
	x0 = (lo[k]-x)/sig;
	x1 = (hi[k]-x)/sig;
	if (x0 < -5 and x1 > 5) y[k] = 1.0;
	else if (x0 >= -5) {
	  x1 = min([x0 + 20.0, x1]);	
	  (blo, bhi) = linear_grid(x0, x1, 100);
	  xt = 0.5*(blo+bhi);
	  y[k] = sum(a*exp(-0.5*xt*xt)*(bhi-blo));	
	} else if (x1 <= 5) {
	  x0 = max([x1 - 20.0, x0]);
	  (blo, bhi) = linear_grid(x0, x1, 100);
	  xt = 0.5*(blo+bhi);
	  y[k] = sum(a*exp(-0.5*xt*xt)*(bhi-blo));
	}
      }      
    }
  }

  return y;
}

define ivoigt(lo, hi, x, sig, alpha) {
  variable a, xt, dx, w, y;

  a = 1.0/(1.41421356*sig);
  dx = hi - lo;
  xt = 0.5*(lo + hi);
  xt = (xt - x)*a;
  w = where(xt > -50.0 and xt < 50.0);
  y = Double_Type[length(lo)];
  y[w] = a*voigt(alpha, xt[w])*dx[w];
  
  return y;
}

% general fit function.
define mlines_fit(lo, hi, par) {
  variable nc, nw, na, nb, ni, nlines, mode;
  variable abund, pos0, pos, str;
  variable a, c, w, b, i, j, k;
  variable x, sig, h, yt, y;

  i = 0;

  mode = int(par[i]);
  i++;

  nc = int(par[i]);
  i++;

  nw = int(par[i]);
  i++;

  na = int(par[i]);
  i++;

  nb = int(par[i]);
  i++;

  ni = int(par[i]);
  i++;

  nlines = int(par[[i:i+ni-1]]);
  i += ni;

  pos0 = par[i];
  i++;

  if (nc > 0) {
    c = par[[i:i+nc-1]];
  } else {
    c = [pos0, 1.0];
  }
  i += nc;
  
  w = par[[i:i+nw-1]];
  i += nw;

  if (mode < 0) {
    a = par[[i:i+na-1]];
    i += na;
  }

  b = par[[i:i+nb-1]];
  i += nb;

  abund = par[[i:i+ni-1]];
  i += ni;

  y = dpoly(0.5*(lo+hi)-pos0, b);
  for (k = 0; k < ni; k++) {
    for (j = 0; j < nlines[k]; j++) {
      str = abund[k]*par[i];
      i++;
      pos = par[i]-pos0;
      i++;
      sig = dpoly(pos, w);
      if (mode < 0) {
	h = dpoly(pos, a);
      }
      x = dpoly(pos, c);
      if (mode > 0) {
	yt = str*igauss(lo, hi, x, sig);
      } else {
	yt = str*ivoigt(lo, hi, x, sig, h);
      }
      y += yt;
    }
  }

  return y;
}

% add the mlines fit function into the internal table
% and set the frozen parameters.
define add_mlines(mode, nc, nw, na, nb, nlines, p0) {
  variable s, par, i, np, j, k, ni;

  if (mode > 0 and na > 0) {
    vmessage("incompartible parameters\n");
    return;
  }

  ni = length(nlines);
  np = 7 + nc + nw + na + nb + 2*ni + 2*int(sum(nlines));

  par = String_Type[np];

  i = 0;
  par[i] = "mode";
  i++;
  par[i] = "ncalib";
  i++;
  par[i] = "nwidth";
  i++;
  par[i] = "nalpha";
  i++;
  par[i] = "nbkgd";
  i++;
  par[i] = "nion";
  i++;
  for (k = 0; k < ni; k++) {
    s = sprintf("nlines%d", k);
    par[i] = s;
    i++;
  }
  par[i] = "pcenter";
  i++;
  for (k = 0; k < nc; k++) {
    s = sprintf("calib%d", k);
    par[i] = s;
    i++;
  }
  for (k = 0; k < nw; k++) {
    s = sprintf("width%d", k);
    par[i] = s;
    i++;
  }
  for (k = 0; k < na; k++) {
    s = sprintf("alpha%d", k);
    par[i] = s;
    i++;
  }
  for (k = 0; k < nb; k++) {
    s = sprintf("bkgd%d", k);
    par[i] = s;
    i++;
  }
  
  for (k = 0; k < ni; k++) {
    s = sprintf("a%d", k);
    par[i] = s;
    i++;
    for (j = 0; j < nlines[k]; j++) {
      s = sprintf("s%d_%d", k, j);
      par[i] = s;
      i++;
      s = sprintf("p%d_%d", k, j);
      par[i] = s;
      i++;
    }
  }

  del_function("mlines");
  add_slang_function("mlines", par);
  
  fit_fun("mlines(1)");

  i = 1;
  set_par(i, mode, 1);
  i++;
  set_par(i, nc, 1);
  i++;
  set_par(i, nw, 1);
  i++;
  set_par(i, na, 1);
  i++;
  set_par(i, nb, 1);
  i++;
  set_par(i, ni, 1);
  i++;
  for (k = 0; k < ni; k++) {
    set_par(i, nlines[k], 1);
    i++;
  }
  set_par(i, p0, 1);
  i++;
  
  return;
}

#ifnexists MSpec_Type
typedef struct {
  fxtal,
  feff,
  d2, %crystal 2d spacing.
  xi, %xi[0][*]=wavlength grid, %xi[1][*]=(1-R_pi/R_sigma)/(1+R_pi/R_sigma)
  eff, % spectrometer efficiency. do not include xtal reflectivity.
  chan0,
  calib, %wavelength scale conversion from channel numbers.
  xwidth, % width x-points.
  width, %line width in gaussian sigma
  iwidth,
  alpha, %voigt damping parameter.
  ialpha,
  scale,
  iscale,
  xbkgd, %background x-points
  bkgd, %background
  ibkgd,
  abund, %abundance of ions
  iabund,
  shifts,
  ishifts,
  ions,
  lines, %lines[0]=lam,lines[1]=str,lines[2]=ch,lines[3]=P,lines[4]=anisotropy
  ilines,
  pid
} MSpec_Type;
#endif

if (is_defined("_spec") == 0) {
  variable _spec = @MSpec_Type;
}

% read the response file. 
define read_response(fn) {
  variable f, d, a, buf, i, x, s;

  f = fopen(fn, "r");
  buf = fgetslines(f);

  () = fclose(f);
  _spec.d2 = 0.0;
  _spec.calib = NULL;
  _spec.xwidth = NULL;
  _spec.xbkgd = NULL;
  _spec.alpha = Double_Type[0];
  _spec.ialpha = Integer_Type[0];
  _spec.scale = Double_Type[0];
  _spec.iscale = Integer_Type[0];
  _spec.eff = NULL;
  _spec.fxtal = NULL;
  _spec.feff = NULL;
  for (i = 0; i < length(buf); i++) {
    a = strcompress(buf[i][[:-2]], " \t");
    a = strchop(a, ' ', 0);
    if (length(a) == 0) continue;
    switch (a[0])
    {case "chan0": _spec.chan0 = atof(a[1]);}
    {case "calib": 
     _spec.calib = array_map(Double_Type, &atof, a[[1:]]);
     _spec.calib = _reshape(_spec.calib, [2,length(_spec.calib)/2]);
    }
    {case "xwidth": _spec.xwidth = array_map(Double_Type, &atof, a[[1:]]);}
    {case "width": _spec.width = array_map(Double_Type, &atof, a[[1:]]);}
    {case "alpha": _spec.alpha = array_map(Double_Type, &atof, a[[1:]]);}
    {case "xbkgd": _spec.xbkgd = array_map(Double_Type, &atof, a[[1:]]);}
    {case "bkgd": _spec.bkgd = array_map(Double_Type, &atof, a[[1:]]);}
    {case "scale": _spec.scale = array_map(Double_Type, &atof, a[[1:]]);}
    {case "iwidth": _spec.iwidth = array_map(Integer_Type, &integer, a[[1:]]);}
    {case "ialpha": _spec.ialpha = array_map(Integer_Type, &integer, a[[1:]]);}
    {case "ibkgd": _spec.ibkgd = array_map(Integer_Type, &integer, a[[1:]]);}
    {case "iscale": _spec.iscale = array_map(Integer_Type, &integer, a[[1:]]);}
    {case "xtal": 
     _spec.fxtal = a[1];
     d = rxtal(a[1]); 
     _spec.d2 = d.d2;
     _spec.xi = Double_Type[3,length(d.e)];
     _spec.xi[0,*] = reverse(const->hc/d.e);
     _spec.xi[1,*] = reverse((1-d.rprs)/(1+d.rprs));
     _spec.xi[2,*] = reverse(0.5*(d.rp+d.rm));
    }
    {case "eff": 
     _spec.feff = a[1];
     (x,d) = readcol(a[1], 1, 2);
     s = array_sort(x);
     x = x[s];
     d = d[s];
     _spec.eff = Double_Type[2,length(x)];
     _spec.eff[0,*] = x;
     _spec.eff[1,*] = d;     
    }
  }
  if (_spec.d2 > 0 and _spec.eff != NULL) {
    x = interpol_points(_spec.eff[0,*], _spec.xi[0,*], _spec.xi[2,*]);
    _spec.eff[1,*] *= x;
  }
}

define write_response(fn) {
  variable f, i, j;

  f = fopen(fn, "w");
  () = fprintf(f, "chan0\t%10.3E\n", _spec.chan0);
  if (_spec.fxtal != NULL) {
    () = fprintf(f, "xtal\t%s\n", _spec.fxtal);
  }
  if (_spec.feff != NULL) {
    () = fprintf(f, "eff\t%s\n", _spec.feff);
  }
  if (_spec.calib != NULL) {
    () = fprintf(f, "calib");
    for (j = 0; j < 2; j++) {
      for (i = 0; i < length(_spec.calib[j,*]); i++) {
	() = fprintf(f, "\t%12.5E", _spec.calib[j, i]);
      }
    }
    () = fprintf(f, "\n");
  }
  if (length(_spec.width) > 0) {
    if (_spec.xwidth != NULL) {
      () = fprintf(f, "xwidth");
      for (i = 0; i < length(_spec.xwidth); i++) {
	() = fprintf(f, "\t%12.5E", _spec.xwidth[i]);
      }      
      () = fprintf(f, "\n");
    }
    () = fprintf(f, "width");
    for (i = 0; i < length(_spec.width); i++) {
      () = fprintf(f, "\t%12.5E", _spec.width[i]);
    }
    () = fprintf(f, "\n");
    () = fprintf(f, "iwidth");
    for (i = 0; i < length(_spec.iwidth); i++) {
      () = fprintf(f, "\t%d", _spec.iwidth[i]);
    }
    () = fprintf(f, "\n");
  }
  if (length(_spec.alpha) > 0) {
    () = fprintf(f, "alpha");
    for (i = 0; i < length(_spec.alpha); i++) {
      () = fprintf(f, "\t%12.5E", _spec.alpha[i]);
    }
    () = fprintf(f, "\n");
    () = fprintf(f, "ialpha");
    for (i = 0; i < length(_spec.ialpha); i++) {
      () = fprintf(f, "\t%d", _spec.ialpha[i]);
    }
    () = fprintf(f, "\n");
  }
  if (length(_spec.scale) > 0) {
    () = fprintf(f, "scale");
    for (i = 0; i < length(_spec.scale); i++) {
      () = fprintf(f, "\t%12.5E", _spec.scale[i]);
    }
    () = fprintf(f, "\n");
    () = fprintf(f, "iscale");
    for (i = 0; i < length(_spec.iscale); i++) {
      () = fprintf(f, "\t%d", _spec.iscale[i]);
    }
    () = fprintf(f, "\n");
  }
  if (length(_spec.bkgd) > 0) {
    if (_spec.xbkgd != NULL) {
      () = fprintf(f, "xbkgd");
      for (i = 0; i < length(_spec.xbkgd); i++) {
	() = fprintf(f, "\t%12.5E", _spec.xbkgd[i]);
      }
      () = fprintf(f, "\n");
    }
    () = fprintf(f, "bkgd");
    for (i = 0; i < length(_spec.bkgd); i++) {
      () = fprintf(f, "\t%12.5E", _spec.bkgd[i]);
    }
    () = fprintf(f, "\n");
    () = fprintf(f, "ibkgd");
    for (i = 0; i < length(_spec.ibkgd); i++) {
      () = fprintf(f, "\t%d", _spec.ibkgd[i]);
    }
    () = fprintf(f, "\n");
  }
  () = fclose(f);
}

define chan2wave();
define chan2wave(x) {
  variable sp, w, t;

  if (typeof(x) == Struct_Type) {
    sp = struct {bin_lo, bin_hi, value, err};
    w = chan2wave(x.bin_lo);
    if (w[1] > w[0]) {
      sp.bin_lo = @w;
      sp.bin_hi = chan2wave(x.bin_hi);
      sp.value = @x.value;
      if (x.err != NULL) {
	sp.err = @x.err;
      }
    } else {
      sp.bin_hi = reverse(w);
      sp.bin_lo = reverse(chan2wave(x.bin_hi));
      sp.value = reverse(x.value);
      if (x.err != NULL) {
	sp.err = reverse(x.err);
      }
    }
    return sp;
  } else {
    if (_spec.calib == NULL) return NULL;
    if (_spec.d2 > 0) {
      t = asin(_spec.calib[1,*]/_spec.d2);
    } else {
      t = _spec.calib[1,*];
    }
    w = interpol_points(x, _spec.calib[0,*], t);
    if (_spec.d2 > 0) {
      w = _spec.d2*sin(w);
    }
    if (typeof(x) != Array_Type) w = w[0];
    return w;
  }
}

define wave2chan();
define wave2chan(x) {
  variable c, sp, t;

  if (typeof(x) == Struct_Type) {
    sp = struct {bin_lo, bin_hi, value, err};
    c = wave2chan(x.bin_lo);
    if (c[1] > c[0]) {
      sp.bin_lo = @c;
      sp.bin_hi = wave2chan(x.bin_hi);
      sp.value = @x.value;
      sp.err = @x.err;
    } else {
      sp.bin_hi = reverse(c);
      sp.bin_lo = reverse(wave2chan(x.bin_hi));
      sp.value = reverse(x.value);
      sp.err = reverse(x.err);
    }
    return sp;
  } else {       
    if (_spec.calib == NULL) return NULL;
    if (_spec.d2 > 0) {
      c = asin(x/_spec.d2);
      t = asin(_spec.calib[1,*]/_spec.d2);
    } else {
      c = @x;
      t = _spec.calib[1,*];
    }
    c = interpol_points(c, t, _spec.calib[0,*]);
    if (typeof(x) != Array_Type) c = c[0];
    return c;
  }
}

% read the line model from file.
% the file must have at least 8 columns. 
% id ion low up channel strength wavelength anisotropy P [fix_c fix_s]
define read_lines(fn, conv, wlam) {
  variable f, a, buf, i, k, na, nb, nc, cmin, cmax, lines, ion, ilines, w, x;

  f = fopen(fn, "r");
  buf = fgetslines(f);
  () = fclose(f);
  nb = length(buf);
  lines = Double_Type[5,nb];
  ion = Integer_Type[nb];
  ilines = Integer_Type[7,nb];
  for (i = 0; i < length(buf); i++) {
    a = strcompress(buf[i][[:-2]], " \t");
    a = strchop(a, ' ', 0);
    na = length(a);
    if (na <= 1) {
      ion[i] = -1;
      continue;
    }
    if (not isdigit(a[0])) {
      ion[i] = -1;
      continue;
    }
    if (na >= 6) {
      ilines[4,i] = integer(a[0]);
      ion[i] = integer(a[1]);
      ilines[2,i] = integer(a[2]);
      ilines[3,i] = integer(a[3]);
      lines[0,i] = atof(a[4]);
      lines[1,i] = atof(a[5]);
      lines[2,i] = atof(a[6]);
      lines[3,i] = atof(a[7]);
      lines[4,i] = atof(a[8]);
      if (na > 9) {
	ilines[0,i] = integer(a[9]);
	ilines[1,i] = integer(a[10]);
      } else {
	ilines[0,i] = 1;
	ilines[1,i] = 1;
      }
    }
  }
  if (length(wlam) == 2) {
    w = where(lines[2,*] < wlam[0] or lines[2,*] > wlam[1]);
    if (length(w) > 0) ion[w] = -1;
  }
  if (conv) lines[2,*] = const->hc/lines[2,*];
  if (_spec.calib != NULL) {
    x = wave2chan(lines[2,*]);
    lines[0,*] = x;
  }

  w = where(ion > 0);
  cmin = min(ion[w]);
  cmax = max(ion[w]);
  nc = 0;
  for (i = cmin; i <= cmax; i++) {
    w = where(ion == i);
    if (length(w) > 0) nc++;
  }
  _spec.pid = NULL;
  _spec.lines = Array_Type[nc];
  _spec.ilines = Array_Type[nc];
  _spec.abund = Double_Type[nc];
  _spec.iabund = Integer_Type[nc];
  _spec.ions = Integer_Type[nc];
  _spec.abund[*] = 1.0;
  _spec.iabund[*] = 1;
  _spec.shifts = Double_Type[2,nc];
  _spec.ishifts = Integer_Type[nc];
  _spec.ishifts[*] = 1;
  k = 0;
  for (i = cmin; i <= cmax; i++) {
    w = where(ion == i);
    if (length(w) > 0) {
      _spec.lines[k] = lines[*,w];
      _spec.ilines[k] = ilines[*,w];
      _spec.ions[k] = i;
      k++;
    }
  }
  for (i = 0; i < length(buf); i++) {
    a = strcompress(buf[i][[:-2]], " \t");
    a = strchop(a, ' ', 0);
    na = length(a);
    if (na <= 1) {
      continue;
    }
    if (a[0] == "Abund") {
      k = where(_spec.ions == integer(a[1]));
      k = k[0];
      _spec.abund[k] = atof(a[2]);
    } else if (a[0] == "Shift") {
      k = where(_spec.ions == integer(a[1]));
      k = k[0];
      _spec.shifts[0,k] = atof(a[2]);
      _spec.shifts[1,k] = atof(a[3]);
      if (length(a) > 4) {
	_spec.ishifts[k] = integer(a[4]);
      }
    }
  }
  if (_spec.calib != NULL) {
    a = chan2wave(_spec.chan0);
    _spec.shifts[0,*] = wave2chan(a+_spec.shifts[1,*]) - _spec.chan0;
  }
}

define write_lines(fn, m) {
  variable f, k, i, lam, xi, eff;

  f = fopen(fn, "w");
  
  if (m < 2) {
    for (k = 0; k < length(_spec.abund); k++) {
      () = fprintf(f, "Abund %2d %11.4E\n", _spec.ions[k], _spec.abund[k]);
    }
    for (k = 0; k < length(_spec.abund); k++) {
      if (_spec.calib != NULL) {
	lam = chan2wave(_spec.chan0+_spec.shifts[0,k])-chan2wave(_spec.chan0);    
      } else {
	lam = 0.0;
      }
      () = fprintf(f, "Shift %2d %7.2f %11.4E %d\n", _spec.ions[k], 
		   _spec.shifts[0,k], lam, _spec.ishifts[k]);
    }
    () = fprintf(f, "\n");
  }
  for (k = 0; k < length(_spec.lines); k++) {
    if (m < 0) {
      lam = _spec.lines[k][0,*]+_spec.shifts[0,k];
    } else {
      if (_spec.calib != NULL) {
	lam = chan2wave(_spec.lines[k][0,*]+_spec.shifts[0,k]);
      } else {
	lam = _spec.lines[k][2,*];
      }
    }
    for (i = 0; i < length(_spec.lines[k][0,*]); i++) {
      if (m > 0 and _spec.ilines[k][0,i] and _spec.ilines[k][1,i]) {
	continue;
      }
      if (m < 2) {
	() = fprintf(f, 
		     "%5d %2d %4d %4d %13.7E %10.3E %13.7E %10.3E %10.3E %d %d\n",
		     _spec.ilines[k][4,i], _spec.ions[k], 
		     _spec.ilines[k][2,i], _spec.ilines[k][3,i],
		     _spec.lines[k][0,i]+_spec.shifts[0,k], 
		     _spec.lines[k][1,i], 
		     lam[i], _spec.lines[k][3,i], _spec.lines[k][4,i],
		     _spec.ilines[k][0,i], _spec.ilines[k][1,i]);
      } else {
	if (_spec.d2 > 0) {
	  xi = interpol_points(lam[i], _spec.xi[0,*], _spec.xi[1,*]);
	  xi = xi[0];
	} else {
	  xi = 0.0;
	}
	if (_spec.eff != NULL) {
	  eff = interpol_points(lam[i], _spec.eff[0,*], _spec.eff[1,*]);
	  eff = eff[0];
	} else {
	  eff = 1.0;
	}
	() = fprintf(f, 
		     "%5d %2d %4d %4d %13.7E %10.3E %13.7E %10.3E %10.3E %10.3E %10.3E\n",
		     _spec.ilines[k][4,i], _spec.ions[k], 
		     _spec.ilines[k][2,i], _spec.ilines[k][3,i],
		     _spec.lines[k][0,i]+_spec.shifts[0,k],
		     _spec.lines[k][1,i]*_spec.abund[k],
		     lam[i], _spec.lines[k][3,i], _spec.lines[k][4,i],
		     xi, eff);
      }
    }
  }
  () = fclose(f);
}

% mspec fit function
define mspec_fit(lo, hi, par) {
  variable i, j, m, k, w, mode, wi, sig;
  variable str, pos, alpha, y, yt, xi, x0;

  k = 0;
  for (i = 0; i < length(_spec.width); i++) {
    if (_spec.iwidth[i] == 0) {
      _spec.width[i] = par[k];
      k++;
    }
  }  
  for (i = 0; i < length(_spec.alpha); i++) {
    if (_spec.ialpha[i] == 0) {
      _spec.alpha[i] = par[k];
      k++;
    }
  }
  for (i = 0; i < length(_spec.scale); i++) {
    if (_spec.iscale[i] == 0) {
      _spec.scale[i] = par[k];
      k++;
    }
  }
  for (i = 0; i < length(_spec.bkgd); i++) {
    if (_spec.ibkgd[i] == 0) {
      _spec.bkgd[i] = par[k];
      k++;
    }
  }
  for (i = 0; i < length(_spec.abund); i++) {
    if (_spec.iabund[i] == 0) {
      _spec.abund[i] = par[k];
      k++;
    }
  }
  for (i = 0; i < length(_spec.ishifts); i++) {
    if (_spec.ishifts[i] == 0) {
      _spec.shifts[0,i] = par[k];
      k++;
    }
  }
  for (i = 0; i < length(_spec.lines); i++) {
    for (j = 0; j < 2; j++) {
      for (m = 0; m < length(_spec.lines[i][j,*]); m++) {
	if (_spec.ilines[i][j,m] <= 0) {
	  _spec.lines[i][j,m] = par[k];
	  k++;
	}
      }
    }
  }
  
  x0 = 0.5*(lo+hi);
  y = 0.0;
  for (m = 0; m < length(_spec.abund); m++) {
    if (_spec.abund[m] <= 0.0) continue;
    wi = where(_spec.ilines[m][0,*] < 2 and _spec.ilines[m][1,*] < 2);
    if (length(wi) == 0) continue;
    str = _spec.abund[m]*_spec.lines[m][1,wi]*_spec.lines[m][3,wi];
    pos = _spec.lines[m][0,wi]+_spec.shifts[0,m];
    if (length(_spec.scale) > 0) {
      pos = dpoly(pos, _spec.scale);
    }
    if (_spec.d2 > 0) {
      xi = interpol_points(_spec.lines[m][2,wi], _spec.xi[0,wi], _spec.xi[1,wi]);
      str *= (1.0 + xi*_spec.lines[m][4,wi]);
    }
    if (_spec.eff != NULL) {
      xi = interpol_points(_spec.lines[m][2,wi], _spec.eff[0,*], _spec.eff[1,*]);
      str *= xi;
    }
    if (_spec.xwidth == NULL or length(_spec.xwidth) == 1) {      
      sig = dpoly(pos-_spec.chan0, _spec.width);    
    } else {
      sig = interpol_points(pos, _spec.xwidth, _spec.width);
    }
    if (length(_spec.alpha) > 0) {
      if (_spec.xwidth == NULL or length(_spec.xwidth) == 1) {      
	alpha = dpoly(pos-_spec.chan0, _spec.alpha);
      } else {
	alpha = interpol_points(pos, _spec.xwidth, _spec.alpha);
      }
      mode = 0;
    } else {
      mode = 1;
    }
    for (i = 0; i < length(wi); i++) {
      if (sig[i] <= 0) sig[i] = 1e-10;
      if (mode) {
	yt = str[i]*igauss(lo, hi, pos[i], sig[i]);
      } else {
	if (alpha[i] <= 0) alpha[i] = 1e-10;
	yt = str[i]*ivoigt(lo, hi, pos[i], sig[i], alpha[i]);
      }
      y += yt;
    }
  }

  if (_spec.xbkgd == NULL or length(_spec.xbkgd) == 1) {
    yt = dpoly(x0-_spec.chan0, _spec.bkgd)*(hi-lo);
  } else {
    yt = interpol(x0, _spec.xbkgd, _spec.bkgd)*(hi-lo);
  }
  y += yt;
  return y;
}
  
% setup the mspec model function
define set_mspec(dx, dy, xr, yr) {
  variable i, k, j, m, par, x, k1, q, xmin, xmax;

  k = 0;
  for (i = 0; i < length(_spec.width); i++) {
    if (_spec.iwidth[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.alpha); i++) {
    if (_spec.ialpha[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.scale); i++) {
    if (_spec.iscale[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.bkgd); i++) {
    if (_spec.ibkgd[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.abund); i++) {
    if (_spec.iabund[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.ishifts); i++) {
    if (_spec.ishifts[i] == 0) {
      k++;
    }
  }
  for (i = 0; i < length(_spec.lines); i++) {
    for (j = 0; j < 2; j++) {
      for (m = 0; m < length(_spec.lines[i][j,*]); m++) {
	if (_spec.ilines[i][j,m] <= 0) {
	  k++;
	}
      }
    }
  }
  _spec.pid = Integer_Type[2,k];
  par = String_Type[k];
  for (i = 0; i < k; i++) {
    par[i] = sprintf("P%d", i+1);
  }

  add_slang_function("mspec", par);  
  fit_fun("mspec(1)");

  k = 1;
  for (i = 0; i < length(_spec.width); i++) {
    if (_spec.iwidth[i] == 0) {
      set_par(k, _spec.width[i]);
      k++;
    }
  }
  for (i = 0; i < length(_spec.alpha); i++) {
    if (_spec.ialpha[i] == 0) {
      set_par(k, _spec.alpha[i]);
      k++;
    }
  }
  for (i = 0; i < length(_spec.scale); i++) {
    if (_spec.iscale[i] == 0) {
      set_par(k, _spec.scale[i]);
      k++;
    }
  }
  for (i = 0; i < length(_spec.bkgd); i++) {
    if (_spec.ibkgd[i] == 0) {
      set_par(k, _spec.bkgd[i]);
      k++;
    }
  }

  for (i = 0; i < length(_spec.abund); i++) {
    if (_spec.iabund[i] == 0) {
      set_par(k, _spec.abund[i], 0, 0.0, 1E30);
      k++;
    }
  }
  for (i = 0; i < length(_spec.ishifts); i++) {
    if (_spec.ishifts[i] == 0) {
      if (dx > 0) {
	set_par(k, _spec.shifts[0,i], 0, _spec.shifts[0,i]-dx, _spec.shifts[0,i]+dx);
      } else {
	set_par(k, _spec.shifts[0,i]);
      }
      k++;
    }
  }
  
  for (i = 0; i < length(_spec.lines); i++) {
    for (j = 0; j < 2; j++) {
      if (j == 1) k1 = k;
      for (m = 0; m < length(_spec.lines[i][j,*]); m++) {
	if (_spec.ilines[i][j,m] <= 0) {
	  x = _spec.lines[i][j,m];
	  if (j == 0) {
	    if (dx > 0) {
	      xmin = x - dx;
	      xmax = x + dx;
	    } else {
	      xmin = 0.0;
	      xmax = 0.0;
	    }
	    if (length(xr) == 2) {
	      if (xmin < xr[0]) xmin = xr[0];
	      if (xmax > xr[1]) xmax = xr[1];
	    }
	  } else {
	    if (dy > 0) {
	      xmin = x*(1-dy);
	      if (xmin < 0) xmin = 1e-10;
	      xmax = x*(1+dy);
	    } else {
	      xmin = 1E-30;
	      xmax = 1E30;
	    }
	    if (length(yr) == 2) {
	      if (xmin < yr[0]) xmin = yr[0];
	      if (xmax > yr[1]) xmax = yr[1];
	    }
	  }
	  set_par(k, x, 0, xmin, xmax);
	  _spec.ilines[i][5+j,m] = k;
	  _spec.pid[j,k-1] = _spec.ilines[i][4,m];
	  k++;
	}
      }
    }
    j = 1;
    for (m = 0; m < length(_spec.lines[i][j,*]); m++) {
      if (_spec.ilines[i][j,m] > 0) continue;
      if (_spec.ilines[i][j,m] < 0) {
	x = _spec.lines[i][j,m];
	q = where(_spec.pid[1,*] == -_spec.ilines[i][j,m]);
	if (length(q) > 0) {
	  q = q[0] + 1;
	  set_par_fun(k1, sprintf("_par(%d) * %10.3E", q, x/get_par(q)));
	}
      }
      k1++;
    }
  }
}

define bkgd(s, x) {
  variable n, xmin, xmax, dx, xb, yb, w, ny, yx, y, ye, nq;
  variable k, i, j, p, ymin, ymax, dy, nb, b, d, db, bm, nt, i0;

  if (x[0] < 0) {
    n = -int(x[0]);
    xmin = double(x[1]);
    xmax = double(x[2]);
    dx = (xmax-xmin)/(n-1);
    xb = [xmin:xmax+dx*0.1:dx];
  } else {
    xb = @x;
    n = length(xb);
  }
  yb = Double_Type[n];
  
  i0 = 0;
  for (i = 1; i < n; i++) {
    w = where(s.bin_lo >= xb[i-1] and s.bin_hi <= xb[i]);
    ny = length(w);
    yx = 0.5*(s.bin_lo[w] + s.bin_hi[w]);
    y = s.value[w];
    ye = s.err[w];
    k = array_sort(y);
    ymax = y[k[ny/2]];
    ymin = y[k[0]];
    dy = ye[k[ny/4]]*0.05;
    nb = int((ymax-ymin)/dy+1);
    dy = (ymax-ymin)/(nb-1);
    b = [ymin:ymax+dy*0.1:dy];
    nt = Integer_Type[nb];
    nq = Integer_Type[nb];
    db = Double_Type[nb];
    if (i0 == 0) {
      for (p = 0; p < nb; p++) {
	d = y - b[p];
	j = where(d < -3*ye);
	nq[p] = length(j);
	j = where(abs(d) < 3*ye);
	nt[p] = length(j);
	if (nt[p] > 0) {
	  db[p] = sum((d[j]/ye[j])^2)/nt[p];
	} else {
	  db[p] = 1E30;
	}
      }
    } else {
      for (p = 0; p < nb; p++) {
	bm = yb[i-1] + (yx-xb[i-1])*(b[p]-yb[i-1])/(xb[i]-xb[i-1]);
	d = y - bm;
	j = where(d < -3*ye);
	nq[p] = length(j);
	j = where(abs(d) < 3*ye);
	nt[p] = length(j);
	if (nt[p] > 0) {
	  db[p] = sum((d[j]/ye[j])^2)/nt[p];
	} else {
	  db[p] = 1E30;
	}
      }
    }
    j = where(nq == min(nq));
    b = b[j];
    nt = nt[j];
    db = db[j];
    j = where(nt >= max(nt)*0.9);
    b = b[j];
    nt = nt[j];
    db = db[j];
    j = where(db == min(db));
    j = j[0];
    yb[i] = b[j];
    if (i0 == 0) {
      yb[0] = yb[i];
      dy = yb[0]/1e2;
      b = [-yb[0]:yb[1]+0.1*dy:dy];
      nb = length(b);
      nt = Integer_Type[nb];
      nq = Integer_Type[nb];
      db = Double_Type[nb];
      for (p = 0; p < nb; p++) {
	bm = yb[0]+b[p] - (yx-xb[0])*2.0*b[p]/(xb[1]-xb[0]);
	d = y - bm;
	j = where(d < -3*ye);
	nq[p] = length(j);
	j = where(abs(d) < 3*ye);
	nt[p] = length(j);
	if (nt[p] > 0) {
	  db[p] = sum((d[j]/ye[j])^2)/nt[p];
	} else {
	  db[p] = 1E30;
	}
      }  
      j = where(nq == min(nq));
      b = b[j];
      nt = nt[j];
      db = db[j];
      j = where(nt >= max(nt)*0.9);
      b = b[j];
      nt = nt[j];
      db =db[j];
      j = where(db == min(db));
      j = j[0];
      yb[0] += b[j];
      yb[1] -= b[j];
      i0 = 1;
      i--;
    }
  }
    
  return (xb, yb);
}

% estimate the background.
define set_bkgd(s, x) {
  variable i, j;

  if (x == NULL) {
    if (_spec.xbkgd == NULL) {
      x = [-length(_spec.bkgd), s.bin_lo[0], s.bin_hi[-1]];
    } else {
      x = @_spec.xbkgd;
    }
  }
  (_spec.xbkgd, _spec.bkgd) = bkgd(s, x);

  for (i = 0; i < length(_spec.xbkgd); i++) {
    j = where(s.bin_lo <= _spec.xbkgd[i] and s.bin_hi > _spec.xbkgd[i]);
    if (length(j) == 0) {
      if (_spec.xbkgd[i] < s.bin_lo[0]) j = 0;
      else j = length(s.bin_lo)-1;
    } else {
      j = j[0];
    }
    _spec.bkgd[i] /= s.bin_hi[j]-s.bin_lo[j];
  }
  _spec.ibkgd = Integer_Type[length(_spec.bkgd)];
}

% estimate the initial ion abundances.    
define set_abund(nk, c0, c1) {
  variable i, k, d, xd, x, sig, bkgd, t, w, wr, smax, wm, a, b, c, m;

  d = get_data_counts(1);
  xd = 0.5*(d.bin_lo + d.bin_hi);
  x = xd - _spec.chan0;
  if (_spec.xwidth == NULL or length(_spec.xwidth)==1) {
    sig = dpoly(x, _spec.width);
    if (length(_spec.alpha) > 0) {
      sig *= 1.0 + dpoly(x, _spec.alpha);
    }
  } else {
    sig = interpol(xd, _spec.xwidth, _spec.width);
    if (length(_spec.alpha) > 0) {
      sig *= 1.0 + interpol(xd, _spec.xwidth, _spec.alpha);
    }
  }
  sig *= 3.0;
  if (_spec.xbkgd == NULL or length(_spec.xbkgd) == 1) {
    bkgd = dpoly(x, _spec.bkgd);
  } else {
    bkgd = interpol(xd, _spec.xbkgd, _spec.bkgd);
  }
  bkgd *= d.bin_hi-d.bin_lo;
  k = where(_spec.ions == nk);
  if (length(k) != 1) return -1;
  k = k[0];
  wr = where(_spec.lines[k][0,*] > c0 and 
	     _spec.lines[k][0,*] < c1);
  smax = max(_spec.lines[k][1,wr]);
  wm = where(_spec.lines[k][1,wr] > 0.1*smax);
  a = 0.0; 
  b = 0.0;
  for (i = 0; i < length(wm); i++) {
    m = wr[wm[i]];
    w = where(xd < _spec.lines[k][0,m]+sig and 
	      xd > _spec.lines[k][0,m]-sig);
    if (length(w) == 0) continue;      
    t = sum(d.value[w] - bkgd[w]);
    if (t <= 0) continue;
    a += t;
    if (_spec.eff != NULL) {
      c = interpol_points(_spec.lines[k][2,m], _spec.eff[0,*], _spec.eff[1,*]);
      c = c[0];
    } else {
      c = 1.0;
    }
    b += _spec.lines[k][1,m]*_spec.lines[k][3,m]*c;
  } 
  _spec.abund[k] = a/b;
}

define fit_abund() {
  variable ilines, i, w;

  ilines = Array_Type[length(_spec.ilines)];
  for (i = 0; i < length(_spec.ilines); i++) {
    ilines[i] = @_spec.ilines[i];
  }

  _spec.iabund[*] = 0;
  w = where(_spec.ishifts == 0);
  if (length(w) > 0) {
    for (i = 0; i < length(_spec.ilines); i++) {
      _spec.ilines[i][0,*] = 1;
      _spec.ilines[i][1,*] = 1;
    }
    set_mspec(0, 0, 0, 0);
    () = fit_counts;
    for (i = 0; i < length(_spec.ilines); i++) {
      _spec.ilines[i] = @ilines[i];
    }
  }

  for (i = 0; i < length(_spec.ilines); i++) {
    _spec.ilines[i][1,*] = 1;
  }
  set_mspec(0, 0, 0, 0);
  () = fit_counts;

  _spec.iabund[*] = 1;
  for (i = 0; i < length(_spec.ilines); i++) {
    _spec.ilines[i] = @ilines[i];
  }
  return get_model_counts(1);
}

define fit_spec() {  
  set_mspec(0, 0, 0, 0);
  () = fit_counts;
  return get_model_counts(1);
}

define set_range() {
  variable n, x0, y0, x1, y1;

  if (_NARGS == 0) {    
    cursor(&x0, &y0);
    vmessage("%10.3E %10.3E\n", x0, y0);
    cursor(&x1, &y1);
    vmessage("%10.3E %10.3E\n", x1, y1);
    notice(1, min([x0,x1]), max([x0,x1]));
  } else {
    n = _NARGS;
    while (n > 0) {
      x1 = ();
      x0 = ();
      n -= 2;
      notice(1, x0, x1);
    }
  }
}
    
% fit the wavelength calib. coeff. using polynomial with order n-1.
define set_calib(n) {
  variable i, j, k, x, y, p, s;

  if (typeof(n) == Array_Type) {
    _spec.calib = @n;    
    return;
  }
  k = 0;
  for (i = 0; i < length(_spec.lines); i++) {
    for (j = 0; j < length(_spec.lines[i][0,*]); j++) {
      k++;
    }
  }
  x = Double_Type[k];
  y = Double_Type[k];
  _spec.calib = Double_Type[2,k];
  k = 0;
  for (i = 0; i < length(_spec.lines); i++) {
    for (j = 0; j < length(_spec.lines[i][0,*]); j++) {
      x[k] = _spec.lines[i][0,j];
      y[k] = _spec.lines[i][2,j];
      k++;
    }
  }
  k = array_sort(x);
  _spec.calib[0,*] = x[k];
  _spec.calib[1,*] = y[k];
}

define zoom_replot(x, y, xr, yr, xs, ys) {
  variable w, s, i;

  color(1);
  if (typeof(x) == Struct_Type) {
    hplot(x);
    if (xs != NULL) {
      for (i = 0; i < length(xs); i++) {
	if (typeof(xs[i]) == Array_Type) {
	  if (length(xs[i]) > 0) {
	    if (max(xs[i]) >= xr[0] and min(xs[i]) <= xr[1] and
		max(ys[i]) >= yr[0] and min(ys[i]) <= yr[1]) {
	      color(i+2);
	      oplot(xs[i], ys[i]);
	    }
	  }
	} else if (typeof(xs[i]) == Struct_Type){
	  if (xs[i].bin_hi[-1] >= xr[0] and xs[i].bin_lo[0] <= xr[1] and
	      max(xs[i].value) >= yr[0] and 
	      min(xs[i].value) <= yr[1]) {
	    color(i+2);
	    ohplot(xs[i]);
	  }
	} else {
	  color(2);
	  if (xs[i] < xr[1] and xs[i] > xr[0]) {
	    oplot([xs[i], xs[i]], [1e-10*ys[i], ys[i]]);
	  }
	}
      }
    }
  } else {
    connect_points(0);
    plot(x, y);
    if (xs != NULL) {
      if (typeof(xs[0]) != Array_Type) {
	for (i = 0; i < length(xs); i += 5) {
	  connect_points(1);
	  color(3);
	  w = [i:i+4];
	  s = where(xs[w] < xr[1] and xs[w] > xr[0] and
		    ys[w] < yr[1] and ys[w] > yr[0]);
	  if (length(s) > 0) {
	    oplot(xs[[i:i+4]], ys[[i:i+4]]);
	  }
	}
      } else {
	for (i = 0; i < length(xs); i++) {
	  connect_points(1);
	  color(3);
	  if (length(xs[i]) > 0) {
	    if (max(xs[i]) >= xr[0] and min(xs[i]) <= xr[1] and
		max(ys[i]) >= yr[0] and min(ys[i]) <= yr[1]) {
	      oplot(xs[i], ys[i]);
	    }
	  }
	}
      }
    }
  }
}

define zoom_plot(x, y, xr0, yr0, c, xs, ys) {
  variable x0, y0, x1, y1, ch, xr, yr, i, w, s, zm;

  if (c != 'h' and 
      c != 'H' and
      c != 'v' and
      c != 'V' and
      c != 'z' and
      c != 'Z') return -1;
  
  if (c == 'h' or c == 'H') zm = "H-ZOOM";
  if (c == 'v' or c == 'V') zm = "V-ZOOM";
  if (c == 'z' or c == 'Z') zm = "ZOOM";
  xr = @xr0;
  yr = @yr0;
  (xr[0], xr[1], yr[0], yr[1]) = _pgqwin();
  while (1) {
    _pgsci(2);
    color(2);
    xylabel(xr[0]+0.1*(xr[1]-xr[0]), yr[1]+0.035*(yr[1]-yr[0]), zm);
    () = _pgband(7, 0, 0, 0, &x0, &y0, &ch);
    if (ch == 's' or ch == 'S') {
      resize(0.0, 100.0);
      (, , ,s) = _pgqvsz(2);
      resize(0.0, 0.01);
      (, w, ,) = _pgqvsz(2);
      resize(w*0.1, s/w);
      zoom_replot(x, y, xr, yr, xs, ys);
      continue;
    }
    if (ch == 'l' or ch == 'L') {
      limits;
      zoom_replot(x, y, xr, yr, xs, ys);
      (xr0[0], xr0[1], yr0[0], yr0[1]) = _pgqwin();
      xr = @xr0;
      yr = @yr0;
      continue;
    }
    if (ch == 'r' or ch == 'R') {
      xrange(xr0[0], xr0[1]);
      yrange(yr0[0], yr0[1]);
      xr = @xr0;
      yr = @yr0;
    } else {
      if (ch == 'q' or ch == 'Q') break;
      () = _pgband(2, 0, x0, y0, &x1, &y1, &ch);
      if (ch == 'q' or ch == 'Q') break;
      if (x0 > x1) (x0, x1) = (x1, x0);
      if (y0 > y1) (y0, y1) = (y1, y0);
      if (x0 < xr0[0]) x0 = xr0[0];
      if (x1 > xr0[1]) x1 = xr0[1];
      if (y0 < yr0[0]) y0 = yr0[0];
      if (y1 > yr0[1]) y1 = yr0[1];
      if (c != 'v' and c != 'V') {
	xrange(x0, x1);
	xr = [x0, x1];
      } 
      if (c != 'h' and c != 'H') {
	yrange(y0, y1);
	yr = [y0, y1];
      }
    }
    zoom_replot(x, y, xr, yr, xs, ys);
  }
  zoom_replot(x, y, xr, yr, xs, ys);
  return 0;
}

define prepare_lines(fn, s, xr, yr, nele, iylog) {
  variable f, k, x, y, ch, i, xs, ys;

  f = NULL;
  k = 0;
  xs = Double_Type[0];
  ys = @xs;
  if (stat_file(fn) != NULL) {
    () = system(sprintf("cp %s %s.bak", fn, fn));
  }
  while (1) {
    _pgsci(2);
    () = _pgband(7, 0, 0, 0, &x, &y, &ch);
    i = zoom_plot(s, NULL, xr, yr, ch, xs, ys);
    if (i >= 0) {
      () = _pgband(7, 0, 0, 0, &x, &y, &ch);
    }
    if (ch == 'q' or ch == 'Q') break;
    if (ch == 'n' or ch == 'N') {
      xs = Double_Type[0];
      ys = @xs;
      color(1);
      hplot(s);
      k = 0;
      if (stat_file(fn+".bak") != NULL) {
	() = system(sprintf("cp %s.bak %s", fn, fn));
      }
      vmessage("restart line selection");
      () = fflush(stdout);
      continue;
    }
    if (iylog) y = 10^y;
    vmessage("x= %10.3E y= %10.3E", x, y);
    () = fflush(stdout);
    if (k == 0) {
      f = fopen(fn, "w");  
    }
    () = fprintf(f, "%3d %2d %5d %5d %11.4E %11.4E %11.4E 1.0 0.0 0 0\n",
		 k, nele, 0, k+1, x, y, x);
    color(2);
    oplot([x,x], [1e-10*y, y]);
    xs = [xs, x];
    ys = [ys, y];
    k++;
  } 
  color(1);
  if (k > 0) () = fclose(f);
}

define xrsrsp(rfn, w, a, nb) {
  variable f, i;

  f = fopen(rfn, "w");
  () = fprintf(f, "chan0 0.0\n");
  () = fprintf(f, "width %10.3E\n", w);
  () = fprintf(f, "iwidth 0\n");
  if (a >= 0) {
    () = fprintf(f, "alpha %10.3E\n", a);
    () = fprintf(f, "ialpha 0\n");
  }
  if (nb != 2) {
    () = fprintf(f, "bkgd  0.0\n");
    () = fprintf(f, "ibkgd  0\n");
  } else {
    () = fprintf(f, "bkgd  0.0 0.0\n");
    () = fprintf(f, "ibkgd  0 0\n");
  }
  () = fclose(f);
}

define fit_peaks(s, flines, frsp, oflines, ofrsp, xr) {
  variable i, m;

  set_fit_method("marquardt;max_loops=1000;tol=1E-8;delta=1E-10");
  set_fit_statistic("chisqr");

  delete_data(all_data);
  () = define_counts(s);

  ignore(1);
  for (i = 0; i < length(xr); i += 2) {
    notice(1, xr[i], xr[i+1]);
  }
  
  read_response(frsp);
  read_lines(flines, 0, 0);
  set_mspec(0, 0, 0, 0);
  _spec.abund[*] = 1.0;
  _spec.iabund[*] = 1;
 
  () = fit_counts();
  
  write_lines(oflines, 0);
  if (ofrsp != NULL) {
    write_response(ofrsp);
  }
  m = get_model_counts(1);
  hplot(s);
  ohplot(m);
  return m;
}

define set_kconstrains(r0, v0, a, xmin, xmax, ia, ksp) {
  variable s, n, x, y, dx, w, i0, i1, i, j, ik, p, xc, yc;
  variable xp, yp, xp0, yp0, q, k, m, r, v;

  r = Array_Type[length(r0)];
  for (i = 0; i < length(r); i++) {
    r[i] = @(r0[i]);
  }
  v = @v0;
  if (ksp > 1) {
    w = where(_spec.lines[2][1,*] == max(_spec.lines[2][1,*]));
    w = w[0];
    xp0 = _spec.ilines[2][5, w];
    yp0 = _spec.ilines[2][6, w];
    if (_spec.lines[2][0,w] >= xmax or _spec.lines[2][0,w] <= xmin) {
      if (xp0) freeze(xp0);
      if (yp0) freeze(yp0);
    } else {
      if (xp0) thaw(xp0);
      if (yp0 and not ia) thaw(yp0);
    }
    for (i = 0; i < length(_spec.lines[2][1,*]); i++) {
      if (i == w) continue;
      xp = _spec.ilines[2][5, i];
      yp = _spec.ilines[2][6, i];
      if (xp and xp0) {
	xc = (_spec.lines[2][0, i] - _spec.lines[2][0, w]);
	set_par(xp, 0, 0, -1, 1e31);
	set_par_fun(xp, sprintf("_par(%d)+(%11.4E)", xp0, xc));
      }
      if (yp and yp0) {
	yc = _spec.lines[2][1, i] / _spec.lines[2][1, w];
	set_par_fun(yp, sprintf("_par(%d)*(%11.4E)", yp0, yc));
      }
    }

    w = where(r[0] < 3);  
    for (i = 0; i < length(r); i++) {
      r[i] = r[i][w];
    }
    v = v[w];
  }

  s = array_sort(v);
  n = length(s);
  x = v[s];
  y = r[5][s];
  dx = x[[1:]] - x[[0:n-2]];
  w = where(dx > a);
  w = [w, n-1];
  i0 = 0;
  for (i = 0; i < length(w); i++) {
    i1 = w[i];
    k = [i0:i1];
    ik = where(y[k] == max(y[k]));
    ik = k[ik[0]];
    m = r[0][s[ik]];
    if (m > 2) continue;
    m = where(_spec.ions == m);
    if (length(m) != 1) continue;
    m = m[0];
    q = where(_spec.ilines[m][2,*] == r[1][s[ik]] and
	      _spec.ilines[m][3,*] == r[2][s[ik]]);
    if (length(q) != 1) continue;
    q = q[0];
    xp0 = _spec.ilines[m][5, q];
    yp0 = _spec.ilines[m][6, q];
    if (_spec.lines[m][0,q] >= xmax or _spec.lines[m][0,q] <= xmin) {
      if (xp0) freeze(xp0);
      if (yp0) freeze(yp0);
    } else {
      if (xp0) thaw(xp0);
      if (yp0 and not ia) thaw(yp0);
    }
    if (xp0 == 0 and yp0 == 0) continue;
    for (j = i0; j <= i1; j++) {
      if (j == ik) continue;
      m = r[0][s[j]];
      m = where(_spec.ions == m);
      if (length(m) != 1) continue;
      m = m[0];
      q = where(_spec.ilines[m][2,*] == r[1][s[j]] and
		_spec.ilines[m][3,*] == r[2][s[j]]);
      if (length(q) != 1) continue;
      q = q[0];
      xp = _spec.ilines[m][5, q];      
      yp = _spec.ilines[m][6, q];
      if (_spec.lines[m][0,q] >= xmax or _spec.lines[m][0,q] <= xmin) {
	if (xp) freeze(xp);
	if (yp) freeze(yp);
      } else {
	if (xp) thaw(xp);
	if (yp and not ia) thaw(yp);
      }
      if (xp > 0 and xp0 > 0) {
	p = get_par_info(xp);
	if (p.freeze == 0) {
	  xc = v[s[j]] - v[s[ik]];
	  set_par(xp, 0, 0, -1, 1e31);
	  set_par_fun(xp, sprintf("_par(%d)+(%11.4E)", xp0, xc));
	}
      }
    }
    i0 = i1+1;
  }
}
      
define fit_kspec(ofn, s, a, uv1, uv2, rfn, maxiter) {
  variable lfn, afn, ldb, r, k, fn, i, iter, x0, x1, b, xmin, xmax;
  variable f, m, nr, v, ia, w, x, y, c, dx, sm, t, vsp, km;
  variable ue, ymax, ix, uv, us, n1, n2, iu, j, q, ksp, sig, xs, ys; 

  if (typeof(a) == Integer_Type) {
    a = const->AtomicSymbol[a];
  }
  vsp = 0;
  if (uv1 == NULL or uv2 == NULL) {
    fn = a;
    if (stat_file(fn) != NULL) {
      r = rcols(fn, [1,2,3,0,6,5], ["I","I","I","I","F","F"], '#', 0, -1, 1);
      (uv, us, ue) = readcol(fn, 5, 6, 7);
      m = where(uv != ue);
      if (length(m) == 0) {
	maxiter = 1;
	vsp = 1;
      }
    } else {
      vmessage("cannot find the input line list, no fitting done");
      () = fflush(stdout);
      return;
    }
    ksp = 0;
  } else {
    ldb = getenv("HOME")+"/lib/isis/klines/";
    fn = ldb + a + sprintf("%02da.ln", 1);
    ksp = 1;
    if (typeof(uv1) == String_Type) {
      if (stat_file(uv1) != NULL) {
	uv1 = readcol(uv1, 5);
      } else {
	uv1 = Double_Type[0];
      }
    }
    if (typeof(uv2) == String_Type) {
      if (stat_file(uv2) != NULL) {
	uv2 = readcol(uv2, 5);
      } else {
	uv2 = Double_Type[0];
      }
    }
    n1 = length(uv1)+2;
    n2 = length(uv2)+2;
    r = Array_Type[6];
    uv = [uv1, uv2];
    ue = @uv;
    if (length(uv2) == 0) {
      ksp = 1;
    } else {
      ksp = 3;
    }
    for (k = 1; k <= ksp; k++) {
      fn = ldb + a + sprintf("%02da.ln", k);
      x = rcols(fn, [0,1,2,3,4,6], ["I","I","I","I","F","F"], '#', 0, -1, 1);
      if (k == 1) {
	w = where(x[3]/100 < n1);
      } else if (k == 2) {
	w = where(x[3]/100 < n2);
      } else {
	w = where(x[5] >= 0.1*max(x[5]));
      }
      for (i = 0; i < 6; i++) {
	if (k == 1) {
	  r[i] = x[i][w];
	} else if (k == 2) {
	  r[i] = [r[i], x[i][w]];
	} else {
	  r[i] = [r[i], x[i][w]];
	}
      }
    }

    iu = Integer_Type[length(uv)];
    for (i = 0; i < length(uv1); i++) {
      k = where(r[0] == 1 and r[3]/100 == i+2);
      w = where(r[5][k] == max(r[5][k]));
      w = k[w[0]];
      ue[i] = r[4][w];
      iu[i] = w;
    } 
    for (i = 0; i < length(uv2); i++) {
      k = where(r[0] == 2 and r[3]/100 == i+2);
      w = where(r[5][k] == max(r[5][k]));
      w = k[w[0]];
      ue[i+length(uv1)] = r[4][w];
      iu[i+length(uv1)] = w;
    }
  }

  lfn = ofn + ".lni";
  afn = ofn + ".lno";
  x = @uv;
  y = @ue;
  set_fit_method("marquardt;max_loops=1000;tol=1e-8;delta=1E-10");
  set_fit_statistic("chisqr;sigma=data");
  %set_fit_method("marquardt");
  %set_fit_statistic("chisqr");
  if (vsp == 0) {
    (c, t) = array_fit(x, y/x, NULL, [0.0, 0.0], NULL, NULL, &dpoly);
    c = [0.0, c];
  } else {
    c = [0.0, 1.0];
  }
  iter = 0;
  while (1) {
    iter++;
    if (vsp == 0) {
      vmessage("iteration %d %10.3E %10.3E", iter, c[1], c[2]);
    } else {
      vmessage("iteration %d %10.3E", iter, c[1]);
    }
    () = fflush(stdout);
    delete_data(all_data);
    () = define_counts(s); 
    read_response(rfn);
    nr = length(r[0]);
    v = ipoly(r[4], c);
    if (iter == 1) {
      sig = _spec.width[0];
    }
    if (not ksp) {
      r[5] = us*sqrt(2.0*PI)*sig/(s.bin_hi[0]-s.bin_lo[0]);      
    }
    f = fopen(lfn, "w");
    m = 0;
    for (i = 0; i < nr; i++) {
      () = fprintf(f, "%3d %2d %4d %4d %11.4E %11.4E %13.7E 1.0 0.0 0 0\n",
		   m, r[0][i], r[1][i], r[2][i], v[i], r[5][i], r[4][i]);
      m++;
    }
    () = fclose(f);

    x0 = min(v) - sig*10;
    x1 = max(v) + sig*10;
    if (ksp) {
      m = where(r[0] == 2);
      if (length(m) > 0) {
	x0 = min(v[m]) - sig*10;
      }
    }
    xnotice(1, x0, x1);
    limits;
    xrange(x0, x1);
    w = where(s.bin_lo > x0 and s.bin_hi < x1);
    ymax = max(s.value);
    yrange(-ymax*0.1, ymax*1.2);
    color(1);
    hplot(s);
    read_lines(lfn, 0, 0);
    set_fit_method("marquardt;max_loops=100;tol=1e-4;delta=1E-6");
    set_fit_statistic("chisqr;sigma=data");
    _spec.abund[*] = 1.0;
    if (ksp) {
      _spec.iabund[*] = 0;
      dx = sig*5;
      set_abund(1, uv[0]-dx, uv[0]+dx);
      if (length(uv2) > 0) {
	set_abund(2, uv[length(uv1)]-dx, uv[length(uv1)]+dx);
	_spec.abund[2] = _spec.abund[1];
      }
      for (i = 0; i < length(_spec.ilines); i++) {
	_spec.ilines[i][1,*] = 1;
      }
      x = ue[0];
      dx = sig*0.5;
      set_mspec(dx, 0, 0, 0);
      set_par(1, sig, 0, 0.1*sig, 5.0*sig);
      xmin = x0;
      xmax = x1;
      xnotice(1, xmin, xmax);
      set_kconstrains(r, v, sig*2, xmin, xmax, 1, ksp);
      () = fit_counts();    
      _spec.iabund[*] = 1;
      sig = _spec.width[0];
      x = Double_Type[length(r[0])];
      m = 0;
      for (i = 0; i < length(_spec.ilines); i++) {
	for (j = 0; j < length(_spec.ilines[i][0,*]); j++) {
	  q = _spec.ilines[i][5,j];
	  if (q) {
	    q = get_par_info(q);
	    if (q.freeze == 0) {
	      x[m] = _spec.lines[i][0,j];
	      m++;
	    }
	  }
	}
	_spec.ilines[i][1,*] = 0;
      }
      x = x[[0:m-1]];
      dx = sig*10.0;
    } else {
      x = @ue;
      m = length(x);
      dx = sig*3.0;
      _spec.iabund[*] = 1;
    }
    ix = array_sort(x);
    x = x[ix];
    y = x[[1:m-1]] - x[[0:m-2]];
    ix = where(y > sig*5);
    xmin = Double_Type[length(ix)+1];
    xmax = @xmin;
    xmin[0] = x0;
    xmax[-1] = x1;
    for (i = 1; i <= length(ix); i++) {
      b = x[ix[i-1]+1] - x[ix[i-1]];
      b *= 0.5;
      if (b > sig*5) b = sig*5;
      xmax[i-1] = x[ix[i-1]]+b;
      xmin[i] = x[ix[i-1]+1]-b;
    }
    set_mspec(dx, 0, 0, 0);
    set_par(1, sig, 0, 0.1*sig, 5.0*sig);
    b = 0.0;
    xs = Struct_Type[length(xmin)];
    ys = @xs;
    for (ix = 0; ix < length(xmin); ix++) {
      xnotice(1, xmin[ix], xmax[ix]);
      set_kconstrains(r, v, sig*2, xmin[ix], xmax[ix], 0, ksp);      
      () = fit_counts();
      sm = get_model_counts(1);
      m = (ix mod 7) + 2;
      color(m);
      m = where(sm.bin_lo > xmin[ix] and sm.bin_hi < xmax[ix]);
      sm.bin_lo = sm.bin_lo[m];
      sm.bin_hi = sm.bin_hi[m];
      sm.value = sm.value[m];
      ohplot(sm);
      xs[ix] = sm;      
      ys[ix] = NULL;
      b += _spec.width[0];
    }
    b /= length(xmin);
    if (not ksp) sig = b;
    color(1);
    xylabel(x1-0.4*(x1-x0), ymax, sprintf("iteration %d", iter));
    if (vsp == 1) {
      write_lines(afn, -1);
    } else {
      write_lines(afn, 0);
    }
    write_response(rfn+".out");
    v = rcols(afn, [1,4,5,10], ["I","F","F","I"], ' ', 0, -1, 1);
    b = @(v[2]);
    for (i = 0; i < length(_spec.ions); i++) {
      ix = where(v[0] == _spec.ions[i]);
      if (length(ix) > 0) {
	b[ix] = v[2][ix]*_spec.abund[i];
      }
    }
    if (vsp == 0) {
      x = v[1];
      y = r[4];    
      if (ksp) {
	t = iu;
	x = x[t];
	y = y[t];
	b = b[t]/max(b[t]);
      } else {
	r[5] = @b;
	b = b/max(b);
      }
      set_fit_method("marquardt;max_loops=1000;tol=1e-8;delta=1E-10");
      set_fit_statistic("chisqr;sigma=data");
      (c, t) = array_fit(x, y/x, b,
			 [0.0, 0.0], NULL, NULL, &dpoly);
      c = [0.0, c];
    } else {
      c = [0.0, 1.0];
    }
    if (iter >= maxiter) break;
    sleep(1.0);
%    vmessage("Type q to stop iteration, Return to continue ...");
%    () = fread(&t, Char_Type, 1, stdin);
%    if (t == 'q' or t == 'Q') break;
  }
  x0 = [0.0, 0.0];
  x1 = [0.0, 0.0];
  (x0[0],x0[1],x1[0],x1[1]) = _pgqwin();
  while (1) {
    color(1);
    _pgsci(2);
    () = _pgband(7, 0, 0, 0, &x, &y, &c);
    if (c == 'q' or c == 'Q') break;
    () = zoom_plot(s, NULL, x0, x1, c, xs, ys);
  }
}

define comblines(ofn, lnf, lab, cts) {
  variable i, f, a, r, w, nw, j, k;

  f = fopen(ofn, "w");
  for (i = 0; i < length(lnf); i++) {
    a = rcols(lnf[i], [2], ["F"],' ', 0, 2, 1);
    a = a[0]; 
    r = rcols(lnf[i], [1,2,3,4,5,6,10], ["I","I","I","F","F","F","I"],
	      ' ', 0, -1, 1);
    w = array_sort(r[3]);
    for (j = 0; j < 6; j++) {
      r[j] = r[j][w];
    }
    r[4] *= a[r[0]-1];
    w = where(r[4] > cts and r[0] < 3);
    nw = length(w);
    for (j = 0; j < nw; j++) {
      k = w[j];
      () = fprintf(f, "%5s %d.%02d.%02d %13.7E %13.7E %12.5E %-s\n",
		   lab[i], r[0][k], r[1][k], 
		   r[2][k], r[5][k], r[3][k], r[4][k], lnf[i]);
    }
  }
  () = fclose(f);
}

% plot labels on a spectrum.
% tx0, the x coord of line position.
% tn0, line labels.
% xr, the x-axis resolution, lines separated by < xr are treated as one.
% dx, the x- separation of labels.
% sp, the spectrum.
% ys, the y-value where the label lines start.
% y1, the y-value where the label lines start to deviate from vertical.
% y2, the y-value where the label lines stop, and label texts begin.
% wfm, the format string for the wavelengths. NULL if no wavelengths.
% cline, the color of the lines.
% ctext, the color of the texts.
% csz, the charsize.
define plot_labels(tx0, tn0, xr, dx, sp, ys, y1, y2, wfm, cline, ctext, csz) {
  variable s, n0, nm, nt, k, tn, tx, tx1, q, j, dd, w, y0, xs;

  s = array_sort(tx0);
  tn0 = tn0[s];
  tx0 = tx0[s];
  n0 = length(tx0);
  nm = Integer_Type[n0];
  nm[0] = 1;
  nt = 0;
  for (k = 1; k < n0; k++) {
    if (tx0[k] - tx0[k-1] < xr) {      
      nm[nt]++;
      continue;
    }
    nt++;
    nm[nt] = 1;
  }  
  nm = nm[[0:nt]];
  nt++;
  tn = Array_Type[nt];
  tx = Double_Type[nt];
  q = 0;
  for (k = 0; k < nt; k++) {
    tn[k] = String_Type[nm[k]];
    tx[k] = tx0[q];
    for (j = 0; j < nm[k]; j++) {
      tn[k][j] = tn0[q];
      if (j > 0) {
	w = where(tn[k][[0:j-1]] == tn0[q]);	
	if (length(w) > 0) {
	  tn[k][j] = "";
	}
      }
      q++;
    }
    w = where(tn[k] != "");
    tn[k] = tn[k][w];
    nm[k] = length(w);
  }
  tx1 = @tx;
  for (k = 1; k < nt; k++) {
    if (tx1[k] - tx1[k-1] < dx) {
      dd = dx - 0.9*(tx1[k] - tx1[k-1]);
      tx1[k] += 0.7*dd;
      for (q = k-1; q >= 0; q--) {
	tx1[q] -= 0.3*dd;
	if (q > 0) {
	  if (tx1[q] - tx1[q-1] > dx) break;
	}
      }
    }
  }
  
  for (k = 0; k < nt; k++) {
    color(cline);
    w = where(sp.bin_lo >= tx[k]-xr and sp.bin_hi <= tx[k]+xr);	
    if (length(w) == 0) {
      w = where(sp.bin_lo <= tx[k] and sp.bin_hi > tx[k]);
      if (length(w) == 0) continue;
      y0 = sp.value[w[0]] + ys;
    } else {
      y0 = max(sp.value[w]) + ys;
    }

    set_line_width(1);
    oplot([tx[k], tx[k], tx1[k]], [y0, y1, y2]);
    s = tn[k][0];
    for (j = 1; j < nm[k]; j++) {
      s = s + "," + tn[k][j];
    }
    if (wfm != NULL) {
      s = s + sprintf("("+wfm+")", tx[k]);
    }
    color(ctext);
    set_line_width(1);
    _pgsch(csz);
    charsize(csz);
    (,xs) = _pgqcs(4);
    xylabel(tx1[k]+4*xs, y2*1.02, s, 90);
  }
}
    
define lifetime_fit(lo, hi, p) {
  variable np, nt, y, i, a, t;

  np = length(p);
  nt = (np-1)/2;
  y = p[0]*(hi-lo);
  for (i = 0; i < nt; i++) {
    a = p[i+1+nt];
    t = p[i+1];
    y += a*t*exp(-(lo-1.0)/t)*(1.0-exp(-(hi-lo)/t));
  }

  return y;
}

define setup_lifetime_fun(p) {
  variable pn, nt, i;
  
  nt = (length(p)-1)/2;
  pn = String_Type[nt*2+1];
  pn[0] = "bkgd";
  for (i = 0; i < nt; i++) {
    pn[i+1] = sprintf("T%d", i);
    pn[i+1+nt] = sprintf("A%d", i);
  }
  add_slang_function("lifetime", pn);
  fit_fun("lifetime(1)");
  for (i = 0; i < length(p); i++) {
    set_par(i+1, p[i]);
  }
  freeze(1);
}

% return channel relative to lam0
define gcalibc(lam, p) {
  variable t, c, ds, x, ti, t0, a, lam0;

  ds = p[0]; % line spacing
  lam0 = p[1]; % reference wavelength;
  ti = p[2]; % beta0, incident angle
  a = p[3]; % alpha
  x = p[[4:]];
  
  ds *= 1e-7;
  ti *= PI/180.0;
  ti = cos(ti);
  a *= PI/180.0;
  
  t0 = acos(ti - lam0*ds);
  t = acos(ti - lam*ds);
  t = t - t0;

  c = sin(t)/sin(a-t);
  c = c*dpoly(c, x);

  return c;
}

% return lam from channel
define gcalibw(c, p) {
  variable t, lam, ds, x, ti, t0, a, lam0;
  
  ds = p[0]; % line spacing
  lam0 = p[1]; % reference wavelength;
  ti = p[2]; % beta0, incident angle
  a = p[3]; % alpha
  x = p[[4:]];
  
  ds *= 1e-7;
  ti *= PI/180.0;
  ti = cos(ti);
  a *= PI/180.0;

  t0 = acos(ti - lam0*ds);

  x = ipoly(c, [0.0, x]);
  t = x*sin(a)/(1+x*cos(a));
  t = atan(t);
  t = t0 + t;
  lam = (ti - cos(t))/ds;
  return lam;
}

% return channel relative to lam0
define ccalibc(lam, p) {
  variable t, t0, c, d2, a, b, lam0;

  d2 = p[0]; % 2d spacing.
  lam0 = p[1]; % reference wavelength;
  a = p[2]; %alpha
  b = p[[3:]]; 
  t0 = asin(lam0/d2);
  t = asin(lam/d2);
  t = t-t0;
  a *= PI/180.0;
  t = sin(t)/sin(a-t);
  c = t*dpoly(t, b);
  return c;
}

% return lam from channel
define ccalibw(c, p) {
  variable t, t0, lam, d2, a, b, lam0, x;

  d2 = p[0]; % 2d spacing.
  lam0 = p[1]; % reference wavelength;
  a = p[2]; %alpha
  b = p[[3:]]; 
  t0 = asin(lam0/d2);

  a *= PI/180.0;

  x = ipoly(c, [0.0,b]);
  t = x*sin(a)/(1+x*cos(a));
  t = atan(t);
  t = t0 + t;
  lam = d2*sin(t);
  
  return lam;
}

define interpol2d(x0, y0, x, y, z) {
  variable nx, i, a;

  nx = length(x);
  a = Double_Type[nx];  
  for (i = 0; i < nx; i++) {
    if (y0 < y[0]) a[i] = z[i,0];
    else if (y0 > y[-1]) a[i] = z[i,-1];
    else a[i] = interpol(y0, y, z[i,*]);
  }
  if (x0 < x[0]) return a[0];
  else if (x0 > x[-1]) return a[-1];
  else return interpol(x0, x, a);
}
