define mac2unix(fn, ofn) {
  variable f, buf, i, n, a;

  f = fopen(fn, "r");
  buf = fgetslines(f);
  () = fclose(f);
  
  f = fopen(ofn, "w");
  for (i = 0; i < length(buf); i++) {
    a = strtrans(buf[i], "\r", "\n");
    () = fprintf(f, "%s", a);
  }
  () = fclose(f);
}
 
define chtext(fn, ofn, m, q) {
  variable f, buf, i, n, a, k;
  
  f = fopen(fn, "r");
  buf = fgetslines(f);
  () = fclose(f);
  
  f = fopen(ofn, "w");
  k = 0;
  for (i = 0; i < length(buf); i++) {
    if (m == 0) {
      a = strtrans(buf[i], "\r", "\n");
    } else {
      a = strtrans(buf[i], "\r", "");
    }
    if (strlen(a) == 1) {
      k++;
    } else {
      if (k > q) {
	() = fprintf(f, "\n");
      }
      () = fprintf(f, "%s", a);
      k = 0;
    }
  }
  () = fclose(f);
}
 
define rcols(fn, ic, fmt, comment, i0, i1, step) {
  variable f, buf, n, m, y, i, j, s, k, w, a, mc, v, x;

  f = fopen(fn, "r");
  buf = fgetslines(f);
  () = fclose(f);

  n = length(buf);
  if (n == 1) {
    buf = strchop(buf[0], '\r', 0);
    n = length(buf);
  }

  m = length(ic);
  if (length(fmt) != m or typeof(fmt) != Array_Type) {
    fmt = String_Type[m];
    fmt[*] = "F";
  }
  
  w = where(fmt != "R");
  mc = max(ic[w])+1;
  y = Array_Type[m];
  for (i = 0; i < m; i++) {
    if (fmt[i] == "F") {
      y[i] = Double_Type[n];
    } else if (fmt[i] == "A") {
      y[i] = String_Type[n];
    } else if (fmt[i] == "I") {
      y[i] = Integer_Type[n];
    } else if (fmt[i] == "R") {
      y[i] = String_Type[n];
    }
  }

  k = 0;
  if (i0 < 0) i0 = 0;
  if (i1 < 0) i1 = n-1;
  if (step <= 0) step = 1;
  for (i = i0; i <= i1; i += step) {
    buf[i] = strtrim(buf[i]);
    if (buf[i] == "") continue;
    if (buf[i][0] == comment) continue;
    a = strchop(strcompress(buf[i], " \t"), ' ', 0);
    if (length(a) < mc) continue;
    v = 1;
    j = 0;
    for (s = 0; s < m; s++) {
      x = ic[s];
      if (fmt[s] == "A") {
	y[s][k] = a[x];
      } else if (fmt[s] == "R") {
	if (x < length(a)) {
	  y[s][k] = strjoin(a[[x:]], " ");
	} else {
	  y[s][k] = "";
	}
      } else {
	if ((not isdigit(a[x])) and a[x][0] != '-' and a[x][0] != '.') {
	  v = 0;
	  break;
	} else {
	  if (fmt[s] == "F") {
	    y[s][k] = atof(a[x]);
	  } else {
	    y[s][k] = integer(a[x]);
	  }
	}
      }	
    }
    if (v == 1) k++;
  }
  
  if (k > 0) {
    w = [0:k-1];
    for (s = 0; s < m; s++) {
      y[s] = y[s][w];
    }
  } else {
    n = 0;
    for (i = 0; i < m; i++) {
      if (fmt[i] == "F") {
	y[i] = Double_Type[n];
      } else if (fmt[i] == "A") {
	y[i] = String_Type[n];
      } else if (fmt[i] == "I") {
	y[i] = Integer_Type[n];
      }
    }
  }

  return y;
}
