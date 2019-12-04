variable _mfun = struct {fun, arg};

define mfun_fit(lo, hi, p) {
  variable y;
  
  y = Double_Type[1];
  y[0] = (@_mfun.fun)(p, _mfun.arg);
  return y;
}

static define mfun_stat(d, y, w) {
  return (y, y[0]);
}

static define mfun_report(s, n, nv) {
  return sprintf(" mfun stat = %0.4g\n", s);
}

define set_minim(p, mfun, arg) {
  variable id, n, pname, i;

  _mfun.fun = mfun;
  _mfun.arg = arg;

  id = define_counts([1.0], [2.0], [0.0], [0.0]);
  
  n = length(p);
  pname = String_Type[n];
  for (i = 0; i < n; i++) {
    pname[i] = sprintf("P%d", i);
  }
  add_slang_function("mfun", pname);

  fit_fun("mfun(1)");
  for (i = 0; i < n; i++) {
    set_par(i+1, p[i]);
  }

  add_slang_statistic("mfun", &mfun_stat, &mfun_report);
  set_fit_statistic("mfun");
  
  return id;
}

define get_minim() {
  variable p, i, n, s;

  s = struct {statistic, num_variable_params, num_bins};
  () = fit_counts(&s);
  n = get_num_pars();
  p = Double_Type[n];
  for (i = 0; i < n; i++) {
    p[i] = get_par(i+1);
  }

  return (p, s.statistic);
}

