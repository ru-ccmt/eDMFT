auto GetTrialValues3(int ii, int itc, int Norder, double beta,
		     const TinyVector<Array<int,2>,4>& diag_c, const TinyVector<Array<int,1>,4>& pseud_c, const TinyVector<Array<int,1>,4>& dlttyp_c,
		     const Array<double,1>& times, const vector<Spline1D<double>>& Gms, const vector<Spline1D<double>>& Deltas,
		     const Array<double,1>& tau, const Array<double,2>& Gp, const Array<double,3>& Delt,
		     int Ndd){
  struct retVals {
    TinyVector<double,4> dPhi;
    TinyVector<double,4> Gp_trial0;
    TinyVector<double,4> Gp_trial1;
    TinyVector<double,4> Delt_trial;
    TinyVector<int,4> id_which;
    map<int,intpar> ip_ieis;
  };
  retVals r;
  double dt0 = (itc>0) ? times(itc)-times(itc-1) : times(itc)-(times(2*Norder-1)-beta);
  intpar ip0 = Interp(dt0,tau);
  int itcm1 = (itc>0) ? itc-1 : 2*Norder-1;
  intpar ip1 = Interp(times(itc+1)-times(itc),tau);
  //bool Print=(ii==2432);
  for (int ic=0; ic<Ndd; ++ic){
    int atomic_state0 = pseud_c[ic](itcm1); // itc-1 propagator
    r.Gp_trial0[ic] = Gms[atomic_state0](ip0); // interpolate Gp for new time difference
    int atomic_state1 = pseud_c[ic](itc);   // itc propagator
    r.Gp_trial1[ic] = Gms[atomic_state1](ip1); // interpolate Gp for the other interval
    // Here we find out which hybridization is affected by this time change. We get index id_which, such that diag[id_which[ic]][isie] gives hybridization line affected
    auto [isie, id_whichi ] = WhichDelta(itc, diag_c[ic]);
    r.id_which[ic] = id_whichi;
    // hybridization line effected has spin and is,ie give below
    auto [dt, ie, is, spin, sgn] = WhatDelta(times(itc), beta, times, isie, dlttyp_c[ic], diag_c[ic], r.id_which[ic]);
    intpar ip;
    int ieis = ie*(2*max_order)+is;
    auto itt = r.ip_ieis.find(ieis);
    if (itt!=r.ip_ieis.end())
      ip = itt->second;
    else{
      ip = Interp(dt, tau);
      r.ip_ieis[ieis] = ip;
    }
    //cout<<"   At changing itc="<<itc<<" and ic="<<ic<<" need to change ie,is=("<<ie<<","<<is<<") dt="<<times(ie)<<"-"<<times(is)<<"="<<dt<<" "<<ip<<endl;
    r.Delt_trial[ic] = Deltas[spin](ip);
    // Now computing Phi from old Phi
    double dppp = (r.Gp_trial0[ic] * r.Gp_trial1[ic])/(Gp(itcm1,atomic_state0) * Gp(itc,atomic_state1));
    double dpdd = r.Delt_trial[ic] / Delt(spin, ie, is );
    r.dPhi[ic] = (dppp * dpdd);
    //if (Print){
    //  cout<<ic<<" Gm("<<atomic_state0<<",dt0="<<dt0<<")="<<r.Gp_trial0[ic]<<" Gm_old="<<Gp(itcm1,atomic_state0)<<endl;
    //  cout<<ic<<" Gm("<<atomic_state1<<",dt1="<<(times(itc+1)-times(itc))<<")="<<r.Gp_trial1[ic]<<" Gm_old="<<Gp(itc,atomic_state1)<<endl;
    //  cout<<ic<<" Dlt("<<spin<<","<<ie<<","<<is<<")="<<r.Delt_trial[ic]<<"  Dlt_old="<<Delt(spin, ie, is )<<endl;
    //  cout<<"dPhi["<<ic<<"]="<<r.dPhi[ic]<<endl;
    //}
  }
  return r;
}
auto GetTrialValues2(int itc, int Norder, double beta,
		     const TinyVector<Array<int,2>,4>& diag_c, const TinyVector<Array<int,1>,4>& pseud_c, const TinyVector<Array<int,1>,4>& dlttyp_c,
		     const Array<double,1>& times, const vector<Spline1D<double>>& Gms, const vector<Spline1D<double>>& Deltas,
		     const Array<double,1>& tau, const Array<double,2>& Gp, const Array<double,3>& Delt){
  struct retVals {
    TinyVector<double,4> dPhi;
    TinyVector<double,4> Gp_trial0;
    TinyVector<double,4> Gp_trial1;
    TinyVector<double,4> Delt_trial;
    TinyVector<int,4> id_which;
    map<int,intpar> ip_ieis;
  };
  retVals r;
  double dt0 = (itc>0) ? times(itc)-times(itc-1) : times(itc)-(times(2*Norder-1)-beta);
  intpar ip0 = Interp(dt0,tau);
  int itcm1 = (itc>0) ? itc-1 : 2*Norder-1;
  intpar ip1 = Interp(times(itc+1)-times(itc),tau);
  //intpar ip0 = Interp(t_new-t00, tau);            // trial_time(itc)-time(itc-1) but for general case
  //intpar ip1 = Interp(times(itc+1)-t_new, tau);   // time(itc+1)-trial_time(itc)
  for (int ic=0; ic<4; ++ic){
    int atomic_state0 = pseud_c[ic](itcm1); // itc-1 propagator
    r.Gp_trial0[ic] = Gms[atomic_state0](ip0); // interpolate Gp for new time difference
    int atomic_state1 = pseud_c[ic](itc);   // itc propagator
    r.Gp_trial1[ic] = Gms[atomic_state1](ip1); // interpolate Gp for the other interval
    // Here we find out which hybridization is affected by this time change. We get index id_which, such that diag[id_which[ic]][isie] gives hybridization line affected
    auto [isie, id_whichi ] = WhichDelta(itc, diag_c[ic]);
    r.id_which[ic] = id_whichi;
    // hybridization line effected has spin and is,ie give below
    auto [dt, ie, is, spin, sgn] = WhatDelta(times(itc), beta, times, isie, dlttyp_c[ic], diag_c[ic], r.id_which[ic]);
    intpar ip;
    int ieis = ie*(2*max_order)+is;
    auto itt = r.ip_ieis.find(ieis);
    if (itt!=r.ip_ieis.end())
      ip = itt->second;
    else{
      ip = Interp(dt, tau);
      r.ip_ieis[ieis] = ip;
    }
    //cout<<"   At changing itc="<<itc<<" and ic="<<ic<<" need to change ie,is=("<<ie<<","<<is<<") dt="<<times(ie)<<"-"<<times(is)<<"="<<dt<<" "<<ip<<endl;
    r.Delt_trial[ic] = Deltas[spin](ip);
    // Now computing Phi from old Phi
    double dppp = (r.Gp_trial0[ic] * r.Gp_trial1[ic])/(Gp(itcm1,atomic_state0) * Gp(itc,atomic_state1));
    double dpdd = r.Delt_trial[ic] / Delt(spin, ie, is );
    r.dPhi[ic] = (dppp * dpdd);
  }
  return r;
}
/*
auto GetTrialValues(double t_new, int itc, int itcm1, double t00, double beta,
		    const TinyVector<Array<int,2>,4>& diag_c, const TinyVector<Array<int,1>,4>& pseud_c, const TinyVector<Array<int,1>,4>& dlttyp_c,
		    const Array<double,1>& times, const vector<Spline1D<double>>& Gms, const vector<Spline1D<double>>& Deltas,
		    const Array<double,1>& tau, const Array<double,2>& Gp, const Array<double,3>& Delt){
  struct retVals {
    TinyVector<double,4> dPhi;
    TinyVector<double,4> Gp_trial0;
    TinyVector<double,4> Gp_trial1;
    TinyVector<double,4> Delt_trial;
    TinyVector<int,4> id_which;
    map<int,intpar> ip_ieis;
  };
  retVals r;
  intpar ip0 = Interp(t_new-t00, tau);            // trial_time(itc)-time(itc-1) but for general case
  intpar ip1 = Interp(times(itc+1)-t_new, tau);   // time(itc+1)-trial_time(itc)
  for (int ic=0; ic<4; ++ic){
    int atomic_state0 = pseud_c[ic](itcm1); // itc-1 propagator
    r.Gp_trial0[ic] = Gms[atomic_state0](ip0); // interpolate Gp for new time difference
    int atomic_state1 = pseud_c[ic](itc);   // itc propagator
    r.Gp_trial1[ic] = Gms[atomic_state1](ip1); // interpolate Gp for the other interval
    // Here we find out which hybridization is affected by this time change. We get index id_which, such that diag[id_which[ic]][isie] gives hybridization line affected
    auto [isie, id_whichi ] = WhichDelta(itc, diag_c[ic]);
    r.id_which[ic] = id_whichi;
    // hybridization line effected has spin and is,ie give below
    auto [dt, ie, is, spin, sgn] = WhatDelta(t_new, beta, times, isie, dlttyp_c[ic], diag_c[ic], r.id_which[ic]);
    intpar ip;
    int ieis = ie*(2*max_order)+is;
    auto itt = r.ip_ieis.find(ieis);
    if (itt!=r.ip_ieis.end())
      ip = itt->second;
    else{
      ip = Interp(dt, tau);
      r.ip_ieis[ieis] = ip;
    }
    r.Delt_trial[ic] = Deltas[spin](ip);
    // Now computing Phi from old Phi
    double dppp = (r.Gp_trial0[ic] * r.Gp_trial1[ic])/(Gp(itcm1,atomic_state0) * Gp(itc,atomic_state1));
    double dpdd = r.Delt_trial[ic] / Delt(spin, ie, is );
    r.dPhi[ic] = (dppp * dpdd);
  }
  return r;
}
*/
auto GetTrialValues0(double t_new, int itc, int itcm1, double t00, double beta,
		    const TinyVector<Array<int,2>,4>& diag_c, const TinyVector<Array<int,1>,4>& pseud_c, const TinyVector<Array<int,1>,4>& dlttyp_c,
		    const Array<double,1>& times, const vector<Spline1D<double>>& Gms, const vector<Spline1D<double>>& Deltas,
		    const Array<double,1>& tau, const Array<double,2>& Gp, const Array<double,3>& Delt){
  struct retVals {
    TinyVector<double,4> dPhi;
    TinyVector<double,4> Gp_trial0;
    TinyVector<double,4> Gp_trial1;
    TinyVector<double,4> Delt_trial;
    TinyVector<int,4> id_which;
  };
  retVals r;
  Linintpar ip0 = LInterp(t_new-t00, tau);            // trial_time(itc)-time(itc-1) but for general case
  Linintpar ip1 = LInterp(times(itc+1)-t_new, tau);   // time(itc+1)-trial_time(itc)
  for (int ic=0; ic<4; ++ic){
    int atomic_state0 = pseud_c[ic](itcm1); // itc-1 propagator
    r.Gp_trial0[ic] = Gms[atomic_state0](ip0); // interpolate Gp for new time difference
    int atomic_state1 = pseud_c[ic](itc);   // itc propagator
    r.Gp_trial1[ic] = Gms[atomic_state1](ip1); // interpolate Gp for the other interval
    // Here we find out which hybridization is affected by this time change. We get index id_which, such that diag[id_which[ic]][isie] gives hybridization line affected
    auto [isie, id_whichi ] = WhichDelta(itc, diag_c[ic]);
    r.id_which[ic] = id_whichi;
    // hybridization line effected has spin and is,ie give below
    auto [dt, ie, is, spin, sgn] = WhatDelta(t_new, beta, times, isie, dlttyp_c[ic], diag_c[ic], r.id_which[ic]);
    Linintpar ip = LInterp(dt, tau);
    r.Delt_trial[ic] = Deltas[spin](ip);
    // Now computing Phi from old Phi
    double dppp = (r.Gp_trial0[ic] * r.Gp_trial1[ic])/(Gp(itcm1,atomic_state0) * Gp(itc,atomic_state1));
    double dpdd = r.Delt_trial[ic] / Delt(spin, ie, is );
    r.dPhi[ic] = (dppp * dpdd);
  }
  return r;
}

