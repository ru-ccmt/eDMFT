// @Copyright 2007 Kristjan Haule
// 

double Segment_ComputeTryalExponentSum(function2D<NState>& state_evolution_left, int op_i, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, int dsize, long long istep)
{
  const int VACUUM = 0; // when superstate are indexed from 1, vacuum state is 0
  int nsize = Operators.size()+dsize;
  
  static function1D<double> exp_sum(praStates.size());
  
  double smallest = -std::numeric_limits<double>::max();
  for (int ist=0; ist<praStates.size(); ist++) exp_sum[ist]=smallest;
    
  double max_val=smallest;
  for (int ist=0; ist<praStates.size(); ist++){
    if (op_i>0 && (state_evolution_left[ist].size()<op_i)) continue;// this product matrix element is zero
    
    int istate = praStates[ist].istate;  // istate index superstates from 1 (vacuum = 0)
    bool survived=true;
    double survived_exp_sum=0.0;
    for (int ip=0; ip<nsize; ++ip) {
      //survived_exp_sum += Operators.tr_exp_(ip)(istate,0); // for segments, there is only one exponent
      survived_exp_sum += Operators.tr_exp_(ip)[istate-1]; // for segments, there is only one exponent

      //cout<<"exp for "<<istate<<" = "<<Operators.tr_exp_(ip)[istate-1]<<endl;
      
      int iop = Operators.tr_typ(ip);
      istate = cluster.Fi(iop, istate);
      if (istate == VACUUM){
        survived=false;
        break;
      }
    }
    if (survived && istate==praStates[ist].istate){
      //survived_exp_sum += Operators.tr_exp_last()(istate,0);
      survived_exp_sum += Operators.tr_exp_last()[istate-1];
      exp_sum[ist] = survived_exp_sum;
      if (survived_exp_sum>max_val) max_val=survived_exp_sum;
    }
  }
  double dval=0.0;
  for (int ist=0; ist<praStates.size(); ist++){
    if (exp_sum[ist]-max_val > -300) dval += exp(exp_sum[ist]-max_val);
  }
  return max_val + log(dval);
}


Number Segment_UpdateStateEvolution(function2D<NState>& state_evolution_left, int op_i, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster,
				    function1D<Number>& Trace, function2D<double>& Prob, long long istep)
{
  int npra = praStates.size();
  int nsize = Operators.size();
  for (int ist=0; ist<npra; ist++) Prob(ist,0) = 0.0;
  
  Number ms=0;
  for (int ist=0; ist<npra; ist++){
    Trace[ist]=0;
    if (op_i>0 && (state_evolution_left[ist].size()<op_i)) continue;// this product matrix element is zero
    int istate = praStates[ist].istate;  // istate index superstates from 1 (vacuum = 0)
    bool survived=true;
    double mantisa=1;
    double survived_exp_sum=0.0;
    int ip=0;
    for (; ip<nsize; ++ip){                     // go over the all operators
      //survived_exp_sum += Operators.exp_(ip)(istate,0); // This is like Evolve
      survived_exp_sum += Operators.exp_(ip)[istate-1]; // This is like Evolve
      int iop = Operators.typ(ip);
      int new_istate = cluster.Fi(iop, istate);                 // This is the second part of apply
      state_evolution_left(ist,ip).istate=new_istate;
      if (new_istate == 0){
	survived=false;
	break;
      }
      mantisa *= cluster.FM(iop,istate)(0,0);          // This is like apply
      istate = new_istate;
    }
    if (survived && istate==praStates[ist].istate){
      //survived_exp_sum += Operators.exp_last()(istate,0); // This is like Evolve
      survived_exp_sum += Operators.exp_last()[istate-1]; // This is like Evolve
      Trace[ist] = Number(mantisa, survived_exp_sum);
      ms += Trace[ist];
      ip++;
    }
    state_evolution_left[ist].resize(ip); // we have stored state evolution, remember the new size
  }

  if (nsize==0){
    for (int ist=0; ist<npra; ist++) Prob(ist,0) = cluster.P_atom(ist,0);
  }else{
    for (int ist=0; ist<npra; ist++){
      double weight_of_chain = divide(Trace[ist],ms);
      if (abs(weight_of_chain)<1e-10) continue;
      for (int ip=0; ip<nsize-1; ++ip){                     // go over the all operators
	int istate = state_evolution_left(ist,ip).istate;
	Prob(istate-1,0) += weight_of_chain*(Operators.t(ip+1)-Operators.t(ip))/common::beta;
      }
      Prob(ist,0) +=  weight_of_chain*(1-(Operators.t(nsize-1)-Operators.t(0))/common::beta);
    }
  }
  return ms;
}


bool Segment_Try_Add_Cd_C_(int bfls, double& t_start, int bfle, double& t_end, const nIntervals& interval, long long istep)
{ // Only trial step // migh succeed or might not. This step should not modify the status of the system
  // This is very optimized routine which is called at every trial move when adding a kink.
  // Replaces the simple lookups in the table tries to determin if the matrix element is zero or not.
  if (bfls!=bfle) return false;
  if (interval.size()==0) return true;
  pair<double,double> ts = interval.Closest_Start_Times(bfls, t_start);
  if (ts.first==-common::beta && ts.second==-common::beta) return true;  // There is no kink of this type yet.
  pair<double,double> te = interval.Closest_End_Times(bfls, t_start);
  //if (t_start==t_end){
  //  t_end -= common::smallest_dt;
  //}
  if (t_end > t_start){
    if (te.first>ts.first && te.second>ts.second && t_end<ts.second) return true;
    if (te.first<ts.first && te.second<ts.second && t_end-common::beta>ts.first) return true;
    return false;
  }else{
    if (te.first<ts.first && te.second<ts.second && t_end>ts.first) return true;
    if (te.first>ts.first && te.second>ts.second && t_end+common::beta<ts.second) return true;
    return false;
  }
}

bool Segment_Try_Remove_C_Cd_(int bfle, int ie, int bfls, int is, const nIntervals& interval, long long istep)
{
  if (bfls!=bfle) return false;
  if (interval.size()==2) return true;

  double t_start = interval.time_s(is);
  double t_end = interval.time_e(ie);

  double ts0 = interval.PreviousTime(nIntervals::cd, is, true);
  double ts1 = interval.NextTime(nIntervals::cd, is, true);
  
  if (t_end > t_start){
    if (ts1 > t_end) return true;
    if (ts0 < t_end-common::beta) return true;
    return false;
  }else{
    if (ts0 < t_end) return true;
    if (ts1 > t_end+common::beta) return true;
    return false;
  }
}

bool Segment_Try_Move_(int bfl, int type, int to_move, double t_new, double t_old, double t_prev, double t_next, const nIntervals& interval, long long istep)
{
  pair<double,double> tn = (type==nIntervals::cd) ? interval.Closest_End_Times(bfl, t_old) : interval.Closest_Start_Times(bfl, t_old);
  if (tn.first < t_new && t_new < tn.second){
    return true;
  }else if(tn.first< t_new+common::beta && t_new+common::beta < tn.second){
    return true;
  }else if(tn.first< t_new-common::beta && t_new-common::beta < tn.second){
    return true;
  }else return false;
}

bool Segment_TryFlip(int iop_a, int iop_b, int op_a, double t_a, int op_b, double t_b, const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{
  /*
  const double small=1e-12;
  iop_a=-1, iop_b=-1;
  for (int ip=0; ip<Operators.size(); ++ip){  // go over all operators
    int op = Operators.typ(ip);
    double tc = Operators.t(ip);
    if (abs(tc-t_a)<small && op==op_a) iop_a=ip;
    if (abs(tc-t_b)<small && op==op_b) iop_b=ip;
  }
  if (iop_a<0 || iop_b<0) {cout<<istep<<" ERROR Could not determine position in Segment_TryFlip: iop_a="<<iop_a<<" iop_b="<<iop_b<<endl;}
  */
  
  int npra = praStates.size();
  int nsize = Operators.size();
  bool one_survived=false;
  for (int ist=0; ist<npra; ist++){
    int istate = praStates[ist].istate;  // istate index superstates from 1 (vacuum = 0)
    bool survived=true;
    
    int ip=0;
    for (; ip<nsize; ++ip){                     // go over the all operators
      //survived_exp_sum += Operators.exp_(ip)(istate,0); // This is like Evolve
      int op = Operators.typ(ip);
      double tc = Operators.t(ip);
      if (ip==iop_a) op=op_b;
      if (ip==iop_b) op=op_a;
      int new_istate = cluster.Fi(op, istate);                 // This is the second part of apply
      if (new_istate == 0){
	survived=false;
	break;
      }
      istate = new_istate;
    }
    if (survived && istate==praStates[ist].istate){
      one_survived=true;
    }
  }
  return one_survived;
}

bool ExponentSumForExhangeInside1(double& survived_exp_sum, int ist, int ip_a, int ip_b, const NOperators& Operators,
				  const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{
  const int VACUUM = 0; // when superstate are indexed from 1, vacuum state is 0
  int istate = praStates[ist].istate;  // istate index superstates from 1 (vacuum = 0)
  survived_exp_sum=0.0;
  for (int ip=0; ip<Operators.size(); ++ip){
    //survived_exp_sum += Operators.exp_(ip)(istate,0); // for segments, there is only one exponent
    survived_exp_sum += Operators.exp_(ip)[istate-1]; // for segments, there is only one exponent
    int iop = Operators.typ(ip);
    if (ip==ip_a) iop = Operators.typ(ip_b);
    if (ip==ip_b) iop = Operators.typ(ip_a);
    istate = cluster.Fi(iop, istate);
    if (istate == VACUUM) return false;
  }
  if (istate==praStates[ist].istate){
    //survived_exp_sum += Operators.exp_last()(istate,0);
    survived_exp_sum += Operators.exp_last()[istate-1];
    return true;
  }
  return false;
}

bool ExponentSumForExhangeInside2(double& exp_sum, bool& connects, int ist, int ip_1, int iop_1, int ip_2, int iop_2, const function1D<Number>& Trace,
				   const function2D<NState>& state_evolution_left,
				   const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{
  const int VACUUM = 0; // when superstate are indexed from 1, vacuum state is 0
  int istate = (ip_1!=0) ? state_evolution_left(ist,ip_1-1).istate : praStates[ist].istate;
  if (istate == VACUUM) return false;
  exp_sum=0.0;
  int iop = Operators.typ(ip_1);
  istate = cluster.Fi(iop_1, istate);
  if (istate == VACUUM) return false;
  for (int ip=ip_1+1; ip<ip_2; ++ip){
    //exp_sum += Operators.exp_(ip)(istate,0);
    exp_sum += Operators.exp_(ip)[istate-1];
    iop = Operators.typ(ip);
    istate = cluster.Fi(iop, istate);
    if (istate == VACUUM) return false;
  }
  //exp_sum += Operators.exp_(ip_2)(istate,0);
  exp_sum += Operators.exp_(ip_2)[istate-1];
  iop = Operators.typ(ip_2);
  istate = cluster.Fi(iop_2, istate);
  if (istate == VACUUM) return false;
  
  if (ip_1>0 && istate!=state_evolution_left(ist,ip_2).istate){
    connects=true;
    //cout<<istep<<" ERROR It should not happen ist="<<ist<<" istate="<<istate<<" istate'="<<state_evolution_left(ist,ip_2).istate<<" ip_1="<<ip_1<<" ip_2="<<ip_2<<" exp_sum="<<exp_sum<<endl;
  }else{ connects=false;}
  
  return true;
}




double Segment_ComputeTryalExponentSumForExchange(int ip_a, int ip_b, const function1D<Number>& Trace, const function2D<NState>& state_evolution_left,
						  const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{
  static function1D<double> exp_sum(praStates.size());
  double smallest = -std::numeric_limits<double>::max();
  for (int ist=0; ist<praStates.size(); ist++) exp_sum[ist]=smallest;

  int ip_1 = min(ip_a,ip_b);
  int ip_2 = max(ip_a,ip_b);
  double max_val=smallest;
  for (int ist=0; ist<praStates.size(); ist++){
    double final_exp_sum=smallest;
    bool can_do_it_fast=false;
    if (ip_1>0){
      if (state_evolution_left[ist].size()<ip_1) continue;// this product matrix element is zero
      double old_exp_sum, new_exp_sum;
      bool connects1, connects2;
      bool old_survived = ExponentSumForExhangeInside2(old_exp_sum, connects1, ist, ip_1, Operators.typ(ip_1), ip_2, Operators.typ(ip_2), Trace, state_evolution_left, Operators, praStates, cluster, istep);
      bool new_survived = ExponentSumForExhangeInside2(new_exp_sum, connects2, ist, ip_1, Operators.typ(ip_2), ip_2, Operators.typ(ip_1), Trace, state_evolution_left, Operators, praStates, cluster, istep);
      if (new_survived && old_survived){
	final_exp_sum = Trace[ist].exp_dbl()+new_exp_sum-old_exp_sum;
	can_do_it_fast=true;
      }else if (!new_survived && !old_survived){
	final_exp_sum = smallest;
	can_do_it_fast=true;
      };
    }
    if (ip_1==0 || !can_do_it_fast){
      double survived_exp_sum;
      bool survived = ExponentSumForExhangeInside1(survived_exp_sum,ist,ip_a,ip_b,Operators,praStates,cluster,istep);
      if (survived){
	exp_sum[ist] = survived_exp_sum;
	if (survived_exp_sum>max_val) max_val=survived_exp_sum;
      }else{
	exp_sum[ist] = smallest;
      }
    }else{
      exp_sum[ist] =  final_exp_sum;
      if (final_exp_sum>max_val) max_val=final_exp_sum;
    }
  }
  double dval=0.0;
  for (int ist=0; ist<praStates.size(); ist++){
    if (exp_sum[ist]-max_val > -300) dval += exp(exp_sum[ist]-max_val);
  }
  return max_val + log(dval);
}

void Segment_Get_N_at_Operator2(functionb<double>& Njc, int ip, const function1D<Number>& Trace, const Number& ms, const function2D<NState>& state_evolution_left,
				const NOperators& Operators, const function1D<NState>& praStates, const ClusterData& cluster, long long istep)
{
  Njc=0.0;
  for (int ist=0; ist<praStates.size(); ist++){
    double Zi = divide(Trace[ist],ms);
    if (fabs(Zi)>1e-10 && state_evolution_left[ist].size()>Operators.size()){
      int istate = (ip!=0) ? state_evolution_left[ist][ip-1].istate : praStates[ist].istate;
      for (int fl2=0; fl2<cluster.Nvfl; fl2++){
	switch (cluster.Njm_c(istate,fl2)){
	case 0: // nothing to add because N=0
	  break;
	case 1: // N is just an idenity. No need to multiply.
	  Njc[fl2] += Zi;
	  break;
	case 2: // N is not simple indentity or zero, hence we need to evaluate the product
	  break;
	  Njc[fl2] += Zi*cluster.Njm(istate,fl2)(0,0);
	}
      }
    }
  }
}
