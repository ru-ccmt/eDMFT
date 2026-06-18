int SignSort(const Array<int,2>& diag)
{
  int Norder = diag.extent(0);
  Array<int,1> list(2*Norder);
  int N=2*Norder;
  //for (int i=0; i<2*Norder; ++i) list(i) = diag.data()[i];
  std::copy_n(diag.data(), Norder*2, list.data());
  for (int i=0; i<Norder; ++i)
    if (list(2*i)>list(2*i+1)) swap( list(2*i),list(2*i+1));
				 
  // bubble sort, which determines if the list is an even or odd permutation of ordered list
  int sign = 1;
  // Sorts l in place and returns it.
  for (int passesLeft=N-1; passesLeft>=0; --passesLeft){
    for (int index=0; index < passesLeft; index++){
      if ( list(index) > list(index + 1) ){
          int t = list(index);
          list(index) = list(index+1);
          list(index+1) = t;
          sign *= -1;
        }
    }
  }
  return sign;
}

TinyVector<int,2> SignMatrixElements(const Array<int,2>& cdiag, const Array<int,1>& cpseudo, const Array<int,4>& Dlt)
{
  int two_n = cdiag.extent(0)*2;
  int sgn_mm = 1;
  int sgn_bw = 1;
  for (int j=0; j<cdiag.extent(0); ++j){
    int a = cdiag(j,0), b = cdiag(j,1);
    if (a>b) swap(a,b);
    TinyVector<int,4> ii( (a-1+two_n)%two_n, a, (b-1+two_n)%two_n, b );
    TinyVector<int,4> mm( cpseudo(ii[0]), cpseudo(ii[1]), cpseudo(ii[2]), cpseudo(ii[3]) );
    int dlt = Dlt(mm[0],mm[1],mm[2],mm[3]);
    sgn_mm *= sign(dlt);
    if (abs(dlt) == 1 or abs(dlt) == 3) sgn_bw *= -1; // It turns out that 1&3 signal backward and 2&4 signal forward propagators.
  }
  return TinyVector<int,2>(sgn_mm, sgn_bw);
}
int SgnMElements(const Array<int,2>& cdiag, const Array<int,1>& cpseudo, const Array<int,4>& Dlt)
{
  int two_n = cdiag.extent(0)*2;
  int sgn_mm = 1;
  for (int j=0; j<cdiag.extent(0); ++j){
    int a = cdiag(j,0), b = cdiag(j,1);
    if (a>b) swap(a,b);
    TinyVector<int,4> ii( (a-1+two_n)%two_n, a, (b-1+two_n)%two_n, b );
    TinyVector<int,4> mm( cpseudo(ii[0]), cpseudo(ii[1]), cpseudo(ii[2]), cpseudo(ii[3]) );
    int dlt = Dlt(mm[0],mm[1],mm[2],mm[3]);
    sgn_mm *= sign(dlt);
  }
  return sgn_mm;
}
