#include <vector>
using namespace std;

template<class container>
double romberg(const container& ff, double dh){
  int m = ff.size();
  int N=1;
  int ic = m-1;
  while( ic>1){
    ic = ic>>1;
    N++;
  }
  int m2 = (1<<(N-1)) + 1;
  if (m2 != m){
    cout<<"ERROR : romberg should have number of points 2**k+1"<<endl;
  }
  double h[N+1], r[N+1][N+1];
  for (int i = 1; i < N + 1; ++i) {
    h[i] = dh / static_cast<int>(pow(2.0,i-1));
  }
  r[1][1] = h[1] / 2 * (ff[0] + ff[m-1]);
  for (int i = 2; i < N + 1; ++i) {
    double coeff = 0;
    int dr = static_cast<int>(pow(2.0,N-i));
    for (int k = 1; k <= static_cast<int>(pow(2.0, i-2)); ++k) coeff += ff[(2*k-1)*dr];
    r[i][1] = 0.5 * (r[i-1][1] + h[i-1]*coeff);
  }
  for (int i = 2; i < N + 1; ++i) {
    for (int j = 2; j <= i; ++j) {
      r[i][j] = r[i][j - 1] + (r[i][j-1] - r[i-1][j-1])/(static_cast<int>(pow(4., j-1))-1);
    }
  }
  return r[N][N];
}

