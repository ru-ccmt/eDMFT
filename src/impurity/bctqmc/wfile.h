void PrintContributions(int j, const Array<double,2>& dSm, const Array<double,2>& dG00, const Array<double,2>& dF00, double Phi, const Array<double,1>& tau){
  //cout << "Phi(order="<<j<<")="<<Phi << endl;
  ofstream bout("dSig."s + to_string(j) + ".dat"s);
  for (int it=0; it<tau.size(); it++) bout<<tau(it)<<" "<<dSm(0,it)<<" "<<dSm(1,it)<<" "<<dSm(2,it)<<" "<<dSm(3,it)<<endl;
  ofstream pout("dG00."s + to_string(j) + ".dat"s);
  for (int it=0; it<tau.size(); it++) pout<<tau(it)<<" "<<dG00(0,it)<<" "<<dG00(1,it)<<endl;
  ofstream sout("dF00."s + to_string(j) + ".dat"s);
  for (int it=0; it<tau.size(); it++) sout<<tau(it)<<" "<<dF00(0,it)<<" "<<dF00(1,it)<<endl;
}
void PrintContributions(int j, const Array<double,2>& dSm, double Phi, const Array<double,1>& tau){
  //cout << "Phi(order="<<j<<")="<<Phi << endl;
  ofstream bout("dSig."s + to_string(j) + ".dat"s);
  for (int it=0; it<tau.size(); it++) bout<<tau(it)<<" "<<dSm(0,it)<<" "<<dSm(1,it)<<" "<<dSm(2,it)<<" "<<dSm(3,it)<<endl;
}
std::string ltrim(const std::string &s)
{
  const std::string WHITESPACE = " \n\r\t\f\v";
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}
std::string rtrim(const std::string &s)
{
  const std::string WHITESPACE = " \n\r\t\f\v";
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}
void RemoveComments(std::stringstream& inp2, std::ifstream& inp){
  //std::stringstream inp2;
  std::string line;
  while (inp){
    getline(inp, line);
    if (!inp) break;
    std::size_t found = line.find("#");
    if (found!=std::string::npos) line.erase(found);
    line = trim(line);
    if (not line.empty())
      inp2 << line << std::endl;
  }
  //return inp2;
}
void PrintContributions(int norder, const Array<double,2>& dSm, const Array<double,1>& tau){
  ofstream bout("dSig."s + to_string(norder) + ".dat"s);
  for (int it=0; it<tau.size(); it++) bout<<tau(it)<<" "<<dSm(0,it)<<" "<<dSm(1,it)<<" "<<dSm(2,it)<<" "<<dSm(3,it)<<endl;
}
void PrintContributions(int norder, const Array<double,2>& dG00, const Array<double,2>& dF00, const Array<double,1>& tau){
  ofstream pout("dG00."s + to_string(norder) + ".dat"s);
  for (int it=0; it<tau.size(); it++) pout<<tau(it)<<" "<<dG00(0,it)<<" "<<dG00(1,it)<<endl;
  ofstream sout("dF00."s + to_string(norder) + ".dat"s);
  for (int it=0; it<tau.size(); it++) sout<<tau(it)<<" "<<dF00(0,it)<<" "<<dF00(1,it)<<endl;
}
bool IsRealAxis(string fDelta){
  ifstream inpd(fDelta);
  double w;
  inpd>>w;
  //cout<<"from "<<fDelta<<" file w="<<w<<" and w<0="<< (w<0) << endl;
  return (w<0);
}
std::tuple<Array<double,1>,Array<double,2>> loadtxt(const string& filename, char delimiter = ' ')
{
  vector<vector<double>> tempData;
  ifstream file(filename);
  if (!file) {
    cerr << "Error: Unable to open file " << filename << endl;
    exit(1);
  }
  string line;
  int ii=0;
  while (getline(file, line)) {
    istringstream ss(line);
    vector<double> row;
    double value;
    while (ss >> value) {
      row.push_back(value);
      if (ss.peek() == delimiter) ss.ignore(); // Handle delimiters (comma/space/tab)
    }
    if (row.size()>0){
      tempData.push_back(row);
      //cout<<"row["<<ii<<"].size="<<row.size()<<endl;
      ii+=1;
    }
  }
  file.close();
  //cout<<"Done"<<endl;
  cout.flush();
  // Convert vector to Blitz++ array
  int rows = tempData.size();
  int cols = tempData[0].size();
  //cout<<"rows="<<rows<<" cols="<<cols<<endl;
  Array<double,1> x(rows);
  Array<double,2> ys(cols-1,rows);

  
  for (int i = 0; i < rows; ++i)  x(i) = tempData[i][0];
  for (int i = 0; i < rows; ++i)
    for (int j = 1; j < cols; ++j)
      ys(j-1, i) = tempData[i][j];
  
  //cout<<"From file "<<filename<<" read rows="<<rows<<" columns= "<<cols << endl;
  //cout<<"x.size="<<x.size()<<" ys.shape="<<ys.shape()<<endl;
  return {x,ys};
}
