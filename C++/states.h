// state space dimensions
const int NG = 2;    // number of genders
const int AG = 17;   // number of age groups
const int DS = 8;    // number of disease stages (CD4 stages)
const int TS = 4;    // number of ART treatment stages

class states {

 public:
  double X[NG][AG][DS][TS];

  int hiv_idx;
  int art_idx;

  /* operators */
  states operator+(const states& a);
  states operator*(const double& c);

  states & operator+=(const states& a);

};


    
