#include "./numerics.hpp"

void Numerics::gauleg(const double &x1, const double &x2, vec &x, vec &w, const int &n) {
  double z1, z, pp, p3, p2, p1;
  int m = (n + 1)/2;
  double xm = 0.5*(x2+x1);
  double xl = 0.5*(x2-x1);

  for(int i = 0; i < m; i++) {
    z = cos(M_PI*(i + 0.75)/(n + 0.5));

    do {
      p1 = 1.;
      p2 = 0.;
      for(int j = 0; j < n; j++){
        p3 = p2;
        p2 = p1;
        p1 = ((2.*j + 1.)*z*p2 - j*p3)/(j+1);
      }
      pp = n*(z*p1 - p2)/(z*z - 1.);
      z1 = z;
      z = z1 - p1/pp;
    } while(abs(z - z1) > 1e-15);
    x[i] = xm - xl*z;
    x[n-1-i] = xm+xl*z;
    w[i] = 2.*xl/((1. - z*z)*pp*pp);
    w[n-1-i] = w[i];
  }
}

complex<double> Numerics::gaussIntegration(function<complex<double>(complex<double>)> func) {
  std::vector<double> gaulegX(1201);
  std::vector<double> gaulegW(1201);

  gauleg(-M_PI, M_PI, gaulegX, gaulegW, 1201);

  complex<double> I(0., 1.);
  complex<double> sum(0.0, 0.0);
  for( int i = 0; i < 1201; i++) {
    complex<double> x = -exp(I*gaulegX[i]);
    sum += func(x)*gaulegW[i];
  }
  return sum;
}


// from https://gist.github.com/lilac/2464434
bool Numerics::invertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  matrix<double> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<double>(A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  // testInvMat(input, inverse);

  return true;
}


// method described in: http://math.uww.edu/~mcfarlat/inverse.htm
bool Numerics::invMat(const ublas::matrix<double> &mat, ublas::matrix<double> &invMat) {
  int n = mat.size1();

  // Augment input matrix with an identity matrix
  ublas::matrix<double> augmat(n, 2*n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<2*n; j++) {
      if( j<n) {
        augmat(i, j) = mat(i, j);
      } else if((i+n) == j) {
        augmat(i, j) = 1.0;
      } else {
        augmat(i, j) = 0.0;
      }
    }
  }

  // loop through pivots
  for(int p=0; p<n; p++) {
    double pEl = augmat(p, p); // pivot element;
    // normalize pivot to 1 ( divide whole row through pivot element )
    for(int col=0; col<2*n; col++) {
      augmat(p, col) = augmat(p, col)/pEl;
    }
    // loop through pivot-row
    for(int row=0; row<n; row++) {
      // set every row element ( except pivot element ) to 0
      // by subtracting a multiple of the pivot row
      if(row != p) {
        double facEl = augmat(row, p);
        for(int col=0; col<2*n; col++) {
          augmat(row, col) = augmat(row, col) - facEl*augmat(p, col);
        }
      }
    }
  }

  // get inverse Matrix from right side of augmented matrix
  for(int row=0; row<n; row++) {
    for(int col=0; col<n; col++) {
      invMat(row, col) = augmat(row, col+n);
    }
  }

  testInvMat(mat, invMat);


  return true;
}

void Numerics::testInvMat(
  const ublas::matrix<double> &mat,
  const ublas::matrix<double> &invMat
) {
  const int n = mat.size1();
  ublas::matrix<double> onemat(n, n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      for(int k=0; k<n; k++) {
        onemat(i,j) += invMat(i,k)*mat(k,j);
      }
    }
  }

  // namespace karma = boost::spirit::karma;
  // using namespace karma;
  // std::cout << "inverse Matrix test:" << std::endl;
  // std::cout << format_delimited(columns(onemat.size2()) [auto_], '\t', onemat.data()) << std::endl;
};

const std::vector<double> Numerics::zeta_ = {
  0,
  0,
  0,
  1.2020569031595942,
  0,
  1.036927755143369926,
  0,
  1.008349277381922827
};
