#include "./numerics.hpp"

// from https://gist.github.com/lilac/2464434
static bool invertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<T> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

static bool Numerics::invMat(const ublas::matrix<double> &matrix, ublas::matrix<double> &inverse) {
  bool flag = true;
  int n = matrix.size1();
  double m;

  // Augment input matrix with an identity matrix
  ublas::matrix<T> augmatrix(n, 2*n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<2*n; j++) {
      if( j < n) {
        augmatrix(i, j) = matrix(i, j);
      } else if((i+n) == j) {
        augmatrix(i, j) = 1.0;
      } else {
        augmatrix(i, j) = 0.0;
      }
    }
  }

  // Reduce augmented matrix to upper trangular form
  for( int k=0; k<n-1; k++) {
    if(augmatrix(k, k) == 0.0) {
      flag = false;
      for( int i=k; i<n; i++) {
        if( augmatrix(i, k) != 0.0) {
          for(int j=0; j<2*n; j++) {
            augmatrix(k, j) += augmatrix(i, j);
          }
          flag = true;
          break;
        }
        if (!flag) {
          std::cout << "1. Matrix is non - invertible!" << std::endl;
          return false;
        }
      }
    }
    for( int j=k; j<n; j++) {
      m = augmatrix(j, k)/ augmatrix(k, k);
      for( int i=k; i<2*n; i++) {
        augmatrix(j, i) -= m*augmatrix(k,i);
      }
    }
  }
  // std::cout << "mat" << std::endl << std::setprecision(15);
  // for(int i = 0; i < augmatrix.size1(); i++) {
  //   for(int j = 0; j < augmatrix.size2(); j++) {
  //     std::cout << augmatrix(i,j) << "\t";
  //   }
  //   std::cout << std::endl;
  // }

  // Test for invertibility
  for( int i=0; i < n; i++) {
    if( augmatrix(i, i) == 0.0) {
      std::cout << "Matrix is non - invertible!" << std::endl;
      return false;
    }
  }
  std::cout << "here4" << std::endl;

  // Make diagonal elements as 1
  for( int i=0; i<n; i++) {
    m = augmatrix(i, i);
    for( int j = i; j<2*n; j++) {
      augmatrix(i, j) /= m;
    }
  }
  std::cout << "here5" << std::endl;

  // Reduce right side half of augmented matrix to identity matrix
  for( int k=n-2; k>1; k--) {
    for( int i=0; i<k; i++) {
      m = augmatrix(i, k+1);
      for( int j=k; j<2*n; j++) {
        augmatrix(i, j) -= augmatrix(k+1, j)*m;
      }
    }
  }
  std::cout << "here6" << std::endl;

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      inverse(i, j) = augmatrix(i, j+n);
    }
  }
  std::cout << "here7" << std::endl;

  return true;
}
