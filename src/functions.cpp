/**
 lassosum
 functions.cpp
 Purpose: functions to perform lassosum

 @author Timothy Mak
 @author Robert Porsch

 @version 0.1

 */
// [[Rcpp::interfaces(r, cpp)]]

#include <stdio.h>
#include <string>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/**
 Opens a Plink binary files

 @s file name
 @BIT ifstream
 @return is plink file in major mode

 */

bool openPlinkBinaryFile(const std::string s, std::ifstream &BIT) {
  BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  if (!BIT.is_open()) {
    throw "Cannot open the bed file";
  }

  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  char ch[1];
  BIT.read(ch, 1);
  std::bitset<8> b;
  b = ch[0];
  bool bfile_SNP_major = false;
  bool v1_bfile = true;
  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  // std::cerr << "check magic number" << std::endl;
  if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
    // Next number
    BIT.read(ch, 1);
    b = ch[0];
    if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
      // Read SNP/Ind major coding
      BIT.read(ch, 1);
      b = ch[0];
      if (b[0])
        bfile_SNP_major = true;
      else
        bfile_SNP_major = false;

      // if (bfile_SNP_major) std::cerr << "Detected that binary PED file is
      // v1.00 SNP-major mode" << std::endl;
      // else std::cerr << "Detected that binary PED file is v1.00
      // individual-major mode" << std::endl;

    } else
      v1_bfile = false;

  } else
    v1_bfile = false;
  // Reset file if < v1
  if (!v1_bfile) {
    Rcerr << "Warning, old BED file <v1.00 : will try to recover..."
          << std::endl;
    Rcerr << "  but you should --make-bed from PED )" << std::endl;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    BIT.read(ch, 1);
    b = ch[0];
  }
  // If 0.99 file format
  if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
    Rcerr << std::endl
          << " *** Possible problem: guessing that BED is < v0.99      *** "
          << std::endl;
    Rcerr << " *** High chance of data corruption, spurious results    *** "
          << std::endl;
    Rcerr
      << " *** Unless you are _sure_ this really is an old BED file *** "
      << std::endl;
    Rcerr << " *** you should recreate PED -> BED                      *** "
          << std::endl
          << std::endl;
    bfile_SNP_major = false;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  } else if (!v1_bfile) {
    if (b[0])
      bfile_SNP_major = true;
    else
      bfile_SNP_major = false;
    Rcerr << "Binary PED file is v0.99" << std::endl;
    if (bfile_SNP_major)
      Rcerr << "Detected that binary PED file is in SNP-major mode"
            << std::endl;
      else
        Rcerr << "Detected that binary PED file is in individual-major mode"
              << std::endl;
  }
  return bfile_SNP_major;
}

//' Count number of lines in a text file
//'
//' @param fileName Name of file
//' @keywords internal
//'
// [[Rcpp::export]]
int countlines(const char* fileName) {

  // Stolen from http://stackoverflow.com/questions/3482064/counting-the-number-of-lines-in-a-text-file
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(fileName);

  while (std::getline(myfile, line))
    ++number_of_lines;
  return number_of_lines;
}

//' Multiply genotypeMatrix by a matrix
//'
//' @param fileName location of bam file
//' @param N number of subjects
//' @param P number of positions
//' @param input the matrix
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @return an armadillo genotype matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat multiBed3(const std::string fileName, int N, int P, const arma::mat input,
                    arma::Col<int> col_skip_pos, arma::Col<int> col_skip,
                    arma::Col<int> keepbytes, arma::Col<int> keepoffset,
                    const int trace) {

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName, bedFile);
  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
                               "individual-major mode. Please use the snp-major "
                               "format");

  int i = 0;
  int ii = 0;
  int iii = 0;
  const bool colskip = (col_skip_pos.n_elem > 0);
  unsigned long long int Nbytes = ceil(N / 4.0);
  const bool selectrow = (keepbytes.n_elem > 0);
  int n;
  if (selectrow)
    n = keepbytes.n_elem;
  else
    n = N;
  int jj;

  arma::mat result = arma::mat(n, input.n_cols, arma::fill::zeros);
  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  int chunk;
  double step;
  double Step = 0;
  if(trace > 0) {
    chunk = input.n_rows / pow(10, trace);
    step = 100 / pow(10, trace);
    // Rcout << "Started C++ program \n";
  }

  while (i < P) {
    Rcpp::checkUserInterrupt();
    if (colskip) {
      if (ii < col_skip.n_elem) {
        if (i == col_skip_pos[ii]) {
          bedFile.seekg(col_skip[ii] * Nbytes, bedFile.cur);
          i = i + col_skip[ii];
          ii++;
          continue;
        }
      }
    }

    if(trace > 0) {
      if (iii % chunk == 0) {
        Rcout << Step << "% done\n";
        Step = Step + step;
      }
    }

    bedFile.read(ch, Nbytes); // Read the information
    if (!bedFile)
      throw std::runtime_error(
          "Problem with the BED file...has the FAM/BIM file been changed?");

    int j = 0;
    if (!selectrow) {
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];

        int c = 0;
        while (c < 7 && j < N) { // from the original PLINK: 7 because of 8 bits
          int first = b[c++];
          int second = b[c++];
          if (first == 0) {
            for (int k = 0; k < input.n_cols; k++) {
              if (input(iii, k) != 0.0) {
                result(j, k) += (2 - second) * input(iii, k);
              }
            }
          }
          j++;
        }
      }
    } else {
      for (jj = 0; jj < keepbytes.n_elem; jj++) {
        b = ch[keepbytes[jj]];

        int c = keepoffset[jj];
        int first = b[c++];
        int second = b[c];
        if (first == 0) {
          for (int k = 0; k < input.n_cols; k++) {
            if (input(iii, k) != 0.0) {
              result(j, k) += (2 - second) * input(iii, k);
            }
          }
        }
        j++;
      }
    }

    i++;
    iii++;
  }

  return result;
}


//' Multiply genotypeMatrix by a matrix (sparse)
//'
//' @param fileName location of bam file
//' @param N number of subjects
//' @param P number of positions
//' @param input the matrix
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @return an armadillo genotype matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat multiBed3sp(const std::string fileName, int N, int P,
                      const arma::vec beta,
                      const arma::Col<int> nonzeros,
                      const arma::Col<int> colpos,
                      const int ncol,
                      arma::Col<int> col_skip_pos, arma::Col<int> col_skip,
                      arma::Col<int> keepbytes, arma::Col<int> keepoffset,
                      const int trace) {

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName, bedFile);
  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
                               "individual-major mode. Please use the snp-major "
                               "format");

  int i = 0;
  int ii = 0;
  int iii = 0;
  int k = 0;
  const bool colskip = (col_skip_pos.n_elem > 0);
  unsigned long long int Nbytes = ceil(N / 4.0);
  const bool selectrow = (keepbytes.n_elem > 0);
  int n;
  if (selectrow)
    n = keepbytes.n_elem;
  else
    n = N;
  int jj;

  arma::mat result = arma::mat(n, ncol, arma::fill::zeros);
  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  int chunk;
  double step;
  double Step = 0;
  if(trace > 0) {
    chunk = nonzeros.n_elem / pow(10, trace);
    step = 100 / pow(10, trace);
    // Rcout << "Started C++ program \n";
  }

  while (i < P) {
    Rcpp::checkUserInterrupt();
    if (colskip) {
      if (ii < col_skip.n_elem) {
        if (i == col_skip_pos[ii]) {
          bedFile.seekg(col_skip[ii] * Nbytes, bedFile.cur);
          i = i + col_skip[ii];
          ii++;
          continue;
        }
      }
    }

    if(trace > 0) {
      if (iii % chunk == 0) {
        Rcout << Step << "% done\n";
        Step = Step + step;
      }
    }

    bedFile.read(ch, Nbytes); // Read the information
    if (!bedFile)
      throw std::runtime_error(
          "Problem with the BED file...has the FAM/BIM file been changed?");

    int j = 0;
    if (!selectrow) {
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];

        int c = 0;
        while (c < 7 && j < N) { // from the original PLINK: 7 because of 8 bits
          int first = b[c++];
          int second = b[c++];
          if(nonzeros[iii] > 0) {
            if (first == 0) {
              for (int kk = 0; kk < nonzeros[iii]; kk++) {
                result(j, colpos[k]) += (2 - second) * beta[k];
                k++;
              }
              k -= nonzeros[iii];
            }
          }
          j++;
        }
      }
    } else {
      for (jj = 0; jj < keepbytes.n_elem; jj++) {
        b = ch[keepbytes[jj]];

        int c = keepoffset[jj];
        int first = b[c++];
        int second = b[c];
        if(nonzeros[iii] > 0) {
          if (first == 0) {
            for (int kk = 0; kk < nonzeros[iii]; kk++) {
              result(j, colpos[k]) += (2 - second) * beta[k];
              k++;
            }
            k -= nonzeros[iii];
          }
        }
        j++;
      }
    }

    k += nonzeros[iii];
    i++;
    iii++;
  }

  return result;
}



//' Performs elnet
//'
//' @param lambda1 lambda
//' @param lambda2 lambda
//' @param X genotype Matrix
//' @param r correlations
//' @param inv_Sb the inverse of the variance-covariance matrix of genetic effects
//' @param inv_Ss the inverse of the residual variance matrix
//' @param x beta coef
//' @param thr threshold
//' @param yhat a vector
//' @param trace if >1 displays the current iteration
//' @param maxiter maximal number of iterations
//' @param sample_size
//' @return conv
//' @keywords internal
//'
// [[Rcpp::export]]

int elnet(double lambda1, double lambda2, const arma::vec& diag, const arma::mat& X,
          const arma::vec& r, const arma ::mat& inv_Sb,const arma ::mat& inv_Ss ,double thr, arma::vec& x, arma::vec& yhat, int trace, int maxiter,
          const arma::vec& sample_size)
{


  // diag is basically diag(X'X)
  // Also, ensure that the yhat=X*x in the input. Usually, both x and yhat are preset to 0.
  // They are modified in place in this function.

  // we modify this since X is one phenotype genotype matrix now :
  //int nq =X.n_rows;
 // int pq =X.n_cols;
 // int q = inv_Sb.n_cols;
  //int p = (pq)/(q);
  
  int n = X.n_rows;
  int p = X.n_cols;
  int q = inv_Sb.n_cols;
  int pq = p*q;

  // avant je comparais r.n_elem avec pq, mais on me donnait warning, donc je compare directement
  // r.n_elem avec X.n_cols, pareil pour les autres
  if(r.n_elem != X.n_cols*q) stop("r.n_elem != X.n_cols*q");
  if(x.n_elem != X.n_cols*q) stop("x.n_elem != X.n_cols*q");
  if(yhat.n_elem != X.n_rows*q) stop("yhat.n_elem != X.n_rows*q");
  if(diag.n_elem != X.n_cols) stop("diag.n_elem != X.n_cols*q");

  double dlx,del,t1,t2,t3, A,S,C ;

  // j : indice des SNPs, k : indice des traits, m: indice des it?rations, u indice pour parcourir le vecteur x ( des betas ),
  // h indice utilis? pour d?finir t1 , l incide utilis? pour d?finir t3
  int j,k,m,u,h,l;

  // On d?finit le vecteur x_before ( qui contient les valeurs des betas ? l'it?ration t-1 )
  arma :: vec x_before(pq,arma::fill::zeros);

  // We modify too because one phenotype matrix genotype
  //arma::vec Lambda2(pq);
  arma::vec Lambda2(p);

  Lambda2.fill(lambda2);
  arma::vec denom=diag + Lambda2;

  int conv=0;

  for(m=0;m<maxiter ;m++) {
    dlx=0.0;
    // Mon beta est x : c'est un vecteur de taille pq : x = ( q betas pour le SNP1 , q betas pour le SNP 2, .., q betas pour le SNP p )

    x_before = x;

    // boucle sur les SNPS :
    for (j= 0; j<p; j++) {

      // Pour chaque SNP, on fait une boucle sur les traits :
      for(k=0; k < q; k++) {

        x.at(q*j+k)=0.0;


        // On initialise les termes dont on aura besoin : t1, t2 et t3 ( ce sont les 3 composantes de A comme d?finie dans la partie th?orique )
        // ainsi que A et S

        t1 = 0.0 ;
        t2 = 0.0 ;
        t3 = 0.0 ;
        A = 0.0 ;
        S = 0.0 ;
        C = 0.0 ;

        // RMQ : c++ commence ? indicer ? partir de 0 ( le premier ?l?ment d'un vecteur ? l'indice 0),
        // alors que R commence ? indicer ? partir de 1 ( le premier ?l?ment d'un vecteur ? l'indice 1 )

        // On d?finit le terme t1 :

        for (h =0 ; h<q;h++){
          if (h!=k) t1=t1+ inv_Sb.at(k,h)*x.at(q*j+h);
        }

        // Rmq quand j'?crivais (-1/2)*t1, on me donnais t1=0 car il consid?re
        // 1/2 comme 0 ( divison enti?re )
        t1=-(0.5)*t1;

        // On d?finit le terme t2 :

        t2 = sample_size.at(k)*(arma::dot(inv_Ss.row(k),r.subvec(q*j,q*(j+1)-1)))  ;

        // On d?finit le terme t3 :
       
        // Vu qu'on va travailler avec une matrice one genotype, on modifie ceci 
        // for(l=0; l< p;l++){
         // C = trans(X.cols(q*j,q*(j+1)-1))*X.cols(q*l,q*(l+1)-1)*x.subvec(q*l,q*(l+1)-1);
         // if(l!=j) S = S + C.at(k);
        //}
        
        for(l=0; l< p;l++){
          C = arma::dot(X.col(j), X.col(l))*x.at(q*l+k);
          if(l!=j) S = S + C;
        }

        t3 = -sample_size.at(k)*(inv_Ss.at(k,k))*S;

        A=t1+t2+t3;

        // On d?finit maintenant la solution Beta

        if (A < 0){
          if (A + lambda1 <=0 ) {
            x.at(q*j+k) = (A+ lambda1)/(inv_Ss.at(k,k)*sample_size.at(k)*denom(j)+inv_Sb.at(k,k));
          }
        }

        if (A > 0){
          if (A - lambda1>= 0){
            x.at(q*j+k) = (A- lambda1)/(inv_Ss.at(k,k)*sample_size.at(k)*denom(j)+inv_Sb.at(k,k));
          }
        }

        if (x.at(q*j+k)==x_before.at(q*j+k) ) continue;
        del = x.at(q*j+k)-x_before.at(q*j+k);
        dlx=std::max(dlx,std::abs(del));
        
        // On remplace ceci vu qu'on travaille avec matrice one phenotype X 
        //yhat += del*X.col(q*j+k);
        
        arma::vec SS(X.n_rows*inv_Sb.n_cols,arma::fill::zeros);
        int d = (q*j+k)%q;
        for (int b=0;b< X.n_rows;b++){
          SS.at(q*b+d) = X(b,j)*del;
          }
        yhat += SS;
      }
    }

    checkUserInterrupt();
    if(trace > 0) Rcout << "Iteration: " << m << "\n";

    if(dlx < thr) {
      conv=1;
      break;
    }
  }

  return conv;
}


//' Performs elnet with s = 1
//'
//' @param lambda1 lambda
//' @param r correlations
//' @param inv_Sb the inverse of the variance-covariance matrix of genetic effects
//' @param inv_Ss the inverse of the residual variance matrix
//' @param x beta coef
//' @param thr threshold
//' @param yhat a vector
//' @param trace if >1 displays the current iteration
//' @param maxiter maximal number of iterations
//' @param sample_size
//' @return conv


// [[Rcpp::export]]

int elnet_s1(double lambda1,const arma::vec& r,int p, int q, int pq, const arma ::mat& inv_Sb,
          const arma ::mat& inv_Ss ,double thr, arma::vec& x, 
          int trace, int maxiter,const arma::vec& sample_size)
{
    


  double dlx,del,t1,t2,t3, A,S,C ;

  // j : indice des SNPs, k : indice des traits, m: indice des it?rations, u indice pour parcourir le vecteur x ( des betas ),
  // h indice utilis? pour d?finir t1 , l incide utilis? pour d?finir t3
  int j,k,m,u,h,l;

  // On d?finit le vecteur x_before ( qui contient les valeurs des betas ? l'it?ration t-1 )
  arma :: vec x_before(pq,arma::fill::zeros);

  int conv=0;

  for(m=0;m<maxiter ;m++) {
    dlx=0.0;
    // Mon beta est x : c'est un vecteur de taille pq : x = ( q betas pour le SNP1 , q betas pour le SNP 2, .., q betas pour le SNP p )

    x_before = x;

    // boucle sur les SNPS :
    for (j= 0; j<p; j++) {

      // Pour chaque SNP, on fait une boucle sur les traits :
      for(k=0; k < q; k++) {

        x.at(q*j+k)=0.0;


        // On initialise les termes dont on aura besoin : t1, t2 et t3 ( ce sont les 3 composantes de A comme d?finie dans la partie th?orique )
        // ainsi que A et S

        t1 = 0.0 ;
        t2 = 0.0 ;
        A = 0.0 ;

        // On d?finit le terme t1 :

        for (h =0 ; h<q;h++){
          if (h!=k) t1=t1+ inv_Sb.at(k,h)*x.at(q*j+h);
        }

        t1=-(0.5)*t1;

        // On d?finit le terme t2 :

        t2 = sample_size.at(k)*(arma::dot(inv_Ss.row(k),r.subvec(q*j,q*(j+1)-1)))  ;

        A=t1+t2;

        // On d?finit maintenant la solution Beta

        if (A < 0){
          if (A + lambda1 <=0 ) {
            x.at(q*j+k) = (A+ lambda1)/(sample_size.at(k)*inv_Ss.at(k,k)+inv_Sb.at(k,k));
          }
        }

        if (A > 0){
          if (A - lambda1>= 0){
            x.at(q*j+k) = (A- lambda1)/(inv_Ss.at(k,k)+inv_Sb.at(k,k));
          }
        }

        if (x.at(q*j+k)==x_before.at(q*j+k) ) continue;
        del = x.at(q*j+k)-x_before.at(q*j+k);
        dlx=std::max(dlx,std::abs(del));
      }
    }

    checkUserInterrupt();
    if(trace > 0) Rcout << "Iteration: " << m << "\n";

    if(dlx < thr) {
      conv=1;
      break;
    }
  }
 return conv;
}

// [[Rcpp::export]]
int repelnet(double lambda1, double lambda2, arma::vec& diag, arma::mat& X, arma::vec& r, arma ::mat& inv_Sb, arma ::mat& inv_Ss,
             double thr, arma::vec& x, arma::vec& yhat, int trace, int maxiter, const arma::vec& sample_size,
             arma::Col<int>& startvec, arma::Col<int>& endvec)
{
  int q = inv_Sb.n_cols;
  // Repeatedly call elnet by blocks...
  int nreps=startvec.n_elem;
  int out=1;
  for(int i=0;i < startvec.n_elem; i++) {
    arma::vec xtouse=x.subvec(startvec(i)*q, endvec(i)*q+(q-1));
    // this was valid when we had multiple phenotype matrix 
    //arma::vec yhattouse=X.cols(startvec(i)*q, endvec(i)*q+(q-1)) * xtouse;
    // We do this instead now : 
    
    arma::vec yhattouse(X.n_rows*inv_Sb.n_cols,arma::fill::zeros);
    int k = 0; 
    for(int j=startvec(i);j<endvec(i)+1;j++){
      arma::vec S(X.n_rows*inv_Sb.n_cols,arma::fill::zeros);
      for (int l=0;l< X.n_rows;l++){
        arma::mat M = arma::mat(q,q);
        M.fill(X(l,j));
        M = diagmat(M);
        S.subvec(l*q,l*q+(q-1)) = M*xtouse.subvec(k*q,k*q+(q-1));
      }
      k += 1;
      yhattouse+= S;
    }
    int out2=elnet(lambda1, lambda2,
                   diag.subvec(startvec(i), endvec(i)),
                   // Vu qu'on travaille avec la matrice one genotype, on modifie
                   // ceci : 
                   //X.cols(startvec(i)*q, endvec(i)*q+(q-1)),
                   X.cols(startvec(i), endvec(i)),
                   r.subvec(startvec(i)*q, endvec(i)*q+(q-1)),
                   inv_Sb,inv_Ss,
                   thr, xtouse,
                   yhattouse, trace - 1, maxiter,sample_size);
    x.subvec(startvec(i)*q, endvec(i)*q+(q-1))=xtouse;
    yhat += yhattouse;
    if(trace > 0) Rcout << "Block: " << i << "\n";
    out=std::min(out, out2);
  }
  return out;
}

//' imports genotypeMatrix
//'
//' @param fileName location of bam file
//' @param N number of subjects
//' @param P number of positions
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @return an armadillo genotype matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat genotypeMatrix(const std::string fileName, int N, int P,
                         arma::Col<int> col_skip_pos, arma::Col<int> col_skip,
                         arma::Col<int> keepbytes, arma::Col<int> keepoffset,
                         const int fillmissing) {

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName, bedFile);

  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
                               "individual-major mode. Please use the snp-major "
                               "format");

  int i = 0;
  int ii = 0;
  const bool colskip = (col_skip_pos.n_elem > 0);
  unsigned long long int Nbytes = ceil(N / 4.0);
  const bool selectrow = (keepbytes.n_elem > 0);
  int n, p, nskip;
  if (selectrow)
    n = keepbytes.n_elem;
  else
    n = N;

  if (colskip) {
    nskip = arma::accu(col_skip);
    p = P - nskip;
  }  else
    p = P;

  int j, jj, iii;

  arma::mat genotypes = arma::mat(n, p, arma::fill::zeros);
  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  iii=0;
  while (i < P) {
    // Rcout << i << std::endl;
    Rcpp::checkUserInterrupt();
    if (colskip) {
      if (ii < col_skip.n_elem) {
        if (i == col_skip_pos[ii]) {
          bedFile.seekg(col_skip[ii] * Nbytes, bedFile.cur);
          i = i + col_skip[ii];
          ii++;
          continue;
        }
      }
    }

    bedFile.read(ch, Nbytes); // Read the information
    if (!bedFile)
      throw std::runtime_error(
          "Problem with the BED file...has the FAM/BIM file been changed?");

    j = 0;
    if (!selectrow) {
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];

        int c = 0;
        while (c < 7 &&
               j < N) { // from the original PLINK: 7 because of 8 bits
          int first = b[c++];
          int second = b[c++];
          if (first == 0) {
            genotypes(j, iii) = (2 - second);
          }
          if(fillmissing == 0 && first == 1 && second == 0) genotypes(j, iii) =arma::datum::nan;
          j++;
        }
      }
    } else {
      for (jj = 0; jj < keepbytes.n_elem; jj++) {
        b = ch[keepbytes[jj]];

        int c = keepoffset[jj];
        int first = b[c++];
        int second = b[c];
        if (first == 0) {
          genotypes(j, iii) = (2 - second);
        }
        if(fillmissing == 0 && first == 1 && second == 0) genotypes(j, iii) =arma::datum::nan;
        j++;
      }
    }
    i++;
    iii++;
  }
  return genotypes;
}


//' normalize genotype matrix
//'
//' @param genotypes a armadillo genotype matrix
//' @return standard deviation
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec normalize(arma::mat &genotypes)
{
  int k = genotypes.n_cols;
  int n = genotypes.n_rows;
  arma::vec sd(k);
  for (int i = 0; i < k; ++i) {
    double m = arma::mean(genotypes.col(i));
    arma::vec mm(n); mm.fill(m);
    sd(i) = arma::stddev(genotypes.col(i));
    // sd(i) = 1.0;
    genotypes.col(i) = arma::normalise(genotypes.col(i) - mm);
  }
  return sd;
}

//' We build a function that gives us the correlation matrix of SNPs
//'
//' Correlation matrix of genotype matrix
//'
//' @param genotypes a armadillo genotype matrix
//' @return correlation matrix
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat Correlation(arma::mat &genotypes)
{
  int k = genotypes.n_cols;
  arma::mat Correlation = arma::mat(k, k, arma::fill::zeros);
  Correlation = cor(genotypes);
  return Correlation;
}


//' Runs elnet with various parameters
//'
//' @param lambda1 a vector of lambdas (lambda2 is 0)
//' @param fileName the file name of the reference panel
//' @param cor a matrix of correlations, rows represent phenotypes, and columns represent SNPs
//' @param inv_Sb the inverse of the variance-covariance matrix of genetic effects
//' @param inv_Ss the inverse of the residual variance matrix
//' @param N number of subjects
//' @param P number of position in reference file
//' @param col_skip_posR which variants should we skip
//' @param col_skipR which variants should we skip
//' @param keepbytesR required to read the PLINK file
//' @param keepoffsetR required to read the PLINK file
//' @param thr threshold
//' @param init a numeric matrix of beta coefficients
//' @param trace if >1 displays the current iteration
//' @param maxiter maximal number of iterations
//' @param Constant a constant to multiply the standardized genotype matrix
//' @return a list of results
//' @keywords internal
//'

// [[Rcpp::export]]



List runElnet(arma::vec& lambda, double shrink, const std::string fileName,
              arma::mat& cor, arma ::mat& inv_Sb ,arma ::mat& inv_Ss,int N, int P,
              arma::Col<int>& col_skip_pos, arma::Col<int>& col_skip,
              arma::Col<int>& keepbytes, arma::Col<int>& keepoffset,
              double thr, arma::mat& init, int trace, int maxiter,const arma::vec& sample_size,
              arma::Col<int>& startvec, arma::Col<int>& endvec) {
  // a) read bed file
  // b) standardize genotype matrix
  // c) multiply by constatant factor
  // d) perfrom elnet

  // Rcout << "ABC" << std::endl;

  int i,j;

  arma::mat genotypes = genotypeMatrix(fileName, N, P, col_skip_pos, col_skip, keepbytes,
                                                     keepoffset, 1);

  // On commence par normaliser la matrice g?notype pour un seul ph?notype

  arma::vec sd = normalize(genotypes);

  genotypes *= sqrt(1.0 - shrink);

  // Rcout << "DEF" << std::endl;

  // Ici, je transforme cor et init en des vecteurs pour pouvoir travailler avec :

  arma::vec r(cor.n_rows*cor.n_cols,arma::fill::zeros);
  int b = 0;
  for(int a=0; a < r.n_elem; a+=cor.n_rows) {
    r.subvec(a,a+cor.n_rows-1) = cor.col(b);
    b=b+1;
  }

  arma::vec x(init.n_rows*init.n_cols,arma::fill::zeros);
  int q = 0;
  for(int e=0; e < x.n_elem; e+=init.n_rows) {
    x.subvec(e,e+init.n_rows-1) = init.col(q);
    q=q+1;
  }

  if (genotypes.n_cols != cor.n_cols) {
    throw std::runtime_error("Number of positions in reference file is not "
                               "equal the number of regression coefficients");
  }

  arma::Col<int> conv(lambda.n_elem);
  int len = r.n_elem;

  arma::mat beta(len, lambda.n_elem);
  arma::mat pred(genotypes.n_rows*inv_Sb.n_cols, lambda.n_elem); pred.zeros();
  arma::vec out(lambda.n_elem);
  arma::vec loss(lambda.n_elem);
  // On definit diag avec le nombre de SNP car on ne travaille pas avec matrice
  // multiple phenotype
  //arma::vec diag(r.n_elem); diag.fill(1.0 - shrink);
  arma::vec diag(genotypes.n_cols); diag.fill(1.0 - shrink);
  // Rcout << "HIJ" << std::endl;

  for(j=0; j < diag.n_elem; j++) {
    if(sd(j) == 0.0)
        {
       diag(j) = 0.0;
        // j'ajoute cette condition pour le traitement des SNPs monomorphiques
        // selon Mak, pour les SNPs monomorphiques, on suppose que Xi ( i etant 
        // l'indice du SNP i) est un vecteur de rempli de 0. 
       genotypes.col(j).fill(0);
       }
  }
  // Rcout << "LMN" << std::endl;

  arma::vec fbeta(lambda.n_elem);
  // On initialise le vecteur yhat avec des 0, sinon on obtient des valeurs incorrects
  // de y_hat, et donc de pred
  arma::vec yhat(genotypes.n_rows*inv_Sb.n_cols,arma::fill::zeros);
  //yhat = genotypes * x;


  // Rcout << "Starting loop" << std::endl;
  for (i = 0; i < lambda.n_elem; ++i) {
    if (trace > 0)
      Rcout << "lambda: " << lambda(i) << "\n" << std::endl;
    out(i) =
      repelnet(lambda(i), shrink, diag,genotypes, r,inv_Sb,inv_Ss, thr, x, yhat, trace-1, maxiter,sample_size,
               startvec, endvec);
    beta.col(i) = x;
    for(j=0; j < beta.n_rows; j++) {
      int k = j/inv_Sb.n_cols;
      if(sd(k) == 0.0) beta(j,i)=beta(j,i) * shrink;
    }

    if (out(i) != 1) {
      throw std::runtime_error("Not converging.....");
    }

    pred.col(i) = yhat;

  }

  return List::create(Named("lambda") = lambda,
                      Named("beta") = beta,
                      Named("conv") = out,
                      Named("pred") = pred,
                      Named("loss") = loss,
                      Named("fbeta") = fbeta,
                      Named("sd")= sd);
}

//' Runs elnet_s1 with various parameters
//'
//' @param lambda1 a vector of lambdas 
//' @param cor a matrix of correlations, rows represent phenotypes, and columns represent SNPs
//' @param inv_Sb the inverse of the variance-covariance matrix of genetic effects
//' @param inv_Ss the inverse of the residual variance matrix
//' @param thr threshold
//' @param init a numeric matrix of beta coefficients
//' @param trace if >1 displays the current iteration
//' @param maxiter maximal number of iterations
//' @return a list of results
//' @keywords internal
//'

// [[Rcpp::export]]



List runElnet_s1(arma::vec& lambda,arma::mat& cor, arma ::mat& inv_Sb ,arma ::mat& inv_Ss,
              double thr, arma::mat& init, int trace, int maxiter,const arma::vec& sample_size) {


    int i,nbr_snps,nbr_trait,nbr_snps_trait ;
    
    nbr_snps = cor.n_cols ;
    nbr_trait = cor.n_rows;
    nbr_snps_trait = nbr_snps*nbr_trait;

 
  // Ici, je transforme cor et init en des vecteurs pour pouvoir travailler avec :

  arma::vec r(cor.n_rows*cor.n_cols,arma::fill::zeros);
  int b = 0;
  for(int a=0; a < r.n_elem; a+=cor.n_rows) {
    r.subvec(a,a+cor.n_rows-1) = cor.col(b);
    b=b+1;
  }

  arma::vec x(init.n_rows*init.n_cols,arma::fill::zeros);
  int q = 0;
  for(int e=0; e < x.n_elem; e+=init.n_rows) {
    x.subvec(e,e+init.n_rows-1) = init.col(q);
    q=q+1;
  }

  arma::Col<int> conv(lambda.n_elem);
  int len = r.n_elem;

  arma::mat beta(len, lambda.n_elem);
  arma::vec out(lambda.n_elem);

  // Rcout << "Starting loop" << std::endl;
  for (i = 0; i < lambda.n_elem; ++i) {
    if (trace > 0)
      Rcout << "lambda: " << lambda(i) << "\n" << std::endl;
    out.at(i) =
      elnet_s1(lambda(i),r,nbr_snps,nbr_trait,nbr_snps_trait,inv_Sb,inv_Ss, thr, x, trace-1, maxiter,sample_size);
    beta.col(i) = x;

    if (out(i) != 1) {
      throw std::runtime_error("Not converging.....");
    }
  }

  return List::create(Named("lambda") = lambda,
                      Named("beta") = beta,
                      Named("conv") = out);
}
