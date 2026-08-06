// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Score a node's parent set (given an ALREADY-BUILT Cholesky factor) plus all
// its single-node "+1" extensions, via the bordered-Cholesky vectorized
// approach -- the C++ analog of bge_score_plus_from_factor / the second half
// of bge_score_plus_parent's lp>=1 branch. All node ids here are 0-indexed
// (into TN), unlike the R-facing entry point below.
// choltemp: upper-triangular Cholesky factor of TN(parentnodes0,parentnodes0), lp x lp
// c_noplus: solve(choltemp^T, TN(j0,parentnodes0)), length lp
// logdetD:  log|TN(parentnodes0,parentnodes0)|, consistent with choltemp
static arma::rowvec score_from_factor(int j0, const arma::uvec& parentnodes0,
                                       const arma::uvec& plus0, const arma::mat& TN,
                                       const arma::mat& choltemp, const arma::vec& c_noplus,
                                       double logdetD, double A, int lp,
                                       double awpN, int N, const arma::vec& scoreconstvec){
  int lpp = plus0.n_elem;
  double awpNd2 = (awpN - N + lp + 1) / 2.0;

  arma::vec TN_j_plus(lpp);
  for(int b=0; b<lpp; b++) TN_j_plus(b) = TN(j0, plus0(b));

  arma::mat TN_plus_off(lp, lpp);
  for(int a=0; a<lp; a++)
    for(int b=0; b<lpp; b++)
      TN_plus_off(a,b) = TN(parentnodes0(a), plus0(b));

  arma::vec TN_plus_diag(lpp);
  for(int b=0; b<lpp; b++) TN_plus_diag(b) = TN(plus0(b), plus0(b));

  double c_noplus_sq = arma::as_scalar(arma::sum(arma::square(c_noplus)));
  double logdetpart2 = std::log(A - c_noplus_sq);
  // scoreconstvec here is 0-indexed; R's scoreconstvec[lp+1] (1-indexed) == scoreconstvec(lp) here
  double corescore = scoreconstvec(lp) - awpNd2*logdetpart2 - logdetD/2.0;

  arma::mat choltemp_new_12 = arma::solve(arma::trimatl(choltemp.t()), TN_plus_off); // lp x lpp

  arma::vec choltemp_new_22(lpp);
  for(int b=0; b<lpp; b++){
    double s = arma::as_scalar(arma::sum(arma::square(choltemp_new_12.col(b))));
    choltemp_new_22(b) = std::sqrt(TN_plus_diag(b) - s);
  }

  arma::vec c_plusses(lpp);
  for(int b=0; b<lpp; b++){
    double dotpart = arma::as_scalar(arma::dot(choltemp_new_12.col(b), c_noplus));
    c_plusses(b) = (TN_j_plus(b) - dotpart) / choltemp_new_22(b);
  }

  arma::rowvec out(lpp+1);
  out(0) = corescore;
  for(int b=0; b<lpp; b++){
    double logdetD_plus = logdetD + 2*std::log(choltemp_new_22(b));
    double logdetpart2_plus = std::log(A - c_noplus_sq - c_plusses(b)*c_plusses(b));
    out(b+1) = scoreconstvec(lp+1) - (awpNd2+0.5)*logdetpart2_plus - logdetD_plus/2.0;
  }
  return out;
}

// The lp==0 leaf case (empty parent set), mirroring bge_score_plus_parent's
// own "0" branch -- no Cholesky factor needed since D is diagonal-only.
static arma::rowvec score_lp0(int j0, const arma::uvec& plus0, const arma::mat& TN,
                               double A, double awpN, int N, const arma::vec& scoreconstvec){
  int lpp = plus0.n_elem;
  double awpNd2 = (awpN - N + 1) / 2.0;
  arma::rowvec out(lpp+1);
  out(0) = scoreconstvec(0) - awpNd2*std::log(A);
  for(int b=0; b<lpp; b++){
    double D = TN(plus0(b), plus0(b));
    double Bv = TN(j0, plus0(b));
    double logdetpart2 = std::log(A - Bv*Bv/D);
    out(b+1) = scoreconstvec(1) - (awpNd2+0.5)*logdetpart2 - std::log(D)/2.0;
  }
  return out;
}

// Depth-first recursion over the K-candidate inclusion tree, building the
// Cholesky factor incrementally via bordered/rank-1 updates. "Excluding" a
// candidate is free (just recurse with the same factor); the C++ call stack
// restores the parent's factor automatically on backtracking, so no
// downdate is ever computed.
static void recurse(int depth, int K, std::vector<int>& included_local,
                     const arma::mat& choltemp, const arma::vec& c_noplus, double logdetD,
                     int j0, const arma::uvec& cand0, const arma::uvec& plus0,
                     const arma::mat& TN, double A, double awpN, int N,
                     const arma::vec& scoreconstvec, const IntegerVector& maps_backwards,
                     arma::mat& score_matr){
  if(depth > K){
    int lp = included_local.size();
    long long bitmask = 0;
    for(int idx : included_local) bitmask += (1LL << (idx-1));
    int row_idx = maps_backwards[bitmask] - 1;  // maps_backwards is 1-indexed (R row numbers)

    arma::rowvec row_scores;
    if(lp==0){
      row_scores = score_lp0(j0, plus0, TN, A, awpN, N, scoreconstvec);
    } else {
      arma::uvec parentnodes0(lp);
      for(int a=0; a<lp; a++) parentnodes0(a) = cand0(included_local[a]-1);
      row_scores = score_from_factor(j0, parentnodes0, plus0, TN, choltemp, c_noplus,
                                     logdetD, A, lp, awpN, N, scoreconstvec);
    }
    score_matr.row(row_idx) = row_scores;
    return;
  }

  // exclude candidate `depth`
  recurse(depth+1, K, included_local, choltemp, c_noplus, logdetD,
          j0, cand0, plus0, TN, A, awpN, N, scoreconstvec, maps_backwards, score_matr);

  // include candidate `depth`: bordered Cholesky update, then recurse
  int new_node0 = cand0(depth-1);
  int lp_old = included_local.size();
  double d_new = TN(new_node0, new_node0);

  arma::mat new_choltemp;
  arma::vec new_c_noplus;
  double new_logdetD;

  if(lp_old==0){
    double u = std::sqrt(d_new);
    new_choltemp = arma::mat(1,1); new_choltemp(0,0) = u;
    new_c_noplus = arma::vec(1); new_c_noplus(0) = TN(j0, new_node0) / u;
    new_logdetD = 2*std::log(u);
  } else {
    arma::vec v(lp_old);
    for(int a=0; a<lp_old; a++) v(a) = TN(cand0(included_local[a]-1), new_node0);
    arma::vec w = arma::solve(arma::trimatl(choltemp.t()), v);
    double u = std::sqrt(d_new - arma::as_scalar(arma::sum(arma::square(w))));

    new_choltemp = arma::mat(lp_old+1, lp_old+1, arma::fill::zeros);
    new_choltemp.submat(0,0,lp_old-1,lp_old-1) = choltemp;
    new_choltemp.submat(0,lp_old,lp_old-1,lp_old) = w;
    new_choltemp(lp_old,lp_old) = u;

    double b_new = TN(j0, new_node0);
    double new_c_last = (b_new - arma::as_scalar(arma::dot(w, c_noplus))) / u;
    new_c_noplus = arma::vec(lp_old+1);
    new_c_noplus.subvec(0,lp_old-1) = c_noplus;
    new_c_noplus(lp_old) = new_c_last;

    new_logdetD = logdetD + 2*std::log(u);
  }

  included_local.push_back(depth);
  recurse(depth+1, K, included_local, new_choltemp, new_c_noplus, new_logdetD,
          j0, cand0, plus0, TN, A, awpN, N, scoreconstvec, maps_backwards, score_matr);
  included_local.pop_back();
}

//' Build the full base-score + all-plus-1-extension-score table for node j,
//' across all 2^K subsets of its K candidate parents, via depth-first
//' Cholesky-factor sharing (see build_plus_score_table_dfs in
//' BROOD_Functions.R for the R/documentation reference this mirrors exactly).
//' All R-facing arguments are 1-indexed, matching R's own convention;
//' converted to 0-indexed internally.
//' @param j                scoring node (1-indexed)
//' @param cand_nodes        the K candidate parent node ids (1-indexed)
//' @param plus_parentnodes  the "+1" candidate node ids (1-indexed)
//' @param maps_backwards    map_pars$maps[[i]]$backwards (1-indexed row numbers)
//' @param TN                the full TN matrix
//' @param awpN              param$awpN
//' @param scoreconstvec     param$scoreconstvec
// [[Rcpp::export]]
NumericMatrix build_plus_score_table_cpp(int j, IntegerVector cand_nodes, IntegerVector plus_parentnodes,
                                         IntegerVector maps_backwards, NumericMatrix TN_r,
                                         double awpN, NumericVector scoreconstvec_r){
  int K = cand_nodes.size();
  int lpp = plus_parentnodes.size();
  int N = TN_r.nrow();
  long long n_parent_sets = 1LL << K;

  arma::mat TN = as<arma::mat>(TN_r);
  arma::vec scoreconstvec = as<arma::vec>(scoreconstvec_r);
  int j0 = j-1;
  arma::uvec cand0(K);
  for(int k=0; k<K; k++) cand0(k) = cand_nodes[k]-1;
  arma::uvec plus0(lpp);
  for(int k=0; k<lpp; k++) plus0(k) = plus_parentnodes[k]-1;
  double A = TN(j0,j0);

  arma::mat score_matr(n_parent_sets, lpp+1);
  std::vector<int> included_local;
  included_local.reserve(K);

  recurse(1, K, included_local, arma::mat(), arma::vec(), 0.0,
          j0, cand0, plus0, TN, A, awpN, N, scoreconstvec, maps_backwards, score_matr);

  return wrap(score_matr);
}

// ============================================================================
// DAG-Wishart analog of the above, mirroring dagwishart_score_plus_parent().
// The DAG-Wishart score combines TWO parallel terms per subset -- a
// posterior-side term (matrix UN) and a prior-side term (matrix U0) -- added
// together. Both terms use the SAME candidate-parent subset at every point,
// so both Cholesky factors are built along the SAME recursive traversal.
//
// One structural asymmetry, confirmed directly against the R source and
// replicated exactly here (not "corrected"): for the base (lp-parent) score,
// BOTH the UN-side exponent (awpNd2_new) and U0-side exponent (awpd2_new)
// use the base lp. For the +1-extension score, ONLY awpNd2_new gets the
// +1/2 shift; awpd2_new stays fixed at its base-lp value. This matches
// dagwishart_score_plus_parent's own "1" and "otherwise" branches exactly.
//
// logedgepvec (if present) is a plain vector indexed by parent-set SIZE
// (not by specific edge, unlike BGe's logedgepmat), and its "+1" contribution
// is a single scalar applied uniformly to all lpp candidates -- a simpler
// penalty structure than BGe's, matched here accordingly.
// ============================================================================

static arma::rowvec dagwishart_score_from_factors(int j0, const arma::uvec& parentnodes0,
                                                   const arma::uvec& plus0, const arma::mat& UN, const arma::mat& U0,
                                                   const arma::mat& cholUN, const arma::vec& cUN, double logdetUN,
                                                   const arma::mat& cholU0, const arma::vec& cU0, double logdetU0,
                                                   double A, double A0, int lp,
                                                   double awpN_new, int N, double scoreconst_lp,
                                                   double scoreconst_lpplus1,
                                                   bool has_penalty, double pen_base, double pen_plus){
  int lpp = plus0.n_elem;
  double awpNd2_new = (awpN_new - lp)/2.0 - 1.0;
  double awpd2_new  = (awpN_new - N - lp)/2.0 - 1.0;   // NOT shifted for the +1 extension

  arma::vec UN_j_plus(lpp), U0_j_plus(lpp);
  for(int b=0; b<lpp; b++){ UN_j_plus(b) = UN(j0, plus0(b)); U0_j_plus(b) = U0(j0, plus0(b)); }

  arma::mat UN_plus_off(lp, lpp), U0_plus_off(lp, lpp);
  for(int a=0; a<lp; a++)
    for(int b=0; b<lpp; b++){
      UN_plus_off(a,b) = UN(parentnodes0(a), plus0(b));
      U0_plus_off(a,b) = U0(parentnodes0(a), plus0(b));
    }

  arma::vec UN_plus_diag(lpp), U0_plus_diag(lpp);
  for(int b=0; b<lpp; b++){ UN_plus_diag(b) = UN(plus0(b), plus0(b)); U0_plus_diag(b) = U0(plus0(b), plus0(b)); }

  double cUN_sq = arma::as_scalar(arma::sum(arma::square(cUN)));
  double cU0_sq = arma::as_scalar(arma::sum(arma::square(cU0)));
  double logdetpart2_UN = std::log(A - cUN_sq);
  double logdetpart2_U0 = std::log(A0 - cU0_sq);

  double corescore = scoreconst_lp - awpNd2_new*logdetpart2_UN - logdetUN/2.0
                                        + awpd2_new*logdetpart2_U0 + logdetU0/2.0;

  arma::mat cholUN_12 = arma::solve(arma::trimatl(cholUN.t()), UN_plus_off);
  arma::mat cholU0_12 = arma::solve(arma::trimatl(cholU0.t()), U0_plus_off);

  arma::rowvec out(lpp+1);
  out(0) = corescore + (has_penalty ? pen_base : 0.0);

  for(int b=0; b<lpp; b++){
    double sUN = arma::as_scalar(arma::sum(arma::square(cholUN_12.col(b))));
    double sU0 = arma::as_scalar(arma::sum(arma::square(cholU0_12.col(b))));
    double uUN = std::sqrt(UN_plus_diag(b) - sUN);
    double uU0 = std::sqrt(U0_plus_diag(b) - sU0);

    double dotUN = arma::as_scalar(arma::dot(cholUN_12.col(b), cUN));
    double dotU0 = arma::as_scalar(arma::dot(cholU0_12.col(b), cU0));
    double cplusUN = (UN_j_plus(b) - dotUN) / uUN;
    double cplusU0 = (U0_j_plus(b) - dotU0) / uU0;

    double logdetUN_plus = logdetUN + 2*std::log(uUN);
    double logdetU0_plus = logdetU0 + 2*std::log(uU0);
    double logdetpart2_UN_plus = std::log(A - cUN_sq - cplusUN*cplusUN);
    double logdetpart2_U0_plus = std::log(A0 - cU0_sq - cplusU0*cplusU0);

    double val = scoreconst_lpplus1 - (awpNd2_new+0.5)*logdetpart2_UN_plus - logdetUN_plus/2.0
                                        + awpd2_new*logdetpart2_U0_plus + logdetU0_plus/2.0;
    out(b+1) = val + (has_penalty ? pen_plus : 0.0);
  }
  return out;
}

static arma::rowvec dagwishart_score_lp0(int j0, const arma::uvec& plus0, const arma::mat& UN, const arma::mat& U0,
                                         double A, double A0, double awpN_new, int N,
                                         double scoreconst_1, double scoreconst_2,
                                         bool has_penalty, double pen_base, double pen_plus){
  int lpp = plus0.n_elem;
  double awpNd2_new = (awpN_new - 0)/2.0 - 1.0;
  double awpd2_new  = (awpN_new - N - 0)/2.0 - 1.0;

  arma::rowvec out(lpp+1);
  out(0) = scoreconst_1 - awpNd2_new*std::log(A) + awpd2_new*std::log(A0) + (has_penalty ? pen_base : 0.0);

  for(int b=0; b<lpp; b++){
    double DUN = UN(plus0(b), plus0(b));
    double DU0 = U0(plus0(b), plus0(b));
    double BUN = UN(j0, plus0(b));
    double BU0 = U0(j0, plus0(b));
    double logdetpart2_UN = std::log(A - BUN*BUN/DUN);
    double logdetpart2_U0 = std::log(A0 - BU0*BU0/DU0);
    double val = scoreconst_2 - (awpNd2_new+0.5)*logdetpart2_UN - std::log(DUN)/2.0
                                  + awpd2_new*logdetpart2_U0 + std::log(DU0)/2.0;
    out(b+1) = val + (has_penalty ? pen_plus : 0.0);
  }
  return out;
}

static void dagwishart_recurse(int depth, int K, std::vector<int>& included_local,
                                const arma::mat& cholUN, const arma::vec& cUN, double logdetUN,
                                const arma::mat& cholU0, const arma::vec& cU0, double logdetU0,
                                int j0, const arma::uvec& cand0, const arma::uvec& plus0,
                                const arma::mat& UN, const arma::mat& U0, double A, double A0,
                                double awpN_new, int N, const List& scoreconstlist,
                                bool has_penalty, const arma::vec& logedgepvec,
                                const IntegerVector& maps_backwards, arma::mat& score_matr){
  if(depth > K){
    int lp = included_local.size();
    long long bitmask = 0;
    for(int idx : included_local) bitmask += (1LL << (idx-1));
    int row_idx = maps_backwards[bitmask] - 1;

    arma::rowvec row_scores;
    if(lp==0){
      double sc1 = as<NumericVector>(scoreconstlist[0])[j0];
      double sc2 = as<NumericVector>(scoreconstlist[1])[j0];
      // R's "0" branch applies NO penalty to the base score, and uses
      // logedgepvec[lp+1] (not [lp+2]) for the +1-extension penalty --
      // this differs from the "1"/"otherwise" branches' convention and is
      // replicated exactly here, not "corrected".
      double pen_base_lp0 = 0.0;
      double pen_plus_lp0 = has_penalty ? logedgepvec(lp) : 0.0;
      row_scores = dagwishart_score_lp0(j0, plus0, UN, U0, A, A0, awpN_new, N, sc1, sc2,
                                        has_penalty, pen_base_lp0, pen_plus_lp0);
    } else {
      double pen_base = has_penalty ? logedgepvec(lp) : 0.0;       // R: logedgepvec[lp+1], 0-indexed -> lp
      double pen_plus = has_penalty ? logedgepvec(lp+1) : 0.0;     // R: logedgepvec[lp+2], 0-indexed -> lp+1
      arma::uvec parentnodes0(lp);
      for(int a=0; a<lp; a++) parentnodes0(a) = cand0(included_local[a]-1);
      double sc_lp = as<NumericVector>(scoreconstlist[lp])[j0];
      double sc_lpplus1 = as<NumericVector>(scoreconstlist[lp+1])[j0];
      row_scores = dagwishart_score_from_factors(j0, parentnodes0, plus0, UN, U0,
                                                 cholUN, cUN, logdetUN, cholU0, cU0, logdetU0,
                                                 A, A0, lp, awpN_new, N, sc_lp, sc_lpplus1,
                                                 has_penalty, pen_base, pen_plus);
    }
    score_matr.row(row_idx) = row_scores;
    return;
  }

  dagwishart_recurse(depth+1, K, included_local, cholUN, cUN, logdetUN, cholU0, cU0, logdetU0,
                     j0, cand0, plus0, UN, U0, A, A0, awpN_new, N, scoreconstlist,
                     has_penalty, logedgepvec, maps_backwards, score_matr);

  int new_node0 = cand0(depth-1);
  int lp_old = included_local.size();
  double dUN = UN(new_node0, new_node0);
  double dU0 = U0(new_node0, new_node0);

  arma::mat new_cholUN, new_cholU0;
  arma::vec new_cUN, new_cU0;
  double new_logdetUN, new_logdetU0;

  if(lp_old==0){
    double uUN = std::sqrt(dUN);
    double uU0 = std::sqrt(dU0);
    new_cholUN = arma::mat(1,1); new_cholUN(0,0) = uUN;
    new_cholU0 = arma::mat(1,1); new_cholU0(0,0) = uU0;
    new_cUN = arma::vec(1); new_cUN(0) = UN(j0, new_node0) / uUN;
    new_cU0 = arma::vec(1); new_cU0(0) = U0(j0, new_node0) / uU0;
    new_logdetUN = 2*std::log(uUN);
    new_logdetU0 = 2*std::log(uU0);
  } else {
    arma::vec vUN(lp_old), vU0(lp_old);
    for(int a=0; a<lp_old; a++){
      int node0 = cand0(included_local[a]-1);
      vUN(a) = UN(node0, new_node0);
      vU0(a) = U0(node0, new_node0);
    }
    arma::vec wUN = arma::solve(arma::trimatl(cholUN.t()), vUN);
    arma::vec wU0 = arma::solve(arma::trimatl(cholU0.t()), vU0);
    double uUN = std::sqrt(dUN - arma::as_scalar(arma::sum(arma::square(wUN))));
    double uU0 = std::sqrt(dU0 - arma::as_scalar(arma::sum(arma::square(wU0))));

    new_cholUN = arma::mat(lp_old+1, lp_old+1, arma::fill::zeros);
    new_cholUN.submat(0,0,lp_old-1,lp_old-1) = cholUN;
    new_cholUN.submat(0,lp_old,lp_old-1,lp_old) = wUN;
    new_cholUN(lp_old,lp_old) = uUN;

    new_cholU0 = arma::mat(lp_old+1, lp_old+1, arma::fill::zeros);
    new_cholU0.submat(0,0,lp_old-1,lp_old-1) = cholU0;
    new_cholU0.submat(0,lp_old,lp_old-1,lp_old) = wU0;
    new_cholU0(lp_old,lp_old) = uU0;

    double bUN_new = UN(j0, new_node0);
    double bU0_new = U0(j0, new_node0);
    double new_cUN_last = (bUN_new - arma::as_scalar(arma::dot(wUN, cUN))) / uUN;
    double new_cU0_last = (bU0_new - arma::as_scalar(arma::dot(wU0, cU0))) / uU0;
    new_cUN = arma::vec(lp_old+1); new_cUN.subvec(0,lp_old-1) = cUN; new_cUN(lp_old) = new_cUN_last;
    new_cU0 = arma::vec(lp_old+1); new_cU0.subvec(0,lp_old-1) = cU0; new_cU0(lp_old) = new_cU0_last;

    new_logdetUN = logdetUN + 2*std::log(uUN);
    new_logdetU0 = logdetU0 + 2*std::log(uU0);
  }

  included_local.push_back(depth);
  dagwishart_recurse(depth+1, K, included_local, new_cholUN, new_cUN, new_logdetUN,
                     new_cholU0, new_cU0, new_logdetU0,
                     j0, cand0, plus0, UN, U0, A, A0, awpN_new, N, scoreconstlist,
                     has_penalty, logedgepvec, maps_backwards, score_matr);
  included_local.pop_back();
}

//' DAG-Wishart analog of build_plus_score_table_cpp -- see that function's
//' documentation for the general design (depth-first Cholesky sharing).
//' Here, TWO Cholesky factors (for UN and U0) are built in parallel along
//' the same traversal, since both matrices are indexed by the identical
//' candidate-parent subset at every point.
//' @param awpN_new         param$alpha_post[j] (node-specific, scalar for this call)
//' @param scoreconstlist   param$scoreconstlist, an R list of per-size, per-node vectors
//' @param logedgepvec      param$logedgepvec if present, else pass a zero-length vector
// [[Rcpp::export]]
NumericMatrix build_plus_score_table_dagwishart_cpp(int j, IntegerVector cand_nodes, IntegerVector plus_parentnodes,
                                                     IntegerVector maps_backwards, NumericMatrix UN_r, NumericMatrix U0_r,
                                                     double awpN_new, List scoreconstlist, NumericVector logedgepvec_r){
  int K = cand_nodes.size();
  int lpp = plus_parentnodes.size();
  int N = UN_r.nrow();
  long long n_parent_sets = 1LL << K;

  arma::mat UN = as<arma::mat>(UN_r);
  arma::mat U0 = as<arma::mat>(U0_r);
  int j0 = j-1;
  arma::uvec cand0(K);
  for(int k=0; k<K; k++) cand0(k) = cand_nodes[k]-1;
  arma::uvec plus0(lpp);
  for(int k=0; k<lpp; k++) plus0(k) = plus_parentnodes[k]-1;
  double A = UN(j0,j0);
  double A0 = U0(j0,j0);

  bool has_penalty = logedgepvec_r.size() > 0;
  arma::vec logedgepvec = has_penalty ? as<arma::vec>(logedgepvec_r) : arma::vec();

  arma::mat score_matr(n_parent_sets, lpp+1);
  std::vector<int> included_local;
  included_local.reserve(K);

  dagwishart_recurse(1, K, included_local, arma::mat(), arma::vec(), 0.0,
                     arma::mat(), arma::vec(), 0.0,
                     j0, cand0, plus0, UN, U0, A, A0, awpN_new, N, scoreconstlist,
                     has_penalty, logedgepvec, maps_backwards, score_matr);

  return wrap(score_matr);
}
