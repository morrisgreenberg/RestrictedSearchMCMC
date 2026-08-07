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

    double val = scoreconst_lpplus1 - (awpNd2_new-0.5)*logdetpart2_UN_plus - logdetUN_plus/2.0
                                        + (awpd2_new-0.5)*logdetpart2_U0_plus + logdetU0_plus/2.0;
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
    double val = scoreconst_2 - (awpNd2_new-0.5)*logdetpart2_UN - std::log(DUN)/2.0
                                  + (awpd2_new-0.5)*logdetpart2_U0 + std::log(DU0)/2.0;
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
      // R's "0" branch: base score uses logedgepvec[lp+1] (matching every
      // other lp value's own convention), and the +1-extension ALSO uses
      // logedgepvec[lp+1] specifically at lp=0 (not [lp+2], since there's
      // no "smaller" state to shift from) -- both restored/fixed to match
      // dagwishart_score_node's corrected lp=0 branch.
      double pen_base_lp0 = has_penalty ? logedgepvec(lp) : 0.0;
      double pen_plus_lp0 = has_penalty ? logedgepvec(lp+1) : 0.0;
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

// ============================================================================
// Incremental birth extension for BGe. See conversation notes: when node j
// gains a single new candidate parent e, the new (K+1)-candidate table's
// old-half rows and new-half base column are already present in the OLD
// table (e was a "+1" candidate before the birth) -- only the new-half
// rows' "+1" columns (against remaining candidates) are genuinely new,
// computed here via a DFS rooted at e rather than the empty set, visiting
// exactly 2^K leaves (half of a full from-scratch rebuild). Nothing
// persists across MCMC iterations; self-contained within one call.
// ============================================================================

static void extend_recurse(int depth, int K, std::vector<int>& included_local,
                           const arma::mat& choltemp, const arma::vec& c_noplus, double logdetD,
                           int j0, const arma::uvec& old_cand0, const arma::uvec& remaining_plus0,
                           const arma::mat& TN, double A, int e0, double awpN, int N,
                           const arma::vec& scoreconstvec, const IntegerVector& new_maps_backwards,
                           const std::vector<int>& new_pos_of_old, int e_new_pos0,
                           arma::mat& score_matr, const IntegerVector& old_maps_backwards,
                           const NumericMatrix& old_score_matr, int e_old_col){
  if(depth > K){
    long long new_bitmask = (1LL << e_new_pos0);
    long long old_bitmask = 0;
    for(int idx : included_local){ new_bitmask += (1LL << new_pos_of_old[idx-1]); old_bitmask += (1LL << (idx-1)); }
    int row_idx = new_maps_backwards[new_bitmask] - 1;
    int old_row_idx = old_maps_backwards[old_bitmask] - 1;
    double base_score_copy = old_score_matr(old_row_idx, e_old_col+1);

    int lpp = remaining_plus0.n_elem;
    score_matr(row_idx, 0) = base_score_copy;
    if(lpp == 0) return;

    int lp = 1 + (int)included_local.size();  // |S U {e}|, varies per leaf
    double awpNd2 = (awpN - N + lp + 1) / 2.0;
    double scoreconst_plus = scoreconstvec(lp+1); // R's scoreconstvec[lp+2] (1-idx) -> (lp+1) 0-idx
    arma::uvec parentnodes0(lp);
    parentnodes0(0) = e0;
    for(size_t a=0; a<included_local.size(); a++) parentnodes0(a+1) = old_cand0(included_local[a]-1);

    arma::vec TN_j_plus(lpp), TN_plus_diag(lpp);
    arma::mat TN_plus_off(lp, lpp);
    for(int b=0; b<lpp; b++){
      TN_j_plus(b) = TN(j0, remaining_plus0(b));
      TN_plus_diag(b) = TN(remaining_plus0(b), remaining_plus0(b));
      for(int a=0; a<lp; a++) TN_plus_off(a,b) = TN(parentnodes0(a), remaining_plus0(b));
    }

    double c_sq = arma::as_scalar(arma::sum(arma::square(c_noplus)));
    arma::mat choltemp_12 = arma::solve(arma::trimatl(choltemp.t()), TN_plus_off);
    for(int b=0; b<lpp; b++){
      double s = arma::as_scalar(arma::sum(arma::square(choltemp_12.col(b))));
      double u = std::sqrt(TN_plus_diag(b) - s);
      double cplus = (TN_j_plus(b) - arma::as_scalar(arma::dot(choltemp_12.col(b), c_noplus))) / u;
      double logdetD_plus = logdetD + 2*std::log(u);
      double logdetpart2_plus = std::log(A - c_sq - cplus*cplus);
      score_matr(row_idx, b+1) = scoreconst_plus - (awpNd2+0.5)*logdetpart2_plus - logdetD_plus/2.0;
    }
    return;
  }

  extend_recurse(depth+1, K, included_local, choltemp, c_noplus, logdetD,
                j0, old_cand0, remaining_plus0, TN, A, e0, awpN, N, scoreconstvec,
                new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                old_maps_backwards, old_score_matr, e_old_col);

  int new_node0 = old_cand0(depth-1);
  int lp_old = 1 + (int)included_local.size();  // includes e already
  double d = TN(new_node0, new_node0);

  arma::vec v(lp_old);
  v(0) = TN(e0, new_node0);
  for(size_t a=0; a<included_local.size(); a++) v(a+1) = TN(old_cand0(included_local[a]-1), new_node0);
  arma::vec w = arma::solve(arma::trimatl(choltemp.t()), v);
  double u = std::sqrt(d - arma::as_scalar(arma::sum(arma::square(w))));

  arma::mat new_choltemp(lp_old+1, lp_old+1, arma::fill::zeros);
  new_choltemp.submat(0,0,lp_old-1,lp_old-1) = choltemp;
  new_choltemp.submat(0,lp_old,lp_old-1,lp_old) = w;
  new_choltemp(lp_old,lp_old) = u;

  double b_new = TN(j0, new_node0);
  double new_c_last = (b_new - arma::as_scalar(arma::dot(w, c_noplus))) / u;
  arma::vec new_c(lp_old+1);
  new_c.subvec(0,lp_old-1) = c_noplus;
  new_c(lp_old) = new_c_last;
  double new_logdetD = logdetD + 2*std::log(u);

  included_local.push_back(depth);
  extend_recurse(depth+1, K, included_local, new_choltemp, new_c, new_logdetD,
                j0, old_cand0, remaining_plus0, TN, A, e0, awpN, N, scoreconstvec,
                new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                old_maps_backwards, old_score_matr, e_old_col);
  included_local.pop_back();
}

//' Incrementally extends an existing BGe "+1" score table when node j gains
//' a single new candidate parent e, reusing the old table for everything
//' that doesn't require new computation (see extend_recurse above).
//' @param old_cand_nodes K old candidates (1-indexed), OLD local order
//' @param old_plus_nodes old "+1" candidates (1-indexed), OLD column order
//' @param new_candidate   e (1-indexed) -- the newly birthed parent
//' @param new_cand_nodes  K+1 new candidates (1-indexed), SORTED
// [[Rcpp::export]]
NumericMatrix extend_plus_score_table_cpp(int j, IntegerVector old_cand_nodes, IntegerVector old_plus_nodes,
                                          int new_candidate, NumericMatrix old_score_matr_r,
                                          IntegerVector old_maps_backwards, IntegerVector new_cand_nodes,
                                          IntegerVector new_maps_backwards,
                                          NumericMatrix TN_r, double awpN, NumericVector scoreconstvec_r){
  int K = old_cand_nodes.size();
  int j0 = j-1;
  int e0 = new_candidate - 1;

  arma::mat TN = as<arma::mat>(TN_r);
  arma::vec scoreconstvec = as<arma::vec>(scoreconstvec_r);
  double A = TN(j0,j0);
  int N = TN_r.nrow();

  arma::uvec old_cand0(K);
  for(int k=0; k<K; k++) old_cand0(k) = old_cand_nodes[k]-1;

  int old_lpp = old_plus_nodes.size();
  int e_old_col = -1;
  for(int b=0; b<old_lpp; b++) if(old_plus_nodes[b]-1 == e0) e_old_col = b;
  // e_old_col is 0-indexed among "+1" candidates; actual matrix column is e_old_col+1

  std::vector<int> remaining_idx; remaining_idx.reserve(old_lpp-1);
  for(int b=0; b<old_lpp; b++) if(b != e_old_col) remaining_idx.push_back(b);
  int lpp = remaining_idx.size();
  arma::uvec remaining_plus0(lpp);
  std::vector<int> remaining_old_cols(lpp);
  for(int b=0; b<lpp; b++){
    remaining_plus0(b) = old_plus_nodes[remaining_idx[b]]-1;
    remaining_old_cols[b] = remaining_idx[b];
  }

  int Knew = new_cand_nodes.size(); // K+1
  int e_new_pos0 = -1;
  for(int k=0; k<Knew; k++) if(new_cand_nodes[k]-1 == e0) e_new_pos0 = k;
  std::vector<int> new_pos_of_old(K);
  for(int k=0; k<K; k++){
    for(int m=0; m<Knew; m++) if(new_cand_nodes[m]-1 == old_cand0(k)){ new_pos_of_old[k] = m; break; }
  }

  long long n_new_rows = 1LL << Knew;
  arma::mat score_matr(n_new_rows, lpp+1);

  // ---- fill OLD-half rows directly (no computation) ----
  // Build source/target row-index vectors once (index arithmetic only, no
  // data movement), then do the actual copy as two vectorized Armadillo
  // calls: gather the relevant old rows/columns, then scatter into place.
  long long n_old_rows = 1LL << K;
  arma::uvec source_rows(n_old_rows), target_rows(n_old_rows);
  for(long long old_bm=0; old_bm<n_old_rows; old_bm++){
    source_rows(old_bm) = old_maps_backwards[old_bm]-1;
    long long new_bm = 0;
    for(int k=0; k<K; k++) if(old_bm & (1LL << k)) new_bm += (1LL << new_pos_of_old[k]);
    target_rows(old_bm) = new_maps_backwards[new_bm]-1;
  }
  arma::uvec cols_to_copy(lpp+1);
  cols_to_copy(0) = 0;
  for(int b=0; b<lpp; b++) cols_to_copy(b+1) = remaining_old_cols[b]+1;

  arma::mat old_score_matr_arma = as<arma::mat>(old_score_matr_r);
  arma::mat gathered = old_score_matr_arma.submat(source_rows, cols_to_copy);
  score_matr.rows(target_rows) = gathered;

  // ---- fill NEW-half rows via DFS rooted at e ----
  double u_e = std::sqrt(TN(e0,e0));
  arma::mat choltemp_e(1,1); choltemp_e(0,0) = u_e;
  arma::vec c_e(1); c_e(0) = TN(j0,e0)/u_e;
  double logdetD_e = 2*std::log(u_e);

  // awpNd2 and scoreconst for the (K+1)-length ("+1" over K+1 base) case --
  // matches bge_score_plus_parent's convention: awpNd2 uses lp = K+1 (the
  // NEW base size), scoreconstvec index is lp+2 = K+3 (1-indexed) -> K+2 (0-indexed)
  std::vector<int> included_local; included_local.reserve(K);
  extend_recurse(1, K, included_local, choltemp_e, c_e, logdetD_e,
                j0, old_cand0, remaining_plus0, TN, A, e0, awpN, N, scoreconstvec,
                new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                old_maps_backwards, old_score_matr_r, e_old_col);

  return wrap(score_matr);
}

// ============================================================================
// DAG-Wishart analog of extend_plus_score_table_cpp -- same free-reuse
// insight (old-half rows and new-half base column copied straight from the
// old table), extended to track both UN and U0 factors in parallel, and
// reusing dagwishart_score_from_factors for leaf scoring (lp is never 0
// here, since e is always included, so the lp=0 penalty asymmetry never
// applies in this function).
// ============================================================================

static void dagwishart_extend_recurse(int depth, int K, std::vector<int>& included_local,
                                      const arma::mat& cholUN, const arma::vec& cUN, double logdetUN,
                                      const arma::mat& cholU0, const arma::vec& cU0, double logdetU0,
                                      int j0, const arma::uvec& old_cand0, const arma::uvec& remaining_plus0,
                                      const arma::mat& UN, const arma::mat& U0, double A, double A0,
                                      int e0, double awpN_new, int N, const arma::vec& scoreconstlist_lp,
                                      const IntegerVector& new_maps_backwards,
                                      const std::vector<int>& new_pos_of_old, int e_new_pos0,
                                      arma::mat& score_matr, const IntegerVector& old_maps_backwards,
                                      const NumericMatrix& old_score_matr, int e_old_col,
                                      bool has_penalty, const arma::vec& logedgepvec){
  if(depth > K){
    long long new_bitmask = (1LL << e_new_pos0);
    long long old_bitmask = 0;
    for(int idx : included_local){ new_bitmask += (1LL << new_pos_of_old[idx-1]); old_bitmask += (1LL << (idx-1)); }
    int row_idx = new_maps_backwards[new_bitmask] - 1;
    (void)old_maps_backwards; (void)old_score_matr; (void)e_old_col; (void)old_bitmask;
    // NOTE: unlike BGe, DAG-Wishart's existing "+1"-from-empty-base formula
    // does not equal the direct base-score formula for the same expanded
    // parent set (verified against dagwishart_score_node -- a pre-existing
    // property of dagwishart_score_plus_parent, not something introduced
    // here). So the base score for S U {e} is computed fresh via
    // dagwishart_score_from_factors below, matching what
    // build_plus_score_table_dagwishart_cpp would produce, rather than
    // copied from the old table's e-column.

    int lpp = remaining_plus0.n_elem;
    int lp = 1 + (int)included_local.size();  // |S U {e}|, varies per leaf
    double pen_base = has_penalty ? logedgepvec(lp) : 0.0;
    double pen_plus = has_penalty ? (lpp>0 ? logedgepvec(lp+1) : 0.0) : 0.0;

    arma::uvec parentnodes0(lp);
    parentnodes0(0) = e0;
    for(size_t a=0; a<included_local.size(); a++) parentnodes0(a+1) = old_cand0(included_local[a]-1);

    double scoreconst_lp = scoreconstlist_lp(lp);
    double scoreconst_lpplus1 = lpp>0 ? scoreconstlist_lp(lp+1) : 0.0;

    if(lpp == 0){
      arma::uvec empty_plus0;
      arma::rowvec row_scores = dagwishart_score_from_factors(j0, parentnodes0, empty_plus0, UN, U0,
                                                               cholUN, cUN, logdetUN, cholU0, cU0, logdetU0,
                                                               A, A0, lp, awpN_new, N, scoreconst_lp,
                                                               scoreconst_lpplus1, has_penalty, pen_base, pen_plus);
      score_matr(row_idx, 0) = row_scores(0);
      return;
    }

    arma::rowvec row_scores = dagwishart_score_from_factors(j0, parentnodes0, remaining_plus0, UN, U0,
                                                             cholUN, cUN, logdetUN, cholU0, cU0, logdetU0,
                                                             A, A0, lp, awpN_new, N, scoreconst_lp,
                                                             scoreconst_lpplus1, has_penalty, pen_base, pen_plus);
    for(int b=0; b<=lpp; b++) score_matr(row_idx, b) = row_scores(b);
    return;
  }

  dagwishart_extend_recurse(depth+1, K, included_local, cholUN, cUN, logdetUN, cholU0, cU0, logdetU0,
                            j0, old_cand0, remaining_plus0, UN, U0, A, A0, e0, awpN_new, N, scoreconstlist_lp,
                            new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                            old_maps_backwards, old_score_matr, e_old_col, has_penalty, logedgepvec);

  int new_node0 = old_cand0(depth-1);
  int lp_old = 1 + (int)included_local.size();
  double dUN = UN(new_node0, new_node0);
  double dU0 = U0(new_node0, new_node0);

  arma::vec vUN(lp_old), vU0(lp_old);
  vUN(0) = UN(e0, new_node0); vU0(0) = U0(e0, new_node0);
  for(size_t a=0; a<included_local.size(); a++){
    int node0 = old_cand0(included_local[a]-1);
    vUN(a+1) = UN(node0, new_node0);
    vU0(a+1) = U0(node0, new_node0);
  }
  arma::vec wUN = arma::solve(arma::trimatl(cholUN.t()), vUN);
  arma::vec wU0 = arma::solve(arma::trimatl(cholU0.t()), vU0);
  double uUN = std::sqrt(dUN - arma::as_scalar(arma::sum(arma::square(wUN))));
  double uU0 = std::sqrt(dU0 - arma::as_scalar(arma::sum(arma::square(wU0))));

  arma::mat new_cholUN(lp_old+1, lp_old+1, arma::fill::zeros);
  new_cholUN.submat(0,0,lp_old-1,lp_old-1) = cholUN;
  new_cholUN.submat(0,lp_old,lp_old-1,lp_old) = wUN;
  new_cholUN(lp_old,lp_old) = uUN;

  arma::mat new_cholU0(lp_old+1, lp_old+1, arma::fill::zeros);
  new_cholU0.submat(0,0,lp_old-1,lp_old-1) = cholU0;
  new_cholU0.submat(0,lp_old,lp_old-1,lp_old) = wU0;
  new_cholU0(lp_old,lp_old) = uU0;

  double bUN_new = UN(j0, new_node0);
  double bU0_new = U0(j0, new_node0);
  double new_cUN_last = (bUN_new - arma::as_scalar(arma::dot(wUN, cUN))) / uUN;
  double new_cU0_last = (bU0_new - arma::as_scalar(arma::dot(wU0, cU0))) / uU0;
  arma::vec new_cUN(lp_old+1); new_cUN.subvec(0,lp_old-1) = cUN; new_cUN(lp_old) = new_cUN_last;
  arma::vec new_cU0(lp_old+1); new_cU0.subvec(0,lp_old-1) = cU0; new_cU0(lp_old) = new_cU0_last;

  double new_logdetUN = logdetUN + 2*std::log(uUN);
  double new_logdetU0 = logdetU0 + 2*std::log(uU0);

  included_local.push_back(depth);
  dagwishart_extend_recurse(depth+1, K, included_local, new_cholUN, new_cUN, new_logdetUN,
                            new_cholU0, new_cU0, new_logdetU0,
                            j0, old_cand0, remaining_plus0, UN, U0, A, A0, e0, awpN_new, N, scoreconstlist_lp,
                            new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                            old_maps_backwards, old_score_matr, e_old_col, has_penalty, logedgepvec);
  included_local.pop_back();
}

//' DAG-Wishart analog of extend_plus_score_table_cpp -- see that function's
//' documentation for the free-reuse design; here two factors (UN, U0) are
//' tracked in parallel along the DFS rooted at e, mirroring
//' build_plus_score_table_dagwishart_cpp's own dual-factor structure.
// [[Rcpp::export]]
NumericMatrix extend_plus_score_table_dagwishart_cpp(int j, IntegerVector old_cand_nodes, IntegerVector old_plus_nodes,
                                                     int new_candidate, NumericMatrix old_score_matr_r,
                                                     IntegerVector old_maps_backwards, IntegerVector new_cand_nodes,
                                                     IntegerVector new_maps_backwards,
                                                     NumericMatrix UN_r, NumericMatrix U0_r, double awpN_new,
                                                     List scoreconstlist, NumericVector logedgepvec_r){
  int K = old_cand_nodes.size();
  int j0 = j-1;
  int e0 = new_candidate - 1;
  int N = UN_r.nrow();

  arma::mat UN = as<arma::mat>(UN_r);
  arma::mat U0 = as<arma::mat>(U0_r);
  double A = UN(j0,j0);
  double A0 = U0(j0,j0);

  arma::uvec old_cand0(K);
  for(int k=0; k<K; k++) old_cand0(k) = old_cand_nodes[k]-1;

  int old_lpp = old_plus_nodes.size();
  int e_old_col = -1;
  for(int b=0; b<old_lpp; b++) if(old_plus_nodes[b]-1 == e0) e_old_col = b;

  std::vector<int> remaining_idx; remaining_idx.reserve(old_lpp-1);
  for(int b=0; b<old_lpp; b++) if(b != e_old_col) remaining_idx.push_back(b);
  int lpp = remaining_idx.size();
  arma::uvec remaining_plus0(lpp);
  std::vector<int> remaining_old_cols(lpp);
  for(int b=0; b<lpp; b++){
    remaining_plus0(b) = old_plus_nodes[remaining_idx[b]]-1;
    remaining_old_cols[b] = remaining_idx[b];
  }

  int Knew = new_cand_nodes.size();
  int e_new_pos0 = -1;
  for(int k=0; k<Knew; k++) if(new_cand_nodes[k]-1 == e0) e_new_pos0 = k;
  std::vector<int> new_pos_of_old(K);
  for(int k=0; k<K; k++){
    for(int m=0; m<Knew; m++) if(new_cand_nodes[m]-1 == old_cand0(k)){ new_pos_of_old[k] = m; break; }
  }

  long long n_new_rows = 1LL << Knew;
  arma::mat score_matr(n_new_rows, lpp+1);

  // ---- fill OLD-half rows directly (vectorized gather-then-scatter) ----
  long long n_old_rows = 1LL << K;
  arma::uvec source_rows(n_old_rows), target_rows(n_old_rows);
  for(long long old_bm=0; old_bm<n_old_rows; old_bm++){
    source_rows(old_bm) = old_maps_backwards[old_bm]-1;
    long long new_bm = 0;
    for(int k=0; k<K; k++) if(old_bm & (1LL << k)) new_bm += (1LL << new_pos_of_old[k]);
    target_rows(old_bm) = new_maps_backwards[new_bm]-1;
  }
  arma::uvec cols_to_copy(lpp+1);
  cols_to_copy(0) = 0;
  for(int b=0; b<lpp; b++) cols_to_copy(b+1) = remaining_old_cols[b]+1;

  arma::mat old_score_matr_arma = as<arma::mat>(old_score_matr_r);
  arma::mat gathered = old_score_matr_arma.submat(source_rows, cols_to_copy);
  score_matr.rows(target_rows) = gathered;

  // ---- fill NEW-half rows via DFS rooted at e (both UN and U0 factors) ----
  double uUN_e = std::sqrt(UN(e0,e0));
  double uU0_e = std::sqrt(U0(e0,e0));
  arma::mat cholUN_e(1,1); cholUN_e(0,0) = uUN_e;
  arma::mat cholU0_e(1,1); cholU0_e(0,0) = uU0_e;
  arma::vec cUN_e(1); cUN_e(0) = UN(j0,e0)/uUN_e;
  arma::vec cU0_e(1); cU0_e(0) = U0(j0,e0)/uU0_e;
  double logdetUN_e = 2*std::log(uUN_e);
  double logdetU0_e = 2*std::log(uU0_e);

  bool has_penalty = logedgepvec_r.size() > 0;
  arma::vec logedgepvec = has_penalty ? as<arma::vec>(logedgepvec_r) : arma::vec();
  arma::vec scoreconstlist_j(N+3); // generous upper bound on lp+1 indexing; only entries up to K+2 actually used
  {
    int max_idx = std::min((int)scoreconstlist.size(), N+3);
    for(int m=0; m<max_idx; m++) scoreconstlist_j(m) = as<NumericVector>(scoreconstlist[m])[j0];
  }

  std::vector<int> included_local; included_local.reserve(K);
  dagwishart_extend_recurse(1, K, included_local, cholUN_e, cUN_e, logdetUN_e, cholU0_e, cU0_e, logdetU0_e,
                            j0, old_cand0, remaining_plus0, UN, U0, A, A0, e0, awpN_new, N, scoreconstlist_j,
                            new_maps_backwards, new_pos_of_old, e_new_pos0, score_matr,
                            old_maps_backwards, old_score_matr_r, e_old_col, has_penalty, logedgepvec);

  return wrap(score_matr);
}
