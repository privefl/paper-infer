/******************************************************************************/

#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/******************************************************************************/

double compute_add(const std::vector<size_t>& p,
                   const IntegerVector& i,
                   const NumericVector& x,
                   const NumericVector& tagged,
                   int j) {

  double sum_add = 0;

  for (size_t k = p[j]; k < p[j + 1]; k++) {

    int    ind = i[k];
    double val = x[k];

    double add = val - tagged[ind];
    if (add > 0) sum_add += add;
  }

  return sum_add;
}

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericVector> greedy_maxtag(const std::vector<size_t>& p,
                                    const IntegerVector& i,
                                    const NumericVector& x0,
                                    const LogicalVector& select,
                                    const LogicalVector& exclude,
                                    const NumericVector& sqrt_info,
                                    double min_add,
                                    bool remove_diag) {

  NumericVector x = x0 * x0;
  LogicalVector removed = Rcpp::clone(exclude);

  int m = p.size() - 1;
  int m_remain = Rcpp::sum(!removed);

  NumericVector can_add(m, 0.0);
  for (int j = 0; j < m; j++) {
    size_t lo = p[j];
    size_t up = p[j + 1];
    for (size_t k = lo; k != up; k++) {
      if (remove_diag && i[k] == j) {
        x[k] = 0;
      } else {
        x[k] *= sqrt_info[j];
        if (!removed[j]) can_add[j] += x[k];
      }
    }
  }
  can_add[select] = R_PosInf;

  NumericVector tagged(m, 0.0), added(m, NA_REAL);

  std::unordered_set<int> ind_has_changed, ind_could_change;

  do {

    int j_max = which_max(can_add);
    if (can_add[j_max] < min_add) break;

    // Rcout << j_max + 1 << " -> " << can_add[j_max] << std::endl;

    added[j_max] = can_add[j_max];
    can_add[j_max] = 0;
    removed[j_max] = true;
    m_remain--;

    // update tagged and can_add for the variants correlated with j_max

    // update tagged
    size_t lo = p[j_max];
    size_t up = p[j_max + 1];
    for (size_t k = lo; k != up; k++) {

      int    ind = i[k];
      double val = x[k];

      if (val > tagged[ind]) {
        tagged[ind] = val;
        ind_has_changed.insert(ind);
      }
    }

    for (auto& j : ind_has_changed) {
      for (size_t k = p[j]; k < p[j + 1]; k++) {
        int ind = i[k];
        if (!removed[ind] & !select[ind]) ind_could_change.insert(ind);
      }
    }


    // update can_add
    for (auto& ind : ind_could_change) {

      double sum_add = compute_add(p, i, x, tagged, ind);
      if (sum_add < min_add) {
        can_add[ind] = 0;
        removed[ind] = true;
        m_remain--;
      } else {
        can_add[ind] = sum_add;
      }
    }

    ind_has_changed.clear();
    ind_could_change.clear();

  }  while (m_remain > 0);

  return List::create(added, tagged);
}

/******************************************************************************/
