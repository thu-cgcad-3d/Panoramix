#include "factor_graph.hpp"
#include "../misc/eigen.hpp"

namespace pano {
namespace core {
void FactorGraph::reserveFactorCategories(size_t cap) {
  factor_cats_.reserve(cap);
}
void FactorGraph::reserveVarCategories(size_t cap) { var_cats_.reserve(cap); }

int FactorGraph::addVarCategory(size_t nlabels, double c_i) {
  assert(nlabels > 0);
  var_cats_.push_back(VarCategory{nlabels, c_i});
  return static_cast<int>(var_cats_.size()) - 1;
}
int FactorGraph::addFactorCategory(FactorCategory::CostFunction cost,
                                   double c_alpha) {
  assert(cost);
  factor_cats_.push_back(FactorCategory{cost, c_alpha});
  return static_cast<int>(factor_cats_.size()) - 1;
}

void FactorGraph::reserveFactors(size_t cap) {
  factor2cat_.reserve(cap);
  factor2vars_.reserve(cap);
}

void FactorGraph::reserveVars(size_t cap) {
  var2cat_.reserve(cap);
  var2factors_.reserve(cap);
}

int FactorGraph::addVar(int var_cat) {
  var2cat_.push_back(var_cat);
  var2factors_.push_back({});
  return var2cat_.size() - 1;
}
int FactorGraph::addFactor(int factor_cat, const std::vector<int> &vars) {
  int factor = factor2cat_.size();
  for (int var : vars) {
    var2factors_[var].insert(factor);
  }
  factor2cat_.push_back(factor_cat);
  factor2vars_.push_back(vars);
  return factor;
}
int FactorGraph::addFactor(int factor_cat, std::vector<int> &&vars) {
  int factor = factor2cat_.size();
  for (int var : vars) {
    var2factors_[var].insert(factor);
  }
  factor2cat_.push_back(factor_cat);
  factor2vars_.push_back(std::move(vars));
  return factor;
}
int FactorGraph::addFactor(int factor_cat, std::initializer_list<int> vars) {
  return addFactor(factor_cat, std::vector<int>(vars));
}

void FactorGraph::clear() {
  factor_cats_.clear();
  var_cats_.clear();
  var2cat_.clear();
  var2factors_.clear();
  factor2cat_.clear();
  factor2vars_.clear();
}

bool FactorGraph::valid() const {
  if (var2cat_.size() != var2factors_.size()) {
    return false;
  }
  if (factor2cat_.size() != factor2vars_.size()) {
    return false;
  }
  for (int var_cat : var2cat_) {
    if (var_cat < 0 || var_cat >= var_cats_.size()) {
      return false;
    }
  }
  for (auto &factors : var2factors_) {
    for (int factor : factors) {
      if (factor < 0 || factor >= factor2cat_.size()) {
        return false;
      }
    }
  }
  return true;
}

double FactorGraph::cost(const std::vector<int> &var_labels,
                         void *additional_data) const {
  assert(valid());
  double cost = 0.0;
  for (int factor = 0; factor < nfactors(); factor++) {
    const auto &vars = factor2vars_[factor];
    std::vector<int> local_var_labels;
    local_var_labels.reserve(vars.size());
    for (int var : vars) {
      local_var_labels.push_back(var_labels[var]);
    }
    double factor_cost = factor_cats_[factor2cat_[factor]].cost(
        local_var_labels, additional_data);
    cost += factor_cost;
  }
  return cost;
}

std::vector<int> FactorGraph::solve(
    int max_epoch, int inner_loop_num,
    std::function<bool(int epoch, double energy, double denergy,
                       const std::vector<int> &cur_best_var_labels)>
        callback,
    void *additional_data) const {

  assert(valid());
  using VecX = Eigen::VectorXd;

  size_t nconnections = 0;
  for (int factor = 0; factor < nfactors(); factor++) {
    nconnections += factor2vars_[factor].size();
  }

  std::vector<int> con2factor(nconnections, -1);
  std::vector<int> con2var(nconnections, -1);
  std::vector<VecX> fv_messages(nconnections); // messages from factor to var
  std::vector<VecX> vf_messages(nconnections); // messages from var to factor

  std::vector<std::vector<int>> var2cons(nvars());
  std::vector<std::vector<int>> factor2cons(nfactors());

  for (int factor = 0, con = 0; factor < nfactors(); factor++) {
    for (int var : factor2vars_[factor]) {
      con2factor[con] = factor;
      con2var[con] = var;
      size_t var_nlabels = var_cats_[var2cat_[var]].nlabels;
      fv_messages[con] = vf_messages[con] = VecX::Zero(var_nlabels);
      var2cons[var].push_back(con);
      factor2cons[factor].push_back(con);
      con++;
    }
  }

  // precompute c_i_hats for each var
  std::vector<double> c_i_hats(nvars());
  for (int var = 0; var < nvars(); var++) {
    double &c_i_hat = c_i_hats[var];
    c_i_hat = var_cats_[var2cat_[var]].c_i;
    for (auto factor : var2factors_[var]) {
      c_i_hat += factor_cats_[factor2cat_[factor]].c_alpha;
    }
  }

  // precompute var dimensions for each factor
  std::vector<std::vector<size_t>> factor2var_dims(nfactors());
  for (int factor = 0; factor < nfactors(); factor++) {
    auto &var_dims = factor2var_dims[factor];
    var_dims.reserve(factor2vars_[factor].size());
    for (int var : factor2vars_[factor]) {
      var_dims.push_back(var_cats_[var2cat_[var]].nlabels);
    }
  }

  // var marginals
  std::vector<VecX> var2marginals(nvars());
  for (int var = 0; var < nvars(); var++) {
    var2marginals[var] = VecX::Zero(var_cats_[var2cat_[var]].nlabels);
  }

  double last_cost = std::numeric_limits<double>::infinity();
  std::vector<int> results(nvars(), 0);

  for (int epoch = 0; epoch < max_epoch; epoch++) {
    for (int l = 0; l < inner_loop_num; l++) {
      // var -> factor
      for (int vf_con = 0; vf_con < nconnections; vf_con++) {
        int var = con2var[vf_con];
        int factor = con2factor[vf_con];

        vf_messages[vf_con].setZero();
        int oppose_fv_con = -1;
        for (int fv_con : var2cons[var]) {
          vf_messages[vf_con] += fv_messages[fv_con];
          assert(!vf_messages[vf_con].hasNaN());
          if (con2factor[fv_con] == factor) {
            oppose_fv_con = fv_con;
          }
        }

        assert(oppose_fv_con != -1);

        auto &factor_cat_data = factor_cats_[factor2cat_[factor]];
        vf_messages[vf_con] =
            vf_messages[vf_con] * factor_cat_data.c_alpha / c_i_hats[var] -
            fv_messages[oppose_fv_con];
        assert(!vf_messages[vf_con].hasNaN());
      }

      // factor -> var
      for (int factor = 0; factor < nfactors(); factor++) {
        auto &vars = factor2vars_[factor];
        auto &cons_here = factor2cons[factor];
        assert(cons_here.size() == vars.size());

        // reset fv messages
        for (int con : cons_here) {
          fv_messages[con].setConstant(std::numeric_limits<double>::infinity());
        }

        // dispatch messages from this factor
        auto &factor_cat_data = factor_cats_[factor2cat_[factor]];
        auto &cost_fun = factor_cat_data.cost;
        auto &dims = factor2var_dims[factor];
        assert(dims.size() > 0);
        std::vector<int> idx(dims.size(), 0);

        while (true) {
          for (int i = 0; i < dims.size(); i++) {
            // i: output var
            // others: input _varCategories
            double &out_value = fv_messages[cons_here[i]][idx[i]];
            // theta
            double theta = cost_fun(idx, additional_data);
            assert(!std::isnan(theta));
            assert(theta >= 0);

            // sum of other input var_cats
            double sum_of_other_vars = 0.0;
            for (int j = 0; j < dims.size(); j++) {
              if (j == i)
                continue;
              sum_of_other_vars += vf_messages[cons_here[j]][idx[j]];
            }

            assert(!std::isnan(sum_of_other_vars));
            double score = theta + sum_of_other_vars;
            assert(!std::isnan(score));
            assert(score <= std::numeric_limits<double>::infinity());
            if (score < out_value) {
              out_value = score;
            }
          } // for (int i = 0; i < dims.size(); i++)

          // add idx
          idx[0]++;
          int k = 0;
          while (k < idx.size() && idx[k] >= dims[k]) {
            idx[k] = 0;
            if (k + 1 < idx.size())
              idx[k + 1]++;
            k++;
          }
          if (k == idx.size()) {
            break;
          }
        } // while (true)

        for (auto &message : fv_messages) {
          assert(!message.hasNaN());
        }
        for (auto &message : vf_messages) {
          assert(!message.hasNaN());
        }
      }
    } // for (int l = 0; l < inner_loop_num; l++)

    // marginalize on variables
    // reset marginals
    for (auto &marginal : var2marginals) {
      marginal.setZero();
    }
    // update marginals
    for (int con = 0; con < nconnections; con++) {
      int var = con2var[con];
      var2marginals[var] += fv_messages[con];
      assert(var2marginals[var].maxCoeff() <
             std::numeric_limits<double>::infinity());
    }

    // get resulted labels
    for (int var = 0; var < nvars(); var++) {
      int label = -1;
      double cur_cost = std::numeric_limits<double>::infinity();
      const VecX &marginal = var2marginals[var];
      for (int i = 0; i < marginal.size(); i++) {
        if (marginal[i] < cur_cost) {
          label = i;
          cur_cost = marginal[i];
        }
      }
      assert(label != -1);
      results[var] = label;
    }

    // compute current holistic cost and invoke callback
    double c = cost(results, additional_data);
    if (callback && !callback(epoch, c, c - last_cost, results)) {
      break;
    }
    last_cost = c;
  }

  return results;
}
}
}
