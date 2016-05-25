#pragma once

#include "basic_types.hpp"

namespace pano {
namespace core {
struct FactorCategory {
  using CostFunction = std::function<double(const std::vector<int> &var_labels,
                                            void *additional_data)>;
  CostFunction cost;
  double c_alpha;
};
struct VarCategory {
  size_t nlabels;
  double c_i;
};

// factor graph
class FactorGraph {
  std::vector<FactorCategory> factor_cats_;
  std::vector<VarCategory> var_cats_;

  std::vector<int> var2cat_;
  std::vector<std::set<int>> var2factors_;

  std::vector<int> factor2cat_;
  std::vector<std::vector<int>> factor2vars_;

public:
  void reserveFactorCategories(size_t cap);
  void reserveVarCategories(size_t cap);

  // returns var_cat
  int addVarCategory(size_t nlabels, double c_i);
  // returns factor_cat
  int addFactorCategory(FactorCategory::CostFunction cost, double c_alpha);

  void reserveFactors(size_t cap);
  void reserveVars(size_t cap);

  // returns var
  int addVar(int var_cat);
  // returns factor
  int addFactor(int factor_cat, const std::vector<int> &vars);
  int addFactor(int factor_cat, std::vector<int> &&vars);
  int addFactor(int factor_cat, std::initializer_list<int> vars);

  inline size_t nvars() const { return var2cat_.size(); }
  inline size_t nfactors() const { return factor2cat_.size(); }

  void clear();
  bool valid() const;

  double cost(const std::vector<int> &var_labels,
              void *additional_data = nullptr) const;

  // convex belief propagation
  std::vector<int>
  solve(int max_epoch, int inner_loop_num = 10,
        std::function<bool(int epoch, double energy, double denergy,
                           const std::vector<int> &cur_best_var_labels)>
            callback = nullptr,
        void *additional_data = nullptr) const;

  template <class CallbackFunT>
  auto solve(int max_epoch, int inner_loop_num, CallbackFunT callback,
             void *additional_data = nullptr) const
      -> decltype(callback(0, 0.0), std::vector<int>()) {
    return solve(
        max_epoch, inner_loop_num,
        [callback](int epoch, double energy, double denergy,
                   const std::vector<int> &cur_best_var_labels) -> bool {
          return callback(epoch, energy);
        },
        additional_data);
  }
};
}
}
