#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "../core/cons_graph.hpp"
#include "../core/utility.hpp"
#include "factor_graph.hpp"

namespace pano {
namespace ml {

using namespace core;

bool FactorGraph::valid() const {
  for (auto &v : _graph.elements<0>()) {
    auto &vid = v.data;
    if (vid < 0 || vid >= _varCategories.size())
      return false;
  }
  for (auto &f : _graph.elements<1>()) {
    auto &fid = f.data;
    auto &nhs = f.topo.lowers;
    if (fid < 0 || fid >= _factorCategories.size())
      return false;
    auto &costFun = _factorCategories[fid].costs;
    if (!costFun)
      return false;
  }
  return true;
}

double FactorGraph::energy(const ResultTable &labels, void *givenData) const {
  assert(valid());
  double e = 0.0;
  for (auto &f : _graph.elements<1>()) {
    auto &vhs = f.topo.lowers;
    std::vector<int> idx(vhs.size());
    for (int i = 0; i < idx.size(); i++) {
      idx[i] = labels[vhs[i]];
    }
    double factorCost = _factorCategories[f.data].costs(idx.data(), idx.size(),
                                                        f.data, givenData);
    e += factorCost;
  }
  return e;
}

FactorGraph::ResultTable FactorGraph::solve(int maxEpoch, int innerLoopNum,
                                            const CallbackFunction &callback,
                                            void *givenData) const {

  // convex belief propagation

  assert(valid());
  using namespace Eigen;

  struct VHandleWrapper {
    VarHandle h;
  };
  struct FHandleWrapper {
    FactorHandle h;
  };
  struct F2VMessage {
    VectorXd values;
  };
  struct V2FMessage {
    VectorXd values;
  };

  ConstraintGraph<
      std::tuple<VHandleWrapper, FHandleWrapper>,
      std::tuple<ConstraintConfig<V2FMessage, VHandleWrapper, FHandleWrapper>,
                 ConstraintConfig<F2VMessage, FHandleWrapper, VHandleWrapper>>>
      messages;
  {
    messages.internalComponents<VHandleWrapper>().reserve(
        _graph.internalElements<0>().size());
    messages.internalComponents<FHandleWrapper>().reserve(
        _graph.internalElements<1>().size());
    size_t vfConnectionNum = 0;
    for (auto &f : _graph.elements<1>()) {
      vfConnectionNum += f.topo.lowers.size();
    }
    messages.internalConstraints<V2FMessage>().reserve(vfConnectionNum);
    messages.internalConstraints<F2VMessage>().reserve(vfConnectionNum);
  }

  using MGVHandle = ComponentHandle<VHandleWrapper>;
  using MGFHandle = ComponentHandle<FHandleWrapper>;
  using MGV2FHandle = ConstraintHandle<V2FMessage>;
  using MGF2VHandle = ConstraintHandle<F2VMessage>;

  // install _graph to message _graph and initialize messages
  HandledTable<VarHandle, MGVHandle> vhToMGVH(
      _graph.internalElements<0>().size());
  for (auto &v : _graph.elements<0>()) {
    MGVHandle vh = messages.addComponent(VHandleWrapper{v.topo.hd});
    vhToMGVH[v.topo.hd] = vh;
  }
  for (auto &f : _graph.elements<1>()) {
    MGFHandle fh = messages.addComponent(FHandleWrapper{f.topo.hd});
    for (auto &lh : f.topo.lowers) {
      MGVHandle vh = vhToMGVH[lh];
      messages.addConstraint(
          V2FMessage{VectorXd::Zero(_varCategories[_graph.data(lh)].nlabels)},
          vh, fh);
      messages.addConstraint(
          F2VMessage{VectorXd::Zero(_varCategories[_graph.data(lh)].nlabels)},
          fh, vh);
    }
  }

  // precompute c_i_hats for each var
  HandledTable<VarHandle, double> c_i_hats(_graph.internalElements<0>().size(),
                                           0.0);
  for (auto &v : _graph.elements<0>()) {
    double &c_i_hat = c_i_hats[v.topo.hd];
    c_i_hat = _varCategories[v.data].c_i;
    for (auto fh : v.topo.uppers) {
      c_i_hat += _factorCategories[_graph.data(fh)].c_alpha;
    }
  }

  // store ordered f2v and v2f handles for each factor
  core::HandledTable<MGFHandle, std::vector<MGF2VHandle>> orderedF2Vmsghs(
      _graph.internalElements<1>().size());
  core::HandledTable<MGFHandle, std::vector<MGV2FHandle>> orderedV2Fmsghs(
      _graph.internalElements<1>().size());
  // store var dimensions of each factor
  core::HandledTable<MGFHandle, std::vector<size_t>> orderedVDims(
      _graph.internalElements<1>().size());
  for (auto &f : messages.components<FHandleWrapper>()) {
    MGFHandle mgfh = f.topo.hd;
    FactorHandle fh = f.data.h;
    auto &vhs = _graph.topo(fh).lowers;
    assert(!vhs.empty());
    auto &v2fmsghsSet = messages.topo(mgfh).constraints<V2FMessage>();
    auto &v2fmsghs = orderedV2Fmsghs[mgfh];
    v2fmsghs.resize(vhs.size());
    auto &f2vmsghsSet = messages.topo(mgfh).constraints<F2VMessage>();
    auto &f2vmsghs = orderedF2Vmsghs[mgfh];
    f2vmsghs.resize(vhs.size());
    for (int i = 0; i < vhs.size(); i++) {
      for (MGV2FHandle v2fmsgh : v2fmsghsSet) {
        if (messages.topo(v2fmsgh).component<0>() == vhToMGVH[vhs[i]]) {
          v2fmsghs[i] = v2fmsgh;
          break;
        }
      }
      for (MGF2VHandle f2vmsgh : f2vmsghsSet) {
        if (messages.topo(f2vmsgh).component<1>() == vhToMGVH[vhs[i]]) {
          f2vmsghs[i] = f2vmsgh;
          break;
        }
      }
    }
    orderedVDims[mgfh].resize(vhs.size());
    for (int i = 0; i < vhs.size(); i++) {
      orderedVDims[mgfh][i] = _varCategories[_graph.data(vhs[i])].nlabels;
    }
  }

  // initialize marginals
  core::HandledTable<VarHandle, VectorXd> varMarginals(
      _graph.internalElements<0>().size());
  for (auto &v : _graph.elements<0>()) {
    varMarginals[v.topo.hd] = VectorXd::Zero(_varCategories[v.data].nlabels);
  }

  double lastE = std::numeric_limits<double>::infinity();
  ResultTable results(_graph.internalElements<0>().size());

  for (int epoch = 0; epoch < maxEpoch; epoch++) {
    for (int l = 0; l < innerLoopNum; l++) {
      // var -> factor
      for (auto &v2f : messages.constraints<V2FMessage>()) {
        MGVHandle vh = v2f.topo.component<0>();
        MGFHandle fh = v2f.topo.component<1>();

        auto &fdata = _factorCategories[_graph.data(messages.data(fh).h)];

        v2f.data.values.setZero();
        MGF2VHandle oppose;
        for (MGF2VHandle f2v : messages.topo(vh).constraints<F2VMessage>()) {
          v2f.data.values += messages.data(f2v).values;
          assert(!v2f.data.values.hasNaN());
          if (messages.topo(f2v).component<0>() == fh) {
            oppose = f2v;
          }
        }
        assert(oppose.valid());
        v2f.data.values =
            v2f.data.values * fdata.c_alpha / c_i_hats[messages.data(vh).h] -
            messages.data(oppose).values;
        assert(!v2f.data.values.hasNaN());
      }

      // factor -> var
      for (auto &f : messages.components<FHandleWrapper>()) {
        MGFHandle mgfh = f.topo.hd;
        FactorHandle fh = f.data.h;
        auto &vhs = _graph.topo(fh).lowers;

        const std::vector<MGV2FHandle> &v2fmsghs = orderedV2Fmsghs[mgfh];
        const std::vector<MGF2VHandle> &f2vmsghs = orderedF2Vmsghs[mgfh];
        assert(vhs.size() == v2fmsghs.size());
        assert(vhs.size() == f2vmsghs.size());

        // reset f2v messages
        for (auto &f2v : f2vmsghs) {
          messages.data(f2v)
              .values.setConstant(std::numeric_limits<double>::infinity());
        }

        // dispatch messages from this factor
        int fid = _graph.data(fh);
        auto &costFun = _factorCategories[fid].costs;
        auto &dims = orderedVDims[mgfh];
        assert(dims.size() > 0);
        std::vector<int> idx(dims.size(), 0);
        assert(vhs.size() == dims.size());
        while (true) {
          for (int i = 0; i < dims.size(); i++) {
            // i: output var
            // others: input _varCategories
            double &outValue = messages.data(f2vmsghs[i]).values[idx[i]];
            // theta
            double theta = costFun(idx.data(), idx.size(), fid, givenData);
            assert(!std::isnan(theta));
            assert(theta >= 0);

            // sum of other input _varCategories
            double sumOfOtherVars = 0.0;
            for (int j = 0; j < dims.size(); j++) {
              if (j == i)
                continue;
              sumOfOtherVars += messages.data(v2fmsghs[j]).values[idx[j]];
            }
            assert(!std::isnan(sumOfOtherVars));
            double score = theta + sumOfOtherVars;
            assert(!std::isnan(score));
            assert(score <= std::numeric_limits<double>::infinity());
            if (score < outValue) {
              outValue = score;
            }
          }

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
        }

        for (auto &f2v : f2vmsghs) {
          assert(messages.data(f2v).values.maxCoeff() <=
                 std::numeric_limits<double>::infinity());
          assert(!messages.data(f2v).values.hasNaN());
        }
      }
    }

    // marginalize on variables
    // reset marginals
    for (auto &v : _graph.elements<0>()) {
      varMarginals[v.topo.hd].setZero();
    }
    for (auto &f2v : messages.constraints<F2VMessage>()) {
      MGVHandle vh = f2v.topo.component<1>();
      varMarginals[messages.data(vh).h] += f2v.data.values;
      assert(varMarginals[messages.data(vh).h].maxCoeff() <
             std::numeric_limits<double>::infinity());
    }
    // get result
    for (auto &v : _graph.elements<0>()) {
      int label = -1;
      double curCost = std::numeric_limits<double>::infinity();
      const VectorXd &marginal = varMarginals[v.topo.hd];
      for (int i = 0; i < marginal.size(); i++) {
        assert(!core::IsInfOrNaN(marginal[i]));
        if (marginal[i] < curCost) {
          label = i;
          curCost = marginal[i];
        }
      }
      assert(label != -1);
      results[v.topo.hd] = label;
    }

    double e = energy(results);
    if (callback && !callback(epoch, e, e - lastE, results)) {
      break;
    }
    lastE = e;
  }

  return results;
}

FactorGraph::ResultTable
FactorGraph::solveWithSimpleCallback(int maxEpoch, int innerLoopNum,
                                     const SimpleCallbackFunction &callback,
                                     void *givenData) const {
  return solve(maxEpoch, innerLoopNum,
               [&callback](int epoch, double energy, double denergy,
                           const ResultTable &results) -> bool {
                 return callback(epoch, energy);
               },
               givenData);
}
}
}
