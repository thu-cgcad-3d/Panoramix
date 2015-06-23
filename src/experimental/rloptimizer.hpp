#ifndef PANORAMIX_EXPERIMENTAL_RLOPTIMIZER_HPP
#define PANORAMIX_EXPERIMENTAL_RLOPTIMIZER_HPP

#include "../core/basic_types.hpp"
#include "../core/utility.hpp"
#include "../core/single_view.hpp"
#include "../ml/factor_graph.hpp"
#include "rlgraph.hpp"

namespace panoramix {
    namespace experimental {

        using namespace panoramix::core; 

        // rlgraph path
        namespace {
            template <class HandleT>
            using RLGraphPatchContainer = std::unordered_set<HandleT>;
        }
        using RLGraphPatch = MetaBind<RLGraphPatchContainer,
            RHandle, LHandle, RRHandle, LLHandle, RLHandle, RRLHandle>;

        std::vector<RLGraphPatch> MakePatches(const RLGraph & graph, int n, double angleRadius);

        // rlgraph patch dict
        namespace {
            template <class T>
            struct RLGraphPatchDictStruct {
                template <class HandleT>
                using Container = std::unordered_map<HandleT, T>;
                using type = MetaBind<Container,
                    RHandle, LHandle, RRHandle, LLHandle>;
            };
        }
        template <class T>
        using RLGraphPatchDict = typename RLGraphPatchDictStruct<T>::type;


        // rlgraph handles -> factor graph var handles
        using RLGraphVarHandleTable = MixedHandledTable<ml::FactorGraph::VarHandle,
            RHandle, LHandle, RRHandle, LLHandle>;


        // rloptimizer
        class RLOptimizer {
        public:
            RLOptimizer() {}
            template <class RLGraphT, class VPsT, class PatchesT>
            RLOptimizer(RLGraphT && graph, VPsT && vps, PatchesT && cs) {
                setup(std::forward<RLGraphT>(graph), std::forward<PatchesT>(cs));
            }

            RLOptimizer(const RLOptimizer &) = delete;
            RLOptimizer & operator = (const RLOptimizer &) = delete;

        public:
            template <class RLGraphT, class VPsT, class PatchesT>
            void setup(RLGraphT && graph, VPsT && vps, PatchesT && cs) {
                _g = std::forward<RLGraphT>(graph);
                _vps = std::forward<VPsT>(vps);
                _patches = std::forward<PatchesT>(cs);
            }

            void preprocess();
            std::vector<double> suggestedWeights() const;
            void inference(const std::vector<double> & weights,
                HandledTable<RHandle, Plane3> & planes, HandledTable<LHandle, Line3> & lines) const;
            void inferenceWithLossFunction(const std::vector<double> & weights,
                std::function<double(const HandledTable<RHandle, Plane3> & planes, const HandledTable<LHandle, Line3> & lines)> & lossFun) const;

        private:
            int featureLength() const;

            std::vector<double> featureInR(RHandle rh, int rlabel) const;
            std::vector<double> featureInL(LHandle lh, int llabel) const;
            std::vector<double> featureInRR(RRHandle rrh, int rrlabel) const;
            std::vector<double> featureInLL(LLHandle llh, int lllabel) const;
            std::vector<double> featureInRRRorRRRR(const int * rrlabels, const int * rrascends, int rrnum) const;

            std::vector<double> featureInReconstructedPatch(
                int patchId,
                const int * allLabelsInPatch,
                std::unordered_map<RHandle, Plane3> & planes,
                std::unordered_map<LHandle, Line3> & lines) const;

            void findPeakyRhs(double rangeAngle = DegreesToRadians(6));
            void findHorizonRhs(double rangeAngle = DegreesToRadians(10));
            void countIncidenceAndIntersectionLLhs();
            void appliedRLCons(const RLGraphPatch & patch, const RLGraphPatchDict<int> & pathLabels,
                std::vector<std::pair<RHandle, LHandle>> & rlcons,
                std::vector<const std::vector<Vec3>*> & nanchorsPtrTable) const;

            const Vec3 & up() const { return _vps[0]; }
            Vec3 down() const { return -_vps[0]; }

        private:
            RLGraph _g;
            std::vector<Vec3> _vps; // _vps[0] MUST BE UP
            
            ml::FactorGraph _fg;
            RLGraphVarHandleTable _varhs;

            double _rAreaSum;
            double _lLengthSum;
            double _rrLengthSum;
            std::vector<std::unordered_set<RHandle>> _peakyRhs;
            std::unordered_set<RHandle> _horizonRhs;
            int _incidenceLLNum, _intersectionLLNum;

            template <class HandleT>
            using FeatureTable = HandledTable<HandleT, std::vector<std::vector<double>>>;
            FeatureTable<RHandle> _rFeatureTable;
            FeatureTable<LHandle> _lFeatureTable;
            FeatureTable<RRHandle> _rrFeatureTable;
            FeatureTable<LLHandle> _llFeatureTable;
            FeatureTable<RRRHandle> _rrrFeatureTable;
            FeatureTable<RRRRHandle> _rrrrFeatureTable;

            std::vector<RLGraphPatch> _patches;
            std::vector<std::vector<std::vector<double>>> _patchFeatureTable; // [patch feature size x [\Pi {related var dims}]
        };



    }
}


#endif