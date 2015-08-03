#pragma once


#include "../core/basic_types.hpp"
#include "../core/utility.hpp"
#include "../core/containers.hpp"
#include "../core/single_view.hpp"
#include "../ml/factor_graph.hpp"
#include "rlgraph.hpp"

namespace pano {
    namespace experimental {

        using namespace pano::core; 

        // rlgraph path
        namespace {
            template <class HandleT>
            using RLGraphPatchContainer = std::unordered_set<HandleT>;
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

        struct RLGraphPatch {
            using ContentType = MetaBind<RLGraphPatchContainer,
                RHandle, LHandle, RRHandle, LLHandle, RLHandle, RRLHandle>;
            ContentType data;
            RHandle centerRh;

            template <class HandleT> void insert(HandleT h) { data.insert(h); }
            template <class HandleT> bool contains(HandleT h) const { return data.contains(h); }
            template <class HandleT> typename ContentType::EntryProperty<HandleT>::ContainerType & container() { return data.container<HandleT>(); }
            template <class HandleT> const typename ContentType::EntryProperty<HandleT>::ContainerType & container() const { return data.container<HandleT>(); }
            size_t size() const { return data.size(); }

            enum HandleType { RH, LH, RRH, LLH };
            std::pair<std::vector<HandleType>, std::vector<int>> handleInfo() const;
            RLGraphPatchDict<int> handlePositions() const;

            template <class Archiver>
            void serialize(Archiver & ar) { ar(data, centerRh); }
        };

        RLGraphPatch SamplePatch(const RLGraph & graph, RHandle centerRh, int maxRNum, double angleRadius);





        // rlgraph handles -> factor graph var handles
        using RLGraphVarHandleTable = Table<ml::FactorGraph::VarHandle,
            RHandle, LHandle, RRHandle, LLHandle>;
        using RLGraphLabel = int;
        using RLGraphAllowedLabelsTable = Table<std::vector<RLGraphLabel>,
            RHandle, LHandle, RRHandle, LLHandle>;


        using RLFeatureVec = std::vector<double>;
        using RLFeatureDictTable = Table<Dictionary<RLFeatureVec>,
            RHandle, LHandle, RRHandle, LLHandle, RRRHandle, RRRRHandle>;
        using RLPatchFeatureData =
            std::tuple<RLFeatureVec, std::unordered_map<RHandle, Plane3>, std::unordered_map<LHandle, Line3>>;
        

        void Visualize(const RLGraph & g, const PanoramicView & view,
            const std::unordered_map<RHandle, Plane3> & planes, 
            const std::unordered_map<LHandle, Line3> & lines);





        struct RLGraphLabelsProvider {
            size_t vpsNum;
            
            RLGraphLabel rFree() const { return 0; }
            RLGraphLabel rTowardVP(int vpid) const { assert(vpid >= 0 && vpid < vpsNum); return vpid + 1; }
            RLGraphLabel rVertical() const { return vpsNum + 1; }
            RLGraphLabel rNotPlanar() const { return vpsNum + 2; }
            size_t rNum() const { return vpsNum + 3; }

            bool rIsToVP(RLGraphLabel label) const { return label > 0 && label <= vpsNum; }
            int rToVPId(RLGraphLabel label) const { assert(rIsToVP(label)); return label - 1; }

            RLGraphLabel lFree() const { return 0; }
            RLGraphLabel lTowardVP(int vpid) const { assert(vpid >= 0 && vpid < vpsNum); return vpid + 1; }
            size_t lNum() const { return vpsNum + 1; }

            bool lIsToVP(RLGraphLabel label) const { return label > 0 && label <= vpsNum; }
            int lToVPId(RLGraphLabel label) const { assert(rIsToVP(label)); return label - 1; }

            RLGraphLabel rrConnected() const { return 0; }
            RLGraphLabel rrFirstIsCloser() const { return 1; }
            RLGraphLabel rrSecondIsCloser() const { return 2; }
            RLGraphLabel rrDisconnected() const { return 3; }
            size_t rrNum() const { return 4; }

            RLGraphLabel rrReverse(RLGraphLabel label) const { return ((3 - (label)) % 3); }

            RLGraphLabel llConnected() const { return 0; }
            RLGraphLabel llDisconnected() const { return 1; }
            size_t llNum() const { return 2; }
        };


       
        
        RLGraphPatch MakeASinglePatch(const RLGraph & g);
        void FillSingleLabels(const RLGraph & g, const std::vector<Vec3> & vps,
            const RLGraphPatch & patch,
            RLGraphPatchDict<RLGraphLabel> & labels);
        bool Reconstruct(const RLGraph & g, const std::vector<Vec3> & vps,
            const RLGraphPatch & patch,
            const RLGraphPatchDict<RLGraphLabel> & labels,
            std::unordered_map<RHandle, Plane3> & planes,
            std::unordered_map<LHandle, Line3> & lines);




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
                _labelsProvider = RLGraphLabelsProvider{ vps.size() };
            }

            void preprocess();
            std::vector<double> suggestedWeights() const;
            void inference(const std::vector<double> & weights,
                HandledTable<RHandle, Plane3> & planes, HandledTable<LHandle, Line3> & lines) const;
            void inferenceWithLossFunction(const std::vector<double> & weights,
                std::function<double(const HandledTable<RHandle, Plane3> & planes, 
                const HandledTable<LHandle, Line3> & lines)> & lossFun) const;

            const RLFeatureDictTable & graphFeatureDictTable() const { return _graphFeatureDict; }
            RLFeatureDictTable & graphFeatureDictTable() { return _graphFeatureDict; }
            const std::vector<Dictionary<RLPatchFeatureData>> & patchFeatureDictTable() const { return _patchFeatureDict; }
            std::vector<Dictionary<RLPatchFeatureData>> & patchFeatureDictTable() { return _patchFeatureDict; }

        private:
            void determinFactorGraphVars();
            void determinFactorGraphFactors();

            int featureLength() const;
            Failable<RLPatchFeatureData> featureInReconstructedPatch(const RLGraphPatch & patch, 
                const std::vector<RLGraphLabel> & labels) const;

            void findPeakyRhs(double rangeAngle = DegreesToRadians(6));
            void findHorizonRhs(double rangeAngle = DegreesToRadians(10));
            void appliedRLCons(const RLGraphPatch & patch, const RLGraphPatchDict<RLGraphLabel> & patchLabels,
                std::vector<std::pair<RHandle, LHandle>> & rlcons,
                std::vector<const std::vector<Vec3>*> & nanchorsPtrTable) const;

            const Vec3 & up() const { return _vps[0]; }
            Vec3 down() const { return -_vps[0]; }

            const RLGraphLabelsProvider & provideLabels() const { return _labelsProvider; }

        private:
            RLGraph _g;
            std::vector<Vec3> _vps; // _vps[0] MUST BE UP
            RLGraphLabelsProvider _labelsProvider;
            
            ml::FactorGraph _fg;
            RLGraphVarHandleTable _varhs;
            RLGraphAllowedLabelsTable _allowedLabels;

            double _rAreaSum;
            double _lLengthSum;
            double _rrLengthSum;
            std::vector<std::unordered_set<RHandle>> _peakyRhs;
            std::unordered_set<RHandle> _horizonRhs;

            // feature table for each factor
            RLFeatureDictTable _graphFeatureDict;

            std::vector<RLGraphPatch> _patches;
            std::vector<Dictionary<RLPatchFeatureData>> _patchFeatureDict; // [patch feature size x [\Pi {related var dims}]


        };



    }
}
