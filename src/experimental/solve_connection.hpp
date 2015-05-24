#ifndef PANORAMIX_EXPERIMENTAL_SOLVE_CONNECTION_HPP
#define PANORAMIX_EXPERIMENTAL_SOLVE_CONNECTION_HPP

#include "../core/basic_types.hpp"
#include "../core/utility.hpp"
#include "../core/generic_topo.hpp"
#include "../core/cons_graph.hpp"
#include "../core/cameras.hpp"
#include "../core/generic_topo.hpp"

namespace panoramix {
    namespace experimental {

        using namespace core;

        class ProjectiveComponent {
        public:
            explicit ProjectiveComponent(int nparams) : _nparams(nparams) {}
            ProjectiveComponent(const ProjectiveComponent &) = delete;
            ProjectiveComponent & operator = (const ProjectiveComponent &) = delete;

            virtual DenseMatd coefficients(const Vec3 & direction) const = 0;
            int nparams() const { return _nparams; }

        private:
            int _nparams;
        };

        // [1.0/norm(line.center())]
        class LineDoF1 : public ProjectiveComponent {
        public:
            LineDoF1(const Line3 & l);
            DenseMatd coefficients(const Vec3 & direction) const override;
            Line3 line;
        };

        // [1.0/norm(line.first), 1.0/norm(line.second)]
        class LineDoF2 : public ProjectiveComponent {
        public:
            LineDoF2(const Line3 & l);
            DenseMatd coefficients(const Vec3 & direction) const override;
            Line3 line;
        };

        // [1.0/norm(center)]
        class RegionDoF1 : public ProjectiveComponent {
        public:
            RegionDoF1(const Plane3 & p);
            DenseMatd coefficients(const Vec3 & direction) const override;
            Plane3 plane;
        };

        // [a, b] in ax+by+cz=1
        class RegionDoF2 : public ProjectiveComponent {
        public:
            RegionDoF2(const Vec3 & ad);
            DenseMatd coefficients(const Vec3 & direction) const override;
            Vec3 axisDirection;
        };

        // [a, b, c] in ax+by+cz=1
        class RegionDoF3 : public ProjectiveComponent {
        public:
            RegionDoF3();
            DenseMatd coefficients(const Vec3 & direction) const override;
        };

        // [scale of placed mesh / scale of real mesh]
        struct MeshDoF1Data;
        class MeshDoF1 : ProjectiveComponent{
        public:
            MeshDoF1(const Mesh<Point3> & m);
            std::vector<Point3> intersectionsOnPlacedMesh(const Vec3 & direction) const;
            DenseMatd coefficients(const Vec3 & direction) const override;
        private:
            std::unique_ptr<MeshDoF1Data> _data;
        };



        // constraints

        enum class DepthRelation {
            Connected, LeftIsCloser, RightIsCloser
        };

        struct Anchor {
            int objId;
            int anchorId;
            bool operator == (Anchor oa) const { 
                return std::tie(objId, anchorId) == std::tie(oa.objId, anchorId);
            }
            bool operator < (Anchor oa) const {
                return std::tie(objId, anchorId) < std::tie(oa.objId, anchorId);
            }
        };

        struct AnchorOrDepth{
            union {
                Anchor oa;
                double d;
            };
            bool isAnchor;
            AnchorOrDepth() : isAnchor(false), d(0.0) {}
            AnchorOrDepth(Anchor objanchor) : isAnchor(true) { oa = objanchor; }
            AnchorOrDepth(double dp) : isAnchor(false) { d = dp; }
            operator double() const { assert(!isAnchor); return d; }
            operator Anchor () const { assert(isAnchor); return oa; }
        };

        struct ProjectiveConstraint {
            DepthRelation relation;
            Vec3 direction;
            AnchorOrDepth left;
            AnchorOrDepth right;
        };


        using ProjectiveComponentArray = 
            std::vector<std::unique_ptr<ProjectiveComponent>>;
        using ProjectiveConstraintArray =
            std::vector<ProjectiveConstraint>;


        template <class SparseMatElementT>
        std::pair<int, int> PrepareData(
            const ProjectiveComponentArray & comps, 
            const ProjectiveConstraintArray & conss,
            std::vector<int> & compId2varStartPosition,
            std::vector<SparseMatElementT> & Aleft,
            std::vector<SparseMatElementT> & Aright,
            std::vector<double> & Bleft,
            std::vector<double> & Bright,
            std::vector<int> & signs){

            compId2varStartPosition = std::vector<int>(comps.size(), -1);
            int varnum = 0;
            for (int i = 0; i < comps.size(); i++){
                compId2varStartPosition[i] = varnum;
                varnum += comps[i]->nparams();
            }
            
            int consnum = conss.resize();
            Aleft.reserve(consnum * 3);
            Aright.reserve(consnum * 3);

            Bleft = std::vector<double>(consnum, 0.0);
            Bright = std::vector<double>(consnum, 0.0);
            signs = std::vector<int>(consnum, 0);

            for (int i = 0; i < conss.size(); i++){
                auto & cons = conss[i];
                auto & dir = cons.direction;
                if (cons.left.isAnchor){
                    Anchor anchor = cons.left;
                    auto coeffs = comps[anchor.objId]->coefficients(dir).at(anchor.anchorId);
                    int compStartPos = compId2varStartPosition.at(anchor.objId);
                    for (int j = 0; j < coeffs.size(); j++){
                        Aleft.emplace_back(i, compStartPos + j, coeffs[j]);
                    }
                }
                else{
                    Bleft[i] = cons.left;
                }
                if (cons.right.isAnchor){
                    Anchor anchor = cons.right;
                    auto coeffs = comps[anchor.objId]->coefficients(dir).at(anchor.anchorId);
                    int compStartPos = compId2varStartPosition.at(anchor.objId);
                    for (int j = 0; j < coeffs.size(); j++){
                        Aright.emplace_back(i, compStartPos + j, coeffs[j]);
                    }
                }
                else{
                    Bright[i] = cons.right;
                }
                switch (cons.relation) {
                case DepthRelation::Connected: signs[i] = 0; break;
                case DepthRelation::LeftIsCloser: signs[i] = 1; break;
                case DepthRelation::RightIsCloser: signs[i] = 2 : break;
                default:
                    break;
                }
            }

            return std::make_pair(varnum, consnum);
        }



    }
}

#endif