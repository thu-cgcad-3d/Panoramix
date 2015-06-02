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

        using namespace panoramix::core;

        class ProjectiveAnchor;
        class ProjectiveComponent;

        class ProjectiveAnchor {
        public:
            using Ptr = std::unique_ptr<ProjectiveAnchor>;            
            virtual Vec3 direction() const = 0;
            virtual bool isGeneratedFrom(const ProjectiveComponent *) const { return false; }
        };

        class ProjectiveComponent {
        public:
            using Ptr = std::unique_ptr<ProjectiveComponent>;
            using CoeffVec = std::vector<double>;
            explicit ProjectiveComponent(int nparams) : _nparams(nparams) {}
            virtual CoeffVec coefficients(const Vec3 & direction) const = 0;
            virtual CoeffVec coefficients(const ProjectiveAnchor::Ptr & anchor) const {
                return coefficients(anchor->direction());
            }
            virtual void updateUsingParams(const double * params) const = 0;
            int nparams() const { return _nparams; }

        private:
            int _nparams;
        };


        class ProjectiveSolver {
        public:
            ProjectiveSolver() : _nvars(0), _neqs(0) {}

            int nvariables() const { return _nvars; }
            int nequations() const { return _neqs; }

            int bindLineDoF1(Line3 & l);
            int bindLineDoF2(Line3 & l);
            int bindPlaneDoF1(Plane3 & p);
            int bindPlaneDoF2(Plane3 & p, const Vec3 & axisDirection);
            int bindPlaneDoF3(Plane3 & p);
            int bindDoF1(Mesh<Point3> & mesh);

            int makeNormalAnchor(const Vec3 & d);
            int makeMeshVertexAnchor(int meshComponentId, Mesh<Point3>::VertHandle vh);

            int makeAEqualToBAt(int a, int b, int anchorId);
            int makeACloserThanBAt(int a, int b, int anchorId);
            int makeAEqualToDepthAt(int a, double d, int anchorId);
            int makeACloserThanDepthAt(int a, double d, int anchorId);
            int makeAFartherThanDepthAt(int a, double d, int anchorId);

            void solve() const;

        private:
            int append(ProjectiveComponent::Ptr && p);
            int append(ProjectiveAnchor::Ptr && a);

        private:
            std::vector<ProjectiveAnchor::Ptr> _anchors;
            std::vector<ProjectiveComponent::Ptr> _components;
            std::vector<int> _compStartPositionsInX;
            std::vector<SparseMatElementd> _Amat;
            std::vector<double> _Bmat;
            enum Op { Equal, LowerThan };
            std::vector<Op> _ops;
            int _nvars;
            int _neqs;
        };

    }
}

#endif