#ifndef PANORAMIX_VIS_RENDERABLE_OBJECT_TREE_HPP
#define PANORAMIX_VIS_RENDERABLE_OBJECT_TREE_HPP

#include "../core/utilities.hpp"
#include "renderable_object.hpp"
 
namespace panoramix {
    namespace vis {

        // renderable object tree
        class RenderableObjectTree {
            struct RenderableObjectBoundingBoxFunctor {
                inline RenderableObjectBoundingBoxFunctor(RenderableObjectTree * const t) : tree(t) {}
                core::Box3 operator()(RenderableObject * ro) const { return tree->_calculatedBoundingBoxes.at(ro); }
                RenderableObjectTree * const tree;
            };

        public:
            RenderableObjectTree();
            explicit RenderableObjectTree(std::shared_ptr<RenderableObject> root);

            void installFromRoot(std::shared_ptr<RenderableObject> root);
            void deleteAll();

            Box3 boundingBox() const;
            void renderWithCamera(RenderModeFlags mode, const core::PerspectiveCamera & cam) const;

           /* std::vector<RenderableObject*> pickByRay(const core::InfiniteLine3 & ray, double distance) const;
            std::vector<RenderableObject*> pickByPointOnScreen(const QPointF & point, 
                const core::PerspectiveCamera & camera, double distanceOnScreen) const;*/

        private:
            std::shared_ptr<RenderableObject> _root;
            std::vector<RenderableObject*> _objects;
            std::map<RenderableObject*, int> _objectIds;
            std::map<RenderableObject*, core::Mat4> _calculatedModelMatrices;
            std::map<RenderableObject*, core::Box3> _calculatedBoundingBoxes;
            core::RTreeWrapper<RenderableObject*, RenderableObjectBoundingBoxFunctor> _rtree;
        };

    }

    namespace core {
        inline Box3 BoundingBox(const vis::RenderableObjectTree & t) {
            return t.boundingBox();
        }
    }

}
 
#endif