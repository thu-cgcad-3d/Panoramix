#include "../core/misc.hpp"

#include "renderable_object_tree.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        RenderableObjectTree::RenderableObjectTree() : _root(nullptr), _rtree(RenderableObjectBoundingBoxFunctor(this)) {}

        RenderableObjectTree::RenderableObjectTree(std::shared_ptr<RenderableObject> root) : _rtree(RenderableObjectBoundingBoxFunctor(this)) {
            installFromRoot(root);
        }

        void RenderableObjectTree::installFromRoot(std::shared_ptr<RenderableObject> root) {
            if (_root) {
                deleteAll();
            }
            _root = root;
            assert(_root);
            _root->depthFirstTraverse([this](RenderableObject * ro) -> bool {
                _objects.push_back(ro);
                _objectIds[ro] = _objects.size() - 1;
                // update model matrix
                if (!ro->parent()) { // is root
                    assert(!core::Contains(_calculatedModelMatrices, ro));
                    _calculatedModelMatrices[ro] = ro->modelMatrix();
                } else {
                    RenderableObject * parent = ro->parent();
                    assert(core::Contains(_calculatedModelMatrices, parent));
                    _calculatedModelMatrices[ro] = _calculatedModelMatrices[parent] * ro->modelMatrix(); // TODO?????
                }
                // update bounding box in world space
                Box3 primaryBbox = ro->primaryBoundingBox();
                Vec4 corner; corner[3] = 1;
                Box3 bbox;
                for (int i = 0; i < 2; i++) {
                    corner[0] = i == 0 ? primaryBbox.minCorner[0] : primaryBbox.maxCorner[0];
                    for (int j = 0; j < 2; j++) {
                        corner[1] = i == 0 ? primaryBbox.minCorner[1] : primaryBbox.maxCorner[1];
                        for (int k = 0; k < 2; k++) {
                            corner[2] = i == 0 ? primaryBbox.minCorner[2] : primaryBbox.maxCorner[2];
                            // transform this corner
                            Point4 newHCorner = _calculatedModelMatrices[ro] * corner;
                            bbox |= BoundingBox(Point3(newHCorner[0], newHCorner[1], newHCorner[2]) / newHCorner[3]);
                        }
                    }
                }
                _calculatedBoundingBoxes[ro] = bbox;
                // insert into rtree
                _rtree.insert(ro);
                return true;
            });
        }

        void RenderableObjectTree::deleteAll() {
            _root = nullptr;
            _objects.clear();
            _objectIds.clear();
            _calculatedBoundingBoxes.clear();
            _rtree.clear();
        }

        Box3 RenderableObjectTree::boundingBox() const {
            Box3 bbox;
            for (auto & b : _calculatedBoundingBoxes) {
                bbox |= b.second;
            }
            return bbox;
        }

        void RenderableObjectTree::renderWithCamera(RenderModeFlags mode, const core::PerspectiveCamera & cam) const {
            assert(_root);
            _root->depthFirstTraverse([this, &mode, &cam](RenderableObject * ro) -> bool {
                ro->render(mode, cam.viewProjectionMatrix() * _calculatedModelMatrices.at(ro));
            });
        }
    }
}