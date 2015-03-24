#ifndef PANORAMIX_DISCRETIZATION_HPP
#define PANORAMIX_DISCRETIZATION_HPP

#include "../core/basic_types.hpp"
#include "basic_types.hpp"
 
namespace panoramix {
    namespace gui {
 

        // triangular mesh
        struct TriMesh {

            struct Vertex {
                Vertex();
                core::Vec4f position;
                core::Vec3f normal;
                core::Vec4f color; // the intrinsic color
                core::Vec2f texCoord;
                int entityIndex; // index in container
                uint8_t isSelected;
            };

            using VertHandle = uint32_t;
            using LineHandle = uint32_t;
            using TriangleHandle = uint32_t;

            std::vector<Vertex> vertices;
            std::vector<VertHandle> iPoints;
            std::vector<LineHandle> iLines;
            std::vector<TriangleHandle> iTriangles;


            VertHandle addVertex(const Vertex & v);

            LineHandle addLine(VertHandle v1, VertHandle v2);
            LineHandle addIsolatedLine(const Vertex & v1, const Vertex & v2);
            size_t numberOfLines() const;
            void fetchLineVerts(LineHandle l, VertHandle & v1, VertHandle & v2) const;

            TriangleHandle addTriangle(VertHandle v1, VertHandle v2, VertHandle v3);
            TriangleHandle addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3);
            size_t numberOfTriangles() const;
            void fetchTriangleVerts(TriangleHandle t, VertHandle & v1, VertHandle & v2, VertHandle & v3) const;

            void addQuad(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4);
            void addPolygon(const std::vector<VertHandle> & vhs);

            void clear();

            core::Box3 boundingBox() const;

            template <class Archive> inline void serialize(Archive & ar) {
                ar(vertices, iPoints, iLines, iTriangles);
            }

        };





        // discretization
        struct DiscretizeOptions {
            inline DiscretizeOptions()
            : color(0, 0, 0, 1), index(0), isolatedTriangles(false) {
                subdivisionNums[0] = 32;
                subdivisionNums[1] = 64;
            }
            Color color;
            ColorTable colorTable;
            int index;
            bool isolatedTriangles;
            int subdivisionNums[2];
        };

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Point<T, 3> & p, const DiscretizeOptions & o){
            TriMesh::Vertex v;
            v.position = core::Vec4f(p[0], p[1], p[2], 1.0f);
            v.color = o.color;
            v.entityIndex = o.index;
            mesh.addVertex(v);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Line<T, 3> & l, const DiscretizeOptions & o){
            TriMesh::Vertex v1, v2;
            v1.position = core::Concat(core::vec_cast<float>(l.first), 1.0f);
            v1.color = o.color;
            v1.entityIndex = o.index;
            v2.position = core::Concat(core::vec_cast<float>(l.second), 1.0f);
            v2.color = o.color;
            v2.entityIndex = o.index;
            mesh.addIsolatedLine(v1, v2);
        }

        void Discretize(TriMesh & mesh, const core::Sphere3 & s, const DiscretizeOptions & o);

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Sphere<T, 3> & s, const DiscretizeOptions & o){
            Discretize(mesh, core::Sphere3{ vec_cast<double>(s.center), static_cast<double>(s.radius) });
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Polygon<T, 3> & p, const DiscretizeOptions & o) {
            std::vector<TriMesh::VertHandle> vhandles(p.corners.size());
            for (int i = 0; i < p.corners.size(); i++){
                TriMesh::Vertex v;
                v.position = core::Vec4f(p.corners[i][0], p.corners[i][1], p.corners[i][2], 1.0);
                v.normal = core::vec_cast<float>(p.normal);
                v.color = o.color;
                v.entityIndex = o.index;
                vhandles[i] = mesh.addVertex(v);
            }
            mesh.addPolygon(vhandles);
        }

        void Discretize(TriMesh & mesh, const SpatialProjectedPolygon & spp, const DiscretizeOptions & o);

        template <class T, class AllocT>
        inline void Discretize(TriMesh & mesh, const std::vector<T, AllocT> & v, const DiscretizeOptions & o){
            auto oo = o;
            for (auto & e : v){
                Discretize(mesh, e, oo);
                oo.index++;
            }
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Classified<T> & c, const DiscretizeOptions & o){
            auto oo = o;
            oo.color = o.colorTable[c.claz];
            Discretize(mesh, c.component, oo);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const Colored<T> & c, const DiscretizeOptions & o){
            auto oo = o;
            oo.color = c.color;
            Discretize(mesh, c.component, oo);
        }

        template <class AttachedT, class T>
        inline void Discretize(TriMesh & mesh, const std::pair<AttachedT, T> & p, const DiscretizeOptions & o){
            Discretize(mesh, p.second, o);
        }


        // Is discretizable ?
        namespace {
            template <class T>
            struct IsDiscretizableImp {
                template <class TT>
                static auto test(int) -> decltype(
                    gui::Discretize(std::declval<TriMesh &>(), std::declval<TT>(), std::declval<DiscretizeOptions>()),
                    std::true_type()
                    );
                template <class>
                static std::false_type test(...);
                static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
            };
        }

        template <class T>
        struct IsDiscretizable : std::integral_constant<bool, IsDiscretizableImp<T>::value> {};


 
    }



    namespace core {

        inline Box3 BoundingBox(const gui::TriMesh & m) {
            return m.boundingBox();
        }

    }

}
 
#endif