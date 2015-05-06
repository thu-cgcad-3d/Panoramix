#ifndef PANORAMIX_DISCRETIZATION_HPP
#define PANORAMIX_DISCRETIZATION_HPP

#include "../core/any.hpp"
#include "../core/basic_types.hpp"
#include "basic_types.hpp"
 
namespace panoramix {
    namespace gui {
 
        // discretize options
        using EntityPtr = core::AnyPtr;
        struct DiscretizeOptions {
            inline DiscretizeOptions()
            : color(0, 0, 0, 1), isolatedTriangles(false) {
                subdivisionNums[0] = 32;
                subdivisionNums[1] = 64;
            }
            EntityPtr entity;
            Color color;
            ColorTable colorTable;
            bool isolatedTriangles;
            int subdivisionNums[2];
        };


        // triangular mesh
        struct TriMesh {

            struct Vertex {
                Vertex();
                core::Vec4f position;
                core::Vec3f normal;
                core::Vec4f color; // the intrinsic color
                core::Vec2f texCoord;

                template <class Archive>
                inline void serialize(Archive & ar) {
                    ar(position, normal, color, texCoord, isSelected);
                }
            };

            using VertHandle = uint32_t;
            using PointHandle = uint32_t;
            using LineHandle = uint32_t;
            using TriangleHandle = uint32_t;

            std::vector<Vertex> vertices;
            std::vector<VertHandle> iverts;

            std::vector<PointHandle> iPoints;
            std::vector<LineHandle> iLines;
            std::vector<TriangleHandle> iTriangles;            

            std::vector<EntityPtr> entPoints;
            std::vector<EntityPtr> entLines;
            std::vector<EntityPtr> entTriangles;


            const Vertex & vertex(VertHandle vh) const { return vertices[vh]; }
            Vertex & vertex(VertHandle vh) { return vertices[vh]; }

            VertHandle addVertex(const Vertex & v, bool asPoint = false, EntityPtr ent = nullptr);
            VertHandle addVertex(const core::Point3 & p, const DiscretizeOptions & o, bool asPoint = false);
            VertHandle addVertex(const core::Point3 & p, const core::Vec3 & normal, const DiscretizeOptions & o, bool asPoint = false);

            size_t numerOfPoints() const;
            void fetchPointVerts(PointHandle p, VertHandle & v) const;

            template <class T> 
            T & entityAtPoint(PointHandle p) const { return *static_cast<T*>(entPoints[p]); }

            LineHandle addLine(VertHandle v1, VertHandle v2, EntityPtr ent = nullptr);
            LineHandle addIsolatedLine(const Vertex & v1, const Vertex & v2, EntityPtr ent = nullptr);
            size_t numberOfLines() const;
            void fetchLineVerts(LineHandle l, VertHandle & v1, VertHandle & v2) const;

            template <class T>
            T & entityAtLine(LineHandle l) const { return *static_cast<T*>(entLines[l]); }

            TriangleHandle addTriangle(VertHandle v1, VertHandle v2, VertHandle v3, EntityPtr entId = nullptr);
            TriangleHandle addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3, EntityPtr ent = nullptr);
            size_t numberOfTriangles() const;
            void fetchTriangleVerts(TriangleHandle t, VertHandle & v1, VertHandle & v2, VertHandle & v3) const;

            template <class T>
            T & entityAtTriangle(TriangleHandle t) const { return *static_cast<T*>(entTriangles[t]); }

            void addQuad(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4, EntityPtr ent = nullptr);
            void addPolygon(const std::vector<VertHandle> & vhs, EntityPtr ent = nullptr);

            void clear();

            core::Box3 boundingBox() const;

            template <class Archive> inline void serialize(Archive & ar) {
                ar(vertices, iverts, iPoints, iLines, iTriangles, entPoints, entLines, entTriangles);
            }

        };





        // discretization

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Point<T, 3> & p, const DiscretizeOptions & o){
            TriMesh::Vertex v;
            v.position = core::Vec4f(p[0], p[1], p[2], 1.0f);
            v.color = o.color;
            mesh.addVertex(v, true, o.entity);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Line<T, 3> & l, const DiscretizeOptions & o){
            TriMesh::Vertex v1, v2;
            v1.position = core::Concat(core::vec_cast<float>(l.first), 1.0f);
            v1.color = o.color;
            v2.position = core::Concat(core::vec_cast<float>(l.second), 1.0f);
            v2.color = o.color;
            mesh.addIsolatedLine(v1, v2, o.entity);
        }

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Chain<T, 3> & c, const DiscretizeOptions & o){
            if (c.size() == 0)
                return;
            std::vector<TriMesh::VertHandle> vhandles(c.size());
            for (int i = 0; i < c.size(); i++){
                TriMesh::Vertex v;
                v.position = core::Vec4f(c.points[i][0], c.points[i][1], c.points[i][2], 1.0);
                v.color = o.color;
                vhandles[i] = mesh.addVertex(v, false, o.entity);
            }
            for (int i = 0; i + 1 < c.size(); i++){
                mesh.addLine(vhandles[i], vhandles[i + 1], o.entity);
            }
            if (c.closed){
                mesh.addLine(vhandles.back(), vhandles.front(), o.entity);
            }
        }

        void Discretize(TriMesh & mesh, const core::LayeredShape3 & m, const DiscretizeOptions & o);
        void Discretize(TriMesh & mesh, const core::Sphere3 & s, const DiscretizeOptions & o);

        template <class T>
        inline void Discretize(TriMesh & mesh, const core::Sphere<T, 3> & s, const DiscretizeOptions & o){
            Discretize(mesh, core::Sphere3{ vec_cast<double>(s.center), static_cast<double>(s.radius) }, o);
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
                vhandles[i] = mesh.addVertex(v, false, o.entity);
            }
            mesh.addPolygon(vhandles, o.entity);
        }

        void Discretize(TriMesh & mesh, const SpatialProjectedPolygon & spp, const DiscretizeOptions & o);


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

        template <class T, class D>
        inline void Discretize(TriMesh & mesh, const core::Decorated<T, D> & d, const DiscretizeOptions & o){
            Discretize(mesh, d.component, o);
        }

        template <class IteratorT>
        inline void DiscretizeRange(TriMesh & mesh, IteratorT begin, IteratorT end, const DiscretizeOptions & o){
            auto oo = o;
            while(begin != end){
                oo.entity = &(*begin);
                Discretize(mesh, *begin, oo);
                ++begin;
            }
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