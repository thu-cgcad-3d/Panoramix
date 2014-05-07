#ifndef PANORAMIX_CORE_MESH_MAKER_HPP
#define PANORAMIX_CORE_MESH_MAKER_HPP

#include <cmath>

#include "mesh.hpp"
 
namespace panoramix {
    namespace deriv {

        template <class VertDataT, class InputValueT>
        struct Vert3MakerDefault {
            inline VertDataT operator () (const InputValueT & v1, const InputValueT & v2, const InputValueT & v3) const {
                return VertDataT(v1, v2, v3);
            }
        };
 
        template <class VertDataT, class HalfDataT, class FaceDataT, class Vert3MakerT = Vert3MakerDefault<VertDataT, float>>
        void MakeTetrahedron(Mesh<VertDataT, HalfDataT, FaceDataT> & mesh, Vert3MakerT vmt = Vert3MakerT()) {
            mesh.clear();
            auto v1 = mesh.addVertex(vmt(0, 0, 0));
            auto v2 = mesh.addVertex(vmt(0, 0, 1));
            auto v3 = mesh.addVertex(vmt(0, 1, 0));
            auto v4 = mesh.addVertex(vmt(1, 0, 0));

            mesh.addFace({ v1, v2, v3 });
            mesh.addFace({ v1, v4, v2 });
            mesh.addFace({ v1, v3, v4 });
            mesh.addFace({ v2, v4, v3 });
        }
        
        template <class VertDataT, class HalfDataT, class FaceDataT, class Vert3MakerT = Vert3MakerDefault<VertDataT, float>>
        void MakeQuadFacedCube(Mesh<VertDataT, HalfDataT, FaceDataT> & mesh, Vert3MakerT vmt = Vert3MakerT()) {
            /*
                   4 ----- 5
                  /		  /|
                 0 ----- 1 |
                 |	     | |
                 | 7	 | 6  -- x
                 |	     |/
                 3 ----- 2
                /
               y
             
             */ 
            mesh.clear();
            auto v1 = mesh.addVertex(vmt(0, 1, 1));
            auto v2 = mesh.addVertex(vmt(1, 1, 1));
            auto v3 = mesh.addVertex(vmt(1, 1, 0));
            auto v4 = mesh.addVertex(vmt(0, 1, 0));
            
            auto v5 = mesh.addVertex(vmt(0, 0, 1));
            auto v6 = mesh.addVertex(vmt(1, 0, 1));
            auto v7 = mesh.addVertex(vmt(1, 0, 0));
            auto v8 = mesh.addVertex(vmt(0, 0, 0));
            
            mesh.addFace({v1, v2, v3, v4});
            mesh.addFace({v2, v6, v7, v3});
            mesh.addFace({v6, v5, v8, v7});
            mesh.addFace({v5, v1, v4, v8});
            mesh.addFace({v5, v6, v2, v1});
            mesh.addFace({v4, v3, v7, v8});
            
        }

        template <class VertDataT, class HalfDataT, class FaceDataT, class Vert3MakerT = Vert3MakerDefault<VertDataT, double>>
        void MakeQuadFacedSphere(Mesh<VertDataT, HalfDataT, FaceDataT> & mesh, int m = 5, int n = 10, Vert3MakerT vmt = Vert3MakerT()) {
            using ThisMesh = Mesh<VertDataT, HalfDataT, FaceDataT>;
            using ThisVertHandle = typename ThisMesh::VertHandle;
            mesh.clear();
            mesh.internalVertices().reserve(m * n);
            mesh.internalHalfEdges().reserve(4 * m * n);
            mesh.internalFaces().reserve(m * n);
            std::vector<std::vector<ThisVertHandle>> vhs(m, std::vector<ThisVertHandle>(n-1));
            for (int i = 0; i < m; i++){
                for (int j = 0; j < n - 1; j++){
                    double xratio = 1.0f - 1.0f / (n - 1) * j;
                    double yratio = 1.0f / (m - 1) * i;
                    double xangle = M_PI * 2 * xratio;
                    double yangle = M_PI * yratio - M_PI_2;
                    vhs[i][j] = mesh.addVertex(vmt(sin(xangle - M_PI_2)*cos(yangle), cos(xangle - M_PI_2)*cos(yangle), sin(yangle)));
                }
            }
            for (int i = 1; i < m; i++){
                for (int j = 1; j < n - 1; j++){
                    mesh.addFace({ vhs[i][j], vhs[i][j - 1], vhs[i - 1][j - 1], vhs[i - 1][j] });
                }
                mesh.addFace({ vhs[i][0], vhs[i][n - 2], vhs[i - 1][n - 2], vhs[i - 1][0] });
            }
        }
 
    }
}
 
#endif