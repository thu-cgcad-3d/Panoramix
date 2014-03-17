#ifndef PANORAMIX_CORE_MESH_MAKER_HPP
#define PANORAMIX_CORE_MESH_MAKER_HPP

#include <cmath>

#include <Eigen/Core>

#include "mesh.hpp"
 
namespace panoramix {
	namespace core {

		template <class VertDataT, class VectorT>
		struct VertConstructorMaker {
			inline VertDataT operator () (const VectorT & v) const {
				return VertDataT(v);
			}
		};
 
		template <class VertDataT, class HalfDataT, class FaceDataT, class VertMakerT = VertConstructorMaker<VertDataT, Eigen::Vector3f>>
		void MakeTetrahedron(Mesh<VertDataT, HalfDataT, FaceDataT> & mesh, VertMakerT vmt = VertMakerT()) {
			using Eigen::Vector3f;
			mesh.clear();
			auto v1 = mesh.addVertex(vmt(Vector3f(0, 0, 0)));
			auto v2 = mesh.addVertex(vmt(Vector3f(0, 0, 1)));
			auto v3 = mesh.addVertex(vmt(Vector3f(0, 1, 0)));
			auto v4 = mesh.addVertex(vmt(Vector3f(1, 0, 0)));

			mesh.addFace({ v1, v2, v3 });
			mesh.addFace({ v1, v4, v2 });
			mesh.addFace({ v1, v3, v4 });
			mesh.addFace({ v2, v4, v3 });
		}
        
        template <class VertDataT, class HalfDataT, class FaceDataT, class VertMakerT = VertConstructorMaker<VertDataT, Eigen::Vector3f>>
        void MakeQuadFacedCube(Mesh<VertDataT, HalfDataT, FaceDataT> & mesh, VertMakerT vmt = VertMakerT()) {
            /*
             * 4 ----- 5
              /		  /|
             0 ----- 1 |
             |	     | |
             | 7	 | 6  -- x
             |	     |/
             3 ----- 2
             /
             y
             */
            using Eigen::Vector3f;
            
            mesh.clear();
            auto v1 = mesh.addVertex(vmt(Vector3f(0, 1, 1)));
            auto v2 = mesh.addVertex(vmt(Vector3f(1, 1, 1)));
            auto v3 = mesh.addVertex(vmt(Vector3f(1, 1, 0)));
            auto v4 = mesh.addVertex(vmt(Vector3f(0, 1, 0)));
            
            auto v5 = mesh.addVertex(vmt(Vector3f(0, 0, 1)));
            auto v6 = mesh.addVertex(vmt(Vector3f(1, 0, 1)));
            auto v7 = mesh.addVertex(vmt(Vector3f(1, 0, 0)));
            auto v8 = mesh.addVertex(vmt(Vector3f(0, 0, 0)));
            
            mesh.addFace({v1, v2, v3, v4});
            mesh.addFace({v2, v6, v7, v3});
            mesh.addFace({v6, v5, v8, v7});
            mesh.addFace({v5, v1, v4, v8});
            mesh.addFace({v5, v6, v2, v1});
            mesh.addFace({v4, v3, v7, v8});
            
        }


 
	}
}
 
#endif