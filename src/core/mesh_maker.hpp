#ifndef PANORAMIX_CORE_MESH_MAKER_HPP
#define PANORAMIX_CORE_MESH_MAKER_HPP

#include <cmath>

#include <Eigen/Core>

#include "mesh.hpp"
 
namespace panoramix {
	namespace core {

		template <class VertDataT, int dim, class ValueT>
		struct VertConstructorMaker {
			inline VertDataT operator () (const Eigen::Matrix<ValueT, dim, 1> & v) const {
				return VertDataT(v);
			}
		};
 
		template <class VertDataT, class HalfDataT, class FaceDataT, class VertMakerT = VertConstructorMaker<VertDataT, 3, double>>
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

 
	}
}
 
#endif