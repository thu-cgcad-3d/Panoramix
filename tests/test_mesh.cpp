#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "../src/core/mesh.hpp"
#include "../src/core/mesh_maker.hpp"

#include "gtest/gtest.h"

using namespace panoramix;
using TestMesh = core::Mesh<Eigen::Vector3f>;

TEST(MeshTest, Conversion) {
	using CVMesh = core::Mesh<cv::Vec3f>;
	CVMesh mesh;
	core::MakeTetrahedron(mesh);
	EXPECT_EQ(4, mesh.internalVertices().size());
	EXPECT_EQ(12, mesh.internalHalfEdges().size());
	EXPECT_EQ(4, mesh.internalFaces().size());
}

TEST(MeshTest, Tetrahedron) {

	TestMesh mesh;
	core::MakeTetrahedron(mesh);
	EXPECT_EQ(4, mesh.internalVertices().size());
	EXPECT_EQ(12, mesh.internalHalfEdges().size());
	EXPECT_EQ(4, mesh.internalFaces().size());
    
    for (size_t i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(3, nmesh.internalVertices().size());
        EXPECT_EQ(6, nmesh.internalHalfEdges().size());
        EXPECT_EQ(1, nmesh.internalFaces().size());
    }

	for (size_t i = 0; i < mesh.internalHalfEdges().size(); i++){
		TestMesh nmesh = mesh;
		nmesh.remove(TestMesh::HalfHandle(i));
		nmesh.gc();

		EXPECT_EQ(4, nmesh.internalVertices().size());
		EXPECT_EQ(10, nmesh.internalHalfEdges().size());
		EXPECT_EQ(2, nmesh.internalFaces().size());
	}

	for (size_t i = 0; i < mesh.internalFaces().size(); i++) {
		TestMesh nmesh = mesh;
		nmesh.remove(TestMesh::FaceHandle(i));
		nmesh.gc();

		EXPECT_EQ(4, nmesh.internalVertices().size());
		EXPECT_EQ(12, nmesh.internalHalfEdges().size());
		EXPECT_EQ(3, nmesh.internalFaces().size());
	}
    
}

TEST(MeshTest, Cube) {
    
    TestMesh mesh;
    core::MakeQuadFacedCube(mesh);
    EXPECT_EQ(8, mesh.internalVertices().size());
	EXPECT_EQ(24, mesh.internalHalfEdges().size());
	EXPECT_EQ(6, mesh.internalFaces().size());
    
    for (size_t i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(7, nmesh.internalVertices().size());
        EXPECT_EQ(18, nmesh.internalHalfEdges().size());
        EXPECT_EQ(3, nmesh.internalFaces().size());
    }

	for (size_t i = 0; i < mesh.internalHalfEdges().size(); i++){
		TestMesh nmesh = mesh;
		nmesh.remove(TestMesh::HalfHandle(i));
		nmesh.gc();

		EXPECT_EQ(8, nmesh.internalVertices().size());
		EXPECT_EQ(22, nmesh.internalHalfEdges().size());
		EXPECT_EQ(4, nmesh.internalFaces().size());
	}

	for (size_t i = 0; i < mesh.internalFaces().size(); i++) {
		TestMesh nmesh = mesh;
		nmesh.remove(TestMesh::FaceHandle(i));
		nmesh.gc();

		EXPECT_EQ(8, nmesh.internalVertices().size());
		EXPECT_EQ(24, nmesh.internalHalfEdges().size());
		EXPECT_EQ(5, nmesh.internalFaces().size());
	}
    
}

TEST(MeshTest, DISABLED_Sphere) {

	TestMesh mesh;
	core::MakeQuadFacedSphere(mesh, 10, 5);

}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}