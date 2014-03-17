#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "../src/core/mesh.hpp"
#include "../src/core/mesh_maker.hpp"

#include "gtest/gtest.h"

using namespace panoramix;
using TestMesh = core::Mesh<Eigen::Vector3f>;

TEST(MeshTest, Tetrahedron) {

	TestMesh mesh;
	core::MakeTetrahedron(mesh);
	EXPECT_EQ(4, mesh.internalVertices().size());
	EXPECT_EQ(12, mesh.internalHalfEdges().size());
	EXPECT_EQ(4, mesh.internalFaces().size());
    
    for (int i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(3, nmesh.internalVertices().size());
        EXPECT_EQ(6, nmesh.internalHalfEdges().size());
        EXPECT_EQ(1, nmesh.internalFaces().size());
    }
    
}

TEST(MeshTest, Cube) {
    
    TestMesh mesh;
    core::MakeQuadFacedCube(mesh);
    EXPECT_EQ(8, mesh.internalVertices().size());
	EXPECT_EQ(24, mesh.internalHalfEdges().size());
	EXPECT_EQ(6, mesh.internalFaces().size());
    
    for (int i = 0; i < mesh.internalVertices().size(); i++) {
        TestMesh nmesh = mesh;
        nmesh.remove(TestMesh::VertHandle(i));
        nmesh.gc();
        
        EXPECT_EQ(7, nmesh.internalVertices().size());
        EXPECT_EQ(18, nmesh.internalHalfEdges().size());
        EXPECT_EQ(3, nmesh.internalFaces().size());
    }
    
}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
	return 0;
}