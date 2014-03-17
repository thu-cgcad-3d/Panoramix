#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include "../src/core/mesh.hpp"
#include "../src/core/mesh_maker.hpp"

#include "gtest/gtest.h"

using namespace panoramix;
using TestMesh = core::Mesh<Eigen::Vector3f>;

TEST(MESH_TEST, CUBE) {

	TestMesh mesh;
	core::MakeTetrahedron(mesh);
	EXPECT_EQ(8, mesh.internalVertices().size());
	EXPECT_EQ(24, mesh.internalHalfEdges().size());
	EXPECT_EQ(6, mesh.internalFaces().size());

}

int main(int argc, char * argv[], char * envp[])
{
	for (int i = 0; i < argc; i++) {
		std::cout << "[INPUT]:" << argv[i] << std::endl;
	}
	char** env;
	for (env = envp; *env != 0; env++) {
		char* thisEnv = *env;
		std::cout << "[ENV]:" << thisEnv << std::endl;
	}
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
	return 0;
}