file (GLOB SOURCES "." *.cpp *.hpp)
panoramix_add_mex (mexBoundaryGraph BoundaryGraph.m ${SOURCES})
