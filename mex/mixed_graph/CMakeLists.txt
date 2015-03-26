file (GLOB SOURCES "." *.cpp *.hpp)
panoramix_add_mex (mexMixedGraph MixedGraph.m ${SOURCES})