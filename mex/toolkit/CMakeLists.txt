file (GLOB SOURCES "." *.cpp *.hpp)
panoramix_add_mex (mexToolkit Toolkit.m ${SOURCES})