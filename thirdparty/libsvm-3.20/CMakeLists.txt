
set (SOURCE svm.cpp svm.h svm-predict.c svm-scale.c svm-train.c)

source_group (libsvm-3.20 FILES ${SOURCE})
add_library (libsvm-3.20 STATIC ${SOURCE})

set_property(TARGET libsvm-3.20 PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")