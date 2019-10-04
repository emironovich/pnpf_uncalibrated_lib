cmake_minimum_required (VERSION 3.0)
project (datagen)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# include_directories(/home/elizaveta/Matlab/extern/include)

set(SOURCE_LIB datagen.cpp)

add_library(datagen STATIC ${SOURCE_LIB})

target_link_libraries(datagen Eigen3::Eigen)	