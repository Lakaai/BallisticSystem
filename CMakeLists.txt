cmake_minimum_required(VERSION 3.24)
project(BallisticSystem CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)
find_package(autodiff 1.0 REQUIRED)

if(autodiff_FOUND)
    message(STATUS "Autodiff found!")
else()
    message(FATAL_ERROR "Autodiff NOT found.")
endif()

# Find Python3
# find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
# include_directories(${Python3_INCLUDE_DIRS})

# Find NumPy
# execute_process(
#     COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
#     OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# include_directories(${NUMPY_INCLUDE_DIR})


# Get from git clone https://github.com/lava/matplotlib-cpp.git
# Include matplotlib-cpp header
# include_directories("${CMAKE_SOURCE_DIR}/include/matplotlib-cpp")

#add_executable(BallisticSystem src/main.cpp src/System.cpp src/KalmanFilter.cpp)
set(PROJECT_SOURCES
    src/main.cpp
    src/System.cpp
    src/KalmanFilter.cpp
    src/StateEstimator.cpp
)

if(Eigen3_FOUND)
    message(STATUS "Found Eigen3 version: ${EIGEN3_VERSION_STRING}")
    message(STATUS "Eigen3 include directories: ${EIGEN3_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "Eigen3 not found. Please install Eigen3 and ensure it's in CMake's search path.")
endif()

get_target_property(AUTODIFF_INCLUDE_DIRS autodiff::autodiff INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Autodiff include directories: ${AUTODIFF_INCLUDE_DIRS}")

add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE 
    ${EIGEN3_INCLUDE_DIRS}
    ${AUTODIFF_INCLUDE_DIRS}
)
target_link_libraries(${PROJECT_NAME} PRIVATE autodiff::autodiff)

# target_include_directories(your_executable_or_library_name PRIVATE ${EIGEN3_INCLUDE_DIRS})
# target_link_libraries(BallisticSystem PRIVATE ${Python3_LIBRARIES})





