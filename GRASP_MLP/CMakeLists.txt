cmake_minimum_required(VERSION 3.3)
project(GRASP_MLP)


set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS_DEBUG "-O3")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++17 -O3" )

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include_directories(src/linear_split)
include_directories(src/metaheuristic/methods)
include_directories(src/metaheuristic/model)

find_package(Boost 1.57.0)
include_directories(${boost_INCLUDE_DIRS})


if(Boost_FOUND)

    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")

    include_directories(${Boost_INCLUDE_DIRS})

endif()

FIND_PACKAGE(OpenMP)


add_executable(GRASP_MLP
        src/linear_split/Pb_Data.cpp
        src/linear_split/Pb_Data.h
        src/linear_split/Split.h
        src/linear_split/Split_Linear.cpp
        src/linear_split/Split_Linear.h
        src/linear_split/Trivial_Deque.h
        src/metaheuristic/methods/CircleSector.h
        src/metaheuristic/methods/cw_savings.cpp
        src/metaheuristic/methods/cw_savings.h
        src/metaheuristic/methods/input_parser.cpp
        src/metaheuristic/methods/input_parser.h
        src/metaheuristic/methods/giant_tour_creation.cpp
        src/metaheuristic/methods/giant_tour_creation.h
        src/metaheuristic/methods/neighborhood.cpp
        src/metaheuristic/methods/neighborhood.h
        src/metaheuristic/methods/vnd.cpp
        src/metaheuristic/methods/vnd.h
        src/metaheuristic/methods/preprocessing.h
        src/metaheuristic/methods/preprocessing.cpp
        src/metaheuristic/model/instance.cpp
        src/metaheuristic/model/instance.h
        src/metaheuristic/model/read_vrp.cpp
        src/metaheuristic/model/read_vrp.h
        src/metaheuristic/model/solution.cpp
        src/metaheuristic/model/solution.h
        src/metaheuristic/twiggy_core/criterion.cpp
        src/metaheuristic/twiggy_core/criterion.h
        src/metaheuristic/twiggy_core/decision_tree.cpp
        src/metaheuristic/twiggy_core/decision_tree.h
        src/metaheuristic/twiggy_core/random_forest.cpp
        src/metaheuristic/twiggy_core/random_forest.h
        src/metaheuristic/twiggy_core/splitter.cpp
        src/metaheuristic/twiggy_core/splitter.h
        src/metaheuristic/twiggy_core/utils.cpp
        src/metaheuristic/twiggy_core/utils.h
        src/metaheuristic/csv/parser.hpp
        src/grasp_mlp.cpp
        )


target_link_libraries(GRASP_MLP ${Boost_LIBRARIES})
target_link_libraries(GRASP_MLP OpenMP::OpenMP_CXX)
