if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  SET(CMAKE_CXX_STANDARD 11)
endif()

cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

if(${CMAKE_VERSION} VERSION_LESS "3.1.0")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    cmessage(STATUS "OpenMP not found. Threading not available.")
endif()

add_compile_options(-Wall -Wno-unused-variable -Wno-sign-compare -Wno-unused-function -Wno-unused-but-set-variable -Wno-reorder)