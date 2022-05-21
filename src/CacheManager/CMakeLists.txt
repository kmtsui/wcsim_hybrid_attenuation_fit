if( USE_CUDA )
  set(SRCFILES
          src/CacheManager.cpp
          )
  set(HEADERS
          include/CacheManager.h
          include/CacheParameters.h
          include/CacheWeights.h
          include/WeightNormalization.h
          include/WeightAttenuation.h
          include/WeightAttenuationZ.h
          include/WeightPolynomialCosth.h
          include/WeightSourcePhiVar.h
          include/WeightBase.h
          include/CacheIndexedSums.h
          )

  if(CMAKE_CUDA_COMPILER)
    cmessage(STATUS "CUDA utilities being compiled")
    set(SRCFILES ${SRCFILES} src/WeightNormalization.cu)
    set(SRCFILES ${SRCFILES} src/WeightAttenuation.cu)
    set(SRCFILES ${SRCFILES} src/WeightAttenuationZ.cu)
    set(SRCFILES ${SRCFILES} src/WeightPolynomialCosth.cu)
    set(SRCFILES ${SRCFILES} src/WeightSourcePhiVar.cu)
    set(SRCFILES ${SRCFILES} src/CacheParameters.cu)
    set(SRCFILES ${SRCFILES} src/CacheWeights.cu)
    set(SRCFILES ${SRCFILES} src/CacheIndexedSums.cu)
  else(CMAKE_CUDA_COMPILER)
    cmessage(STATUS "CUDA utilities are not being compiled")
    set(SRCFILES ${SRCFILES} src/WeightNormalization.cpp)
    set(SRCFILES ${SRCFILES} src/WeightAttenuation.cpp)
    set(SRCFILES ${SRCFILES} src/WeightAttenuationZ.cpp)
    set(SRCFILES ${SRCFILES} src/WeightPolynomialCosth.cpp)
    set(SRCFILES ${SRCFILES} src/WeightSourcePhiVar.cpp)
    set(SRCFILES ${SRCFILES} src/CacheParameters.cpp)
    set(SRCFILES ${SRCFILES} src/CacheWeights.cpp)
    set(SRCFILES ${SRCFILES} src/CacheIndexedSums.cpp)
  endif(CMAKE_CUDA_COMPILER)

  add_library(CudaCache SHARED ${SRCFILES})

  set_target_properties(CudaCache PROPERTIES 
	                      PUBLIC_HEADER "${HEADERS}")

  # Make sure the current directories are available for the later
  # compilation.
  target_include_directories( CudaCache PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src)

  LinkToROOT(CudaCache)

  install(TARGETS CudaCache DESTINATION lib
          PUBLIC_HEADER DESTINATION include/CacheManager)

endif()
