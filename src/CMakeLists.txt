set(HEADERS 
    Likelihoods.hh
    BinManager.hh
    AnaEvent.hh
    AnaTree.hh
    AnaSample.hh
    Fitter.hh
    AnaFitParameters.hh
    ParameterFunction.hh
    EigenDecomp.hh
    ToyThrower.hh
    ColorOutput.hh
)

add_library(OpticalFit SHARED 
    BinManager.cc
    AnaTree.cc
    AnaSample.cc
    Fitter.cc
    AnaFitParameters.cc
    EigenDecomp.cc
    ToyThrower.cc
)

set_target_properties(OpticalFit PROPERTIES 
	PUBLIC_HEADER "${HEADERS}")

LinkToROOT(OpticalFit)

if( USE_CUDA )
    target_link_libraries( OpticalFit CudaCache )
    target_include_directories( OpticalFit PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/CacheManager/include
        ${CMAKE_SOURCE_DIR}/src)
endif()

install(TARGETS OpticalFit DESTINATION lib/
		PUBLIC_HEADER DESTINATION include/OpticalFit)

install(DIRECTORY toml DESTINATION include/)
execute_process(COMMAND ln -s ${CMAKE_SOURCE_DIR}/src/toml ${CMAKE_BINARY_DIR}/toml OUTPUT_QUIET ERROR_QUIET)

install(DIRECTORY utils DESTINATION include/)
execute_process(COMMAND ln -s ${CMAKE_SOURCE_DIR}/src/utils ${CMAKE_BINARY_DIR}/utils OUTPUT_QUIET ERROR_QUIET)