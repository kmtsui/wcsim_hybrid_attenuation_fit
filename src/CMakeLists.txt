set(HEADERS 
    Likelihoods.hh
    BinManager.hh
    AnaEvent.hh
    AnaTree.hh
    AnaSample.hh
    Fitter.hh
    AnaFitParameters.hh
    ParameterFunction.hh
)

add_library(OpticalFit SHARED 
    BinManager.cc
    AnaTree.cc
    AnaSample.cc
    Fitter.cc
    AnaFitParameters.cc
)

set_target_properties(OpticalFit PROPERTIES 
	PUBLIC_HEADER "${HEADERS}")

LinkToROOT(OpticalFit)

install(TARGETS OpticalFit DESTINATION lib/
		PUBLIC_HEADER DESTINATION include/OpticalFit)

install(DIRECTORY toml DESTINATION include/)
execute_process(COMMAND ln -s ${CMAKE_SOURCE_DIR}/src/toml ${CMAKE_BINARY_DIR}/toml OUTPUT_QUIET ERROR_QUIET)