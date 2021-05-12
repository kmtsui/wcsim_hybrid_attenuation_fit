set(HEADERS 
    Likelihoods.hh
    BinManager.hh
    AnaEvent.hh
    AnaTree.hh
    AnaSample.hh
)

add_library(OpticalFit SHARED 
    BinManager.cc
    AnaTree.cc
    AnaSample.cc
)

set_target_properties(OpticalFit PROPERTIES 
	PUBLIC_HEADER "${HEADERS}")

LinkToROOT(OpticalFit)

install(TARGETS OpticalFit DESTINATION lib/
		PUBLIC_HEADER DESTINATION include/OpticalFit)