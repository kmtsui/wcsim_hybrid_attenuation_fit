IF(USE_WCSIM)
	foreach(app analysis_absorption)
		add_executable(${app} ${app}.cc)
		target_link_libraries(${app} OpticalFit)
		LinkToWCSIM(${app})
		LinkToROOT(${app})

		install(TARGETS ${app} DESTINATION bin)
	endforeach(app)
ENDIF()

foreach(app optical_fit)
	add_executable(${app} ${app}.cc)
	target_link_libraries(${app} OpticalFit)
	LinkToROOT(${app})

	install(TARGETS ${app} DESTINATION bin)
endforeach(app)


install(DIRECTORY config DESTINATION var/OPTICALFIT)