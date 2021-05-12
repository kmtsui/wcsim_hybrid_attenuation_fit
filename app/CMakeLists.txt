foreach(app analysis_absorption)
	add_executable(${app} ${app}.cc)
#	target_link_libraries(${app} NIWGWeightEngines)
    LinkToWCSIM(${app})
	LinkToROOT(${app})

	install(TARGETS ${app} DESTINATION bin)
endforeach(app)


#install(DIRECTORY config DESTINATION var/OPTICALFIT)