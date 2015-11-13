add_definitions (-D_GNU_SOURCE)

find_package( OpenCV REQUIRED )
find_package( CFITSIO REQUIRED )

# add the executable
SET(SOLAR_ANALYZER_FILES
    ${PROJECT_SOURCE_DIR}/src/SolarAnalyzer.cpp
    ${PROJECT_SOURCE_DIR}/src/SolarAnalyzerBlob.cpp
    ${PROJECT_SOURCE_DIR}/src/SolarAnalyzerDetect.cpp
    ${PROJECT_SOURCE_DIR}/src/SolarAnalyzerGroup.cpp
    )

add_executable (sunspot detect.cpp blob.cpp main.cpp mag.cpp group.cpp fitstool.cpp ${SOLAR_ANALYZER_FILES})

target_link_libraries( sunspot ${OpenCV_LIBS} ${CFITSIO_LIBRARY})
