# files specific to LISEM

# Platform-specific configurations
IF(WIN32)
    # QWT configuration for double axis display, note a double axis branch of qwt is used
    SET(QWT_BUILD_DIR "C:/prgc/lisemgit/qwt/git")    # Adjust to your folder names
    SET(MINGW_BUILD_DIR "c:/qt/msys64/mingw64")     # Adjust to your folder names
    SET(OSSL_LIBRARIES "${MINGW_BUILD_DIR}/lib/libgnutls-openssl.dll.a")
    SET(OSSL_INCLUDE_DIRS "${MINGW_BUILD_DIR}/include/openssl")
    SET(GDAL_INCLUDE_DIRS "${MINGW_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${MINGW_BUILD_DIR}/lib/libgdal.dll.a")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/src")
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")

    FIND_PATH(OMP_INCLUDE_DIRS
        NAMES omp.h
        PATHS "${MINGW_BUILD_DIR}/include"
    )
ENDIF()


#linux
IF(UNIX AND NOT CYGWIN)
    #SET(QWT_BUILD_DIR "/usr/local/qwt-6.4.0-ma")
    IF(DEFINED ENV{NIX_QWT})
         SET(QWT_BUILD_DIR "$ENV{NIX_QWT}")
    ENDIF()
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
    #SET(GDAL_INCLUDE_DIRS "/usr/include/gdal")
    #SET(GDAL_LIBRARIES "/usr/lib/x86_64-linux-gnu/libgdal.so")
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/")
ENDIF()

# Include directories
INCLUDE_DIRECTORIES(
    ${GDAL_INCLUDE_DIRS}
    #${OSSL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${OMP_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
    SYSTEM
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

# Find OpenMP
find_package(OpenMP REQUIRED)

#Find GDAL
find_package(GDAL REQUIRED)

find_package(CURL REQUIRED)

# Enable automatic handling of MOC, UIC, and RCC based on file type changes instead of timestamps
set(CMAKE_AUTOMOC_DEPEND_FILTERS "moc" "*.h")
set(CMAKE_AUTOUIC_DEPEND_FILTERS "ui" "*.ui")
set(CMAKE_AUTORCC_DEPEND_FILTERS "qrc" "*.qrc")

# Optionally skip rule dependency checks to avoid timestamp issues
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set_property(DIRECTORY PROPERTY CMAKE_CONFIGURE_DEPENDS "")

# Compiler flags
IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -Wno-var-tracking-assignments -std=c++11 ${OpenMP_CXX_FLAGS}")
    IF(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wl,-rpath=${ORIGIN}./lib")
    ENDIF()
ENDIF()

# Source files
SET(APP_SOURCES
    fixesandbugs.txt
    maps/CsfMap.cpp
    maps/CsfRGBMap.cpp
    maps/error.cpp
    maps/fixture.cpp
    maps/io.cpp
    maps/operation.cpp
    ui_full/LisUIPatch.cpp
    ui_full/LisUIdialogs.cpp
    ui_full/LisUIScreenshot.cpp
    ui_full/LisUItreecheck.cpp
    ui_full/LisUIModel.cpp
    ui_full/LisUIrunfile.cpp
    ui_full/LisUImapnames.cpp
    ui_full/LisUItreeitem.cpp
    ui_full/LisUItreemodel.cpp
    ui_full/LisUIDefaultNames.cpp
    ui_full/lisemqt.cpp
    ui_full/LisUIplot.cpp
    ui_full/LisUImapplot.cpp
    ui_full/LisUImapplot.h
    ui_full/lismpeg.cpp
    ui_full/lisUIStyle.cpp
    ui_full/lismpeg.h
    ui_full/lisemqt.h
    swatre/swatstep.cpp
    swatre/swatinit.cpp
    swatre/lookup.cpp
    swatre/swatinp.cpp
    channel/lisChannelflood.cpp
    channel/lisChannelflow.cpp
    channel/lisExtendedChannel.cpp
    channel/lisSWOF2DChannel.cpp
    channel/lisDischargein.cpp
    flow/lisFlowBarriers.cpp
    flow/lisGWflow.cpp
    flow/lisKinematic.cpp
    flow/lisOverlandflow.cpp
    flow/lisBoundary.cpp
    hydro/lisInfiltration.cpp
    hydro/lisInterception.cpp
    hydro/lisPercolation.cpp
    hydro/lisSurfstor.cpp
    hydro/lisSoilmoisture.cpp
    meteo/lisRainfall.cpp
    meteo/lisSnowmelt.cpp
    meteo/lisEvaporation.cpp
    model/main.cpp
    model/lisReportfile.cpp
    model/lisReportmaps.cpp
    model/lisRunfile.cpp
    model/lisTotalsMB.cpp
    model/lisModel.cpp
    model/lisDataInit.cpp
    model/lisDataFunctions.cpp
    flow/lisSWOF2Daux.cpp
    flow/lisSWOF2Dopen.cpp
    flow/lisSWOF2DopenMUSCL.cpp
    flow/lisTiledrainflow.cpp
    erosion/lisChannelErosion.cpp
    erosion/lisSWOF2DSediment.cpp
    erosion/lisErosion.cpp
    pest/lisPesticide.cpp
    include/array.h
    include/CsfMap.h
    include/CsfRGBMap.h
    include/pcrtypes.h
    include/csf.h
    include/csfattr.h
    include/csftypes.h
    include/csfimpl.h
    include/lerror.h
    include/fixture.h
    include/global.h
    include/io.h
    include/LisUIoutput.h
    include/masked_raster.h
    include/mmath.h
    include/model.h
    include/operation.h
    include/option.h
    include/raster.h
    include/swatre_p.h
    include/TMmapVariables.h
    include/VectormapVariables.h
    include/version.h
    PCR/create2.c
    PCR/mclose.c
    PCR/ruseas.c
    PCR/gvalscal.c
    PCR/gcellrep.c
    PCR/putsomec.c
    PCR/setangle.c
    PCR/kernlcsf.c
    PCR/gproj.c
    PCR/csfglob.c
    PCR/setvtmv.c
    PCR/dumconv.c
    PCR/swapio.c
    openlisemico.rc
)

qt6_wrap_cpp(MOC_FILES
    include/model.h
    ui_full/lisemqt.h
    ui_full/lismpeg.h
    # Add all header files with Q_OBJECT here
)

# Generate UI source files
qt_wrap_ui(UI_SOURCES ui_full/lisemqt.ui ui_full/lismpeg.ui)

# Generate resource source files
qt_add_resources(RCC_SOURCES resources/openlisem.qrc)

# Add executable target
add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
    ${MOC_FILES}
)

# Link the necessary libraries
target_link_libraries(Lisem
    Qt6::Widgets Qt6::Gui Qt6::Core Qt6::Network
    ${GDAL_LIBRARIES} ${QWT_LIBRARIES} ${OSSL_LIBRARIES}
    OpenMP::OpenMP_CXX
)

