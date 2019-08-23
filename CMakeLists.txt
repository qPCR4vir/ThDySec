cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(ThDy_DNAHybrid VERSION      0.02.02
        DESCRIPTION "Thermo Dynamic DNA Hybridizations"
        HOMEPAGE_URL https://github.com/qPCR4vir/ThDySec
        LANGUAGES CXX         )

#######################     Main setting of sources, targets and install

add_executable(   ThDy_DNAHybrid)

target_include_directories(ThDy_DNAHybrid PUBLIC include)
#target_compile_features   (ThDy_DNAHybrid PUBLIC cxx_std_17)
target_sources (ThDy_DNAHybrid PRIVATE
        src/ThDySec/cod_deg.cpp
        src/ThDySec/sec.cpp
        src/ThDySec/sec_basic.cpp
        src/ThDySec/sec_mult.cpp
        src/ThDySec/sec_rang.cpp
        src/ThDySec/th_dy_align.cpp
        src/ThDySec/th_dy_param.cpp

        src/ThDy_programs/init_thdy_prog_param.cpp
        src/ThDy_programs/prog_comm_functions.cpp
        src/ThDy_programs/ThDy_DegTmCalc.cpp
        src/ThDy_programs/ThDy_microArray.cpp
        src/ThDy_programs/ThDy_MultiplexPCR.cpp
        src/ThDy_programs/ThDy_SondeDesign.cpp

        src/ThDy_DNAHybrid.Nana/FindSondenPage.cpp
        src/ThDy_DNAHybrid.Nana/main.Nana.cpp
        src/ThDy_DNAHybrid.Nana/SeqExpl.cpp
        src/ThDy_DNAHybrid.Nana/SetupPage.cpp
        src/ThDy_DNAHybrid.Nana/TableResults.cpp
        src/ThDy_DNAHybrid.Nana/TmCalcPage.cpp
        src/ThDy_DNAHybrid.Nana/uArray.cpp)


add_subdirectory(../ExtLib/nana.ext  ../ExtLib/cmake-nana_ext-build-${CMAKE_BUILD_TYPE}
                 EXCLUDE_FROM_ALL )
target_link_libraries(ThDy_DNAHybrid PRIVATE nana_ext )


install(TARGETS ThDy_DNAHybrid RUNTIME DESTINATION "ThDy")

if (NANA_CMAKE_VERBOSE_PREPROCESSOR)

    include(CMakePrintHelpers)
    cmake_print_properties(TARGETS ThDy_DNAHybrid  PROPERTIES
            COMPILE_DEFINITIONS   COMPILE_OPTIONS    COMPILE_FLAGS  LINK_LIBRARIES
            INCLUDE_DIRECTORIES   INSTALL_NAME_DIR   LINK_FLAGS     VERSION)

endif()