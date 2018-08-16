cmake_minimum_required(VERSION 3.12)
project(ThDy_DNAHybrid.gcc8 VERSION      0.02.00
        DESCRIPTION "Thermo Dynamic DNA Hybridizations"
        HOMEPAGE_URL https://github.com/qPCR4vir/ThDySec
        LANGUAGES CXX         )




# TODO: Enforce compatibility with the already compiled nana
option(NANA_CMAKE_FIND_BOOST_FILESYSTEM "Search: Is Boost filesystem available?" OFF)
#option(NANA_CMAKE_NANA_FILESYSTEM_FORCE "Force nana filesystem over ISO and boost?" ON)
option(NANA_CMAKE_STD_FILESYSTEM_FORCE "Use of STD filesystem?(a compilation error will ocurre if not available)" OFF)
option(NANA_CMAKE_BOOST_FILESYSTEM_FORCE "Force use of Boost filesystem if available (over STD)?" OFF)

set(NANA_LIB_HINT "C:/Prog/ExtLib" CACHE PATH "Hint for the directory where the directories lib and include are")

option(NANA_FORCE_INSTALL_PATH "Use my onw absolute path and not DESTDIR/CMAKE_INSTALL_PREFIX (CLion workaround!)" OFF)

set(CMAKE_DEBUG_POSTFIX "_d")

########### OS

if(WIN32)
    add_definitions(-DWIN32)
    set(NANA_INSTALL_PATH ${NANA_INSTALL_PATH_W})
    #Global MSVC definitions. You may prefer the hand-tuned sln and projects from the nana repository.
    if(MSVC)
        option(WIN32_USE_MP "Set to ON to build nana with the /MP option (Visual Studio 2005 and above)." ON)
        # ??
        if(WIN32_USE_MP)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
        endif(WIN32_USE_MP)
    endif(MSVC)
endif(WIN32)

if(APPLE)
    add_definitions(-DAPPLE)
    include_directories(/opt/X11/include/)
    set(NANA_LINKS "${NANA_LINKS}-L/opt/X11/lib/ -liconv")
    set(ENABLE_AUDIO OFF)
elseif(UNIX)
    add_definitions(-Dlinux)
    message("added -D linux")
endif(APPLE)

if(UNIX)
    set(NANA_LINKS "${NANA_LINKS}-lX11 -lfontconfig")
    set(NANA_INSTALL_PATH ${NANA_INSTALL_PATH_L})
    find_package(Freetype)
    if (FREETYPE_FOUND)
        include_directories( ${FREETYPE_INCLUDE_DIRS})
        set(NANA_LINKS "${NANA_LINKS}-lXft")
    endif(FREETYPE_FOUND)
endif(UNIX)


########### Compilers
#
#
# see at end of:  https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_dynamic_or_shared.html

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")    # Clang || GNU

    if  ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -Wall -g")       # Clang

    else ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -Wall -g")       # GNU

    endif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

endif(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

# enable static linkage     # GNU || CLang not MinGW
if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") #  AND NOT MINGW
    # set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++  -pthread")
    set(NANA_LINKS "${NANA_LINKS}-static-libgcc -static-libstdc++  -pthread")
    # message("Setting NANA_LINKS to -static-libgcc -static-libstdc++  -pthread or ${NANA_LINKS}")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.3)
        # IS_GNUCXX < 5.3
    else(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.3)
        # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lstdc++fs")    # IS_GNUCXX 5.3 or more
        set(NANA_LINKS "${NANA_LINKS}-lstdc++fs")
    endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.3)

endif(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") #  AND NOT MINGW


if (APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")                    # APPLE Clang
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")
    set(NANA_LINKS "${NANA_LINKS}-stdlib=libstdc++")
endif ()


# Find/Select filesystem
if (1) #NANA_CMAKE_NANA_FILESYSTEM_FORCE)
    add_definitions(-DNANA_FILESYSTEM_FORCE)
elseif (NANA_CMAKE_STD_FILESYSTEM_FORCE)
    add_definitions(-DSTD_FILESYSTEM_FORCE)
elseif (NANA_CMAKE_FIND_BOOST_FILESYSTEM OR NANA_CMAKE_BOOST_FILESYSTEM_FORCE)
    if (NANA_CMAKE_BOOST_FILESYSTEM_FORCE)
        add_definitions(-DBOOST_FILESYSTEM_FORCE)
    endif(NANA_CMAKE_BOOST_FILESYSTEM_FORCE)
    # Some CMake versions don't list missing boost dependencies, so we need to add them manually
    # - libboost_system is needed for HelloWorld demo at least
    # - libboost_thread is needed for threading demo
    # - libboost_chrono is needed for thread_pool demo
    # Boost_THREADAPI should be set to correct suffix (empty, or win32, or pthread)
    find_package(Boost COMPONENTS filesystem system thread chrono)
    if (Boost_FOUND)
        add_definitions(-DBOOST_FILESYSTEM_AVAILABLE)
        include_directories(SYSTEM "${Boost_INCLUDE_DIR}")

        #set(NANA_LINKS "${NANA_LINKS} ${Boost_LIBRARIES}")    ######  FIRST !!!!!!!!!!!!!!!!!    add   is not first
        set(NANA_LINKS "${NANA_LINKS}-lboost_filesystem -lboost_system")
    endif (Boost_FOUND)
    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_STATIC_RUNTIME ON)  # ??
endif (NANA_CMAKE_NANA_FILESYSTEM_FORCE)


######## Nana options

add_definitions(-DNANA_IGNORE_CONF)
if(NANA_CMAKE_VERBOSE_PREPROCESSOR)
    add_definitions(-DVERBOSE_PREPROCESSOR)
endif(NANA_CMAKE_VERBOSE_PREPROCESSOR)
if(NANA_CMAKE_AUTOMATIC_GUI_TESTING)
    add_definitions(-DNANA_AUTOMATIC_GUI_TESTING)
    enable_testing ()
endif(NANA_CMAKE_AUTOMATIC_GUI_TESTING)

if (NANA_FORCE_INSTALL_PATH OR NANA_CLION)
    set(NANA_INSTALL_PATH "${NANA_INSTALL_PATH}/")
else()
    set(NANA_INSTALL_PATH)
endif()


#######################     Main setting of Nana lib, demos, sources, targets and install

set(nana_lib_file_name nana)    #${CMAKE_DEBUG_POSTFIX}

find_library(nana_lib ${nana_lib_file_name}
        C:/Prog/ExtLib/nanalib/lib
        "${NANA_LIB_HINT}/lib"
        "${NANA_LIB_HINT}"
        "${NANA_LIB_HINT}/nanalib"
        "${NANA_LIB_HINT}/nanalib/lib"
        "${NANA_LIB_HINT}/nana_lib"
        "${NANA_LIB_HINT}/nana_lib/lib"
        "${NANA_LIB_HINT}/build"
        "${NANA_LIB_HINT}/build/bin"
        "${CMAKE_INSTALL_PREFIX}/../nana_lib/lib"
        "${CMAKE_INSTALL_PREFIX}/../nanalib/lib"
        "${CMAKE_INSTALL_PREFIX}/../nana/lib"
        "${CMAKE_INSTALL_PREFIX}/../lib"
        "${CMAKE_INSTALL_PREFIX}/../nana/build"
        "${CMAKE_INSTALL_PREFIX}/../nana/build/bin"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana_lib/lib"
        "${CMAKE_CURRENT_BINARY_DIR}/../nanalib/lib"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana/lib"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana/build"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana/build/bin"
        "${CMAKE_CURRENT_BINARY_DIR}/../lib"
        "${CMAKE_CURRENT_DIR}/../nana_lib/lib"
        "${CMAKE_CURRENT_DIR}/../nanalib/lib"
        "${CMAKE_CURRENT_DIR}/../nana/lib"
        "${CMAKE_CURRENT_DIR}/../lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana_lib/lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nanalib/lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana/lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/../lib"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana/build"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana/build/bin"
        )
find_path(nana_include nana/config.hpp
        C:/Prog/ExtLib/nanalib/include
        "${NANA_LIB_HINT}/include"
        "${NANA_LIB_HINT}"
        "${NANA_LIB_HINT}/nanalib"
        "${NANA_LIB_HINT}/nanalib/include"
        "${NANA_LIB_HINT}/nana_lib"
        "${NANA_LIB_HINT}/nana_lib/include"
        "${CMAKE_INSTALL_PREFIX}/../nana_lib/include"
        "${CMAKE_INSTALL_PREFIX}/../nanalib/include"
        "${CMAKE_INSTALL_PREFIX}/../nana/include"
        "${CMAKE_INSTALL_PREFIX}/../include"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana_lib/include"
        "${CMAKE_CURRENT_BINARY_DIR}/../nanalib/include"
        "${CMAKE_CURRENT_BINARY_DIR}/../nana/include"
        "${CMAKE_CURRENT_BINARY_DIR}/../include"
        "${CMAKE_CURRENT_DIR}/../nana_lib/include"
        "${CMAKE_CURRENT_DIR}/../nanalib/include"
        "${CMAKE_CURRENT_DIR}/../nana/include"
        "${CMAKE_CURRENT_DIR}/../include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana_lib/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nanalib/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../nana/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../include"
        DOC  "Nana Include directory")
include_directories(${nana_include})

message("nana_lib_file_name = " ${nana_lib_file_name})
message("nana_lib = " ${nana_lib})
message("nana_include= " ${nana_include})




add_executable(${PROJECT_NAME}
        ../ExtLib/nana.ext/src/CompoWidget.cpp
        ../ExtLib/nana.ext/src/EditableForm.cpp
        ../ExtLib/nana.ext/src/EditLayot_Form_Demo.cpp
        ../ExtLib/nana.ext/src/number.cpp

        ../ProgParam/src/init_prog_param.cpp
        ../ProgParam/src/ParamGUIBind.cpp
        ../ProgParam/src/Units.cpp

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
        src/ThDy_DNAHybrid.Nana/uArray.cpp

        )

set_property( TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14 )
#target_include_directories(${PROJECT_NAME} PUBLIC src)
target_link_libraries(${PROJECT_NAME} ${nana_lib} ${NANA_LINKS} )
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "ThDy")

# Just for information:
message ( "CMAKE_CXX_COMPILER_ID   = "  ${CMAKE_CXX_COMPILER_ID})
message ( "COMPILER_IS_CLANG       = "  ${COMPILER_IS_CLANG})
message ( "CMAKE_COMPILER_IS_GNUCXX= "  ${CMAKE_COMPILER_IS_GNUCXX})
message ( "CMAKE_CXX_COMPILER_VERSION= "  ${CMAKE_CXX_COMPILER_VERSION})
message ( "CMAKE_CXX_FLAGS         = "  ${CMAKE_CXX_FLAGS})
message ( "CMAKE_EXE_LINKER_FLAGS  = "  ${CMAKE_EXE_LINKER_FLAGS})
message ( "CMAKE_STATIC_LINKER_FLAGS= "  ${CMAKE_STATIC_LINKER_FLAGS})
message ( "NANA_LINKS              = "  ${NANA_LINKS})
message ( "DESTDIR                 = "  ${DESTDIR})
message ( "CMAKE_INSTALL_PREFIX    = "  ${CMAKE_INSTALL_PREFIX})
message ( "NANA_CMAKE_ENABLE_AUDIO = "  ${NANA_CMAKE_ENABLE_AUDIO})
message ( "CMAKE_CURRENT_SOURCE_DIR= "  ${CMAKE_CURRENT_SOURCE_DIR})
message ( "CMAKE_CURRENT_DIR       = "  ${CMAKE_CURRENT_DIR})
message ( "CMAKE_CURRENT_BINARY_DIR= "  ${CMAKE_CURRENT_BINARY_DIR})
message ( "CMAKE_DEBUG_POSTFIX     = "  ${CMAKE_DEBUG_POSTFIX})

message ( "NANA_CMAKE_NANA_FILESYSTEM_FORCE         = "  ${NANA_CMAKE_NANA_FILESYSTEM_FORCE})
message ( "NANA_CMAKE_FIND_BOOST_FILESYSTEM         = "  ${NANA_CMAKE_FIND_BOOST_FILESYSTEM})
message ( "NANA_CMAKE_BOOST_FILESYSTEM_FORCE        = "  ${NANA_CMAKE_BOOST_FILESYSTEM_FORCE})
message ( "NANA_CMAKE_BOOST_FILESYSTEM_INCLUDE_ROOT = "  ${NANA_CMAKE_BOOST_FILESYSTEM_INCLUDE_ROOT})
message ( "NANA_CMAKE_BOOST_FILESYSTEM_LIB          = "  ${NANA_CMAKE_BOOST_FILESYSTEM_LIB})
message ( "NANA_CMAKE_AUTOMATIC_GUI_TESTING         = "  ${NANA_CMAKE_AUTOMATIC_GUI_TESTING})
message ( "NANA_CMAKE_ADD_DEF_AUTOMATIC_GUI_TESTING = "  ${NANA_CMAKE_ADD_DEF_AUTOMATIC_GUI_TESTING})