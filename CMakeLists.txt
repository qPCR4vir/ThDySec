cmake_minimum_required(VERSION 3.12)
project(ThDy_DNAHybrid VERSION      0.02.00
        DESCRIPTION "Thermo Dynamic DNA Hybridizations"
        HOMEPAGE_URL https://github.com/qPCR4vir/ThDySec
        LANGUAGES CXX         )

#######################     Main setting of sources, targets and install

include(../ExtLib/nana/CMakeLists.txt)


add_executable(   ThDy_DNAHybrid
        ../ProgParam/src/init_prog_param.cpp
        ../ProgParam/src/ParamGUIBind.cpp
        ../ProgParam/src/Units.cpp

        ../ExtLib/nana.ext/src/CompoWidget.cpp
        ../ExtLib/nana.ext/src/EditableForm.cpp
        ../ExtLib/nana.ext/src/number.cpp

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

include_directories(include)
include_directories(../ProgParam/include)


set_property( TARGET ThDy_DNAHybrid PROPERTY CXX_STANDARD 14 )
#target_include_directories(${PROJECT_NAME} PUBLIC src)
target_link_libraries(ThDy_DNAHybrid nana )
install(TARGETS ThDy_DNAHybrid RUNTIME DESTINATION "ThDy")

message ( "NANA_CMAKE_ADD_DEF_AUTOMATIC_GUI_TESTING = "  ${NANA_CMAKE_ADD_DEF_AUTOMATIC_GUI_TESTING})