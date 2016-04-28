/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\TmCalcPage.cpp
*
* @brief 
*
*/

#include "ThDy_DNAHybrid.Nana\TmCalcPage.h"
#include "ThDy_DNAHybrid.Nana\main.Nana.h"

//#include <iostream>    // temp, for debugging
//#include <fstream>     // temp, for debugging
//#include <filesystem>
//
//#include "thdy_programs\init_thdy_prog_param.h"
//#include "matrix.h" 
//#include "common_basics.h" 

//#include <../../nana.ext/include/EditableForm.hpp>
//#include <../../nana.ext/include/Numer.hpp>
//#include "../../nana.ext/include/nanaBind.hpp"


TmCalcPage::TmCalcPage        (ThDyNanaForm& tdForm) try
        : _Pr           (tdForm), 
          CompoWidget  (tdForm, ("Tm Calc"), ("Tm Calc.lay.txt"))
    {
                         sec_.multi_lines(false).editable(true ).tip_string (("forward primer"));
                   sec2align_.multi_lines(false).editable(true ).tip_string (("reverse primer"));
              txtBx_ResultSec.multi_lines(false).editable(false).tip_string (("alingned forward primer"));
        txtBx_ResultSec2Align.multi_lines(false).editable(false).tip_string (("alingned reverse primer"));
        
        using ParamGUIBind::link;

        _TmCalc << link (   _Pr._TmCal.align      ,    chkBx_align    )    
                << link (   _Pr._TmCal._Sec       ,    sec_           )
                << link (   _Pr._TmCal._Sec2Align ,    sec2align_     )
                ;

        run_      .events().click([&](){Run ();});
        copy_f_s_2.events().click([&](){Copy();});      ;   //(*this, ("copy")),   
        copy_s    .events().click([&](){Self();});      ;  //  (*this, ("c")),
        copy_s_a  .events().click([&](){Rev ();});      ;  

        InitMyLayout();
        SelectClickableWidget( *this);
        SelectClickableWidget( error_);
    }
catch (std::exception & e)
{
	throw std::runtime_error(std::string("An error ocurred during initialization of the Tm Calc page window:\n") + e.what());
}
catch (...)
{
	throw std::runtime_error(std::string("An unknonw error ocurred during initialization of the Tm Calc page window"));
}