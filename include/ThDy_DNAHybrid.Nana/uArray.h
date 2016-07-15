/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
** @file  ThDySec\include\ThDy_DNAHybrid.Nana\uArray.h
*
* @brief 
*/

#ifndef uArray_H
#define uArray_H

#include <../../nana.ext/include/EditableForm.hpp>
#include "../../nana.ext/include/nanaBind.hpp"
//
//#include <iostream>    // temp, for debugging
//#include <fstream>     // temp, for debugging
//#include <filesystem>
//
//#include <../../nana.ext/include/Numer.hpp>
//
//#include "thdy_programs\init_thdy_prog_param.h"
//#include "matrix.h" 
//#include "common_basics.h" 
//
//
//using namespace ParamGUIBind;
//
class ThDyNanaForm ;
extern std::string e_mail_firma;
// 
//using List = nana::listbox;

class uArray : public CompoWidget
{ public: 
    ThDyNanaForm      &_Pr;
    nana::button  _do_uArray{*this, (" uArray ! ")};
    ParamGUIBind::BindGroup          _uArray;

    uArray (ThDyNanaForm& tdForm);

    void SetDefLayout   () override
    {
        _DefLayout= "vertical      gap=2                                                \n\t"
	        "  < weight=23>                                                             \n\t "
            "  <<><_do_uArray  vertical min=50 max=200><> weight=50>                    \n\t "
			"             <>                                                            \n\t "
			"             < weight = 21 < > <Firma weight = 180> <weight = 3 > >        \n\t "
            ;
    }

    void AsignWidgetToFields() override
    {
	    _place.field("_do_uArray"         )<<_do_uArray;
		_place.field("Firma") << e_mail_firma;
    }

  private: void buttuArray_Click(); //	  Run      _IPrgPar_mPCR
};


class MplexPCR : public CompoWidget
{ public: 
    ThDyNanaForm      &_Pr;
    nana::button  _do_mPCR{*this, (" PCR ! ")};
    ParamGUIBind::BindGroup          _mPCR;

    MplexPCR (ThDyNanaForm& tdForm);

    void SetDefLayout   () override
    {
        _DefLayout= "vertical      gap=2                                                \n\t"
	        "             <  weight=23>                                   \n\t "
            "             < <><_do_mPCR  vertical min=50 max=200> <> weight=50 >        \n\t "
			"             <>                                                            \n\t "
			"             < weight = 21 < > <Firma weight = 180> <weight = 3 > >        \n\t "
            ;
    }
    void AsignWidgetToFields() override
    {
	    _place.field("_do_mPCR" ) << _do_mPCR;
		_place.field("Firma"    ) << e_mail_firma;
	}

  private: void buttPCR_Click(); //	  Run      _IPrgPar_mPCR
};

#endif