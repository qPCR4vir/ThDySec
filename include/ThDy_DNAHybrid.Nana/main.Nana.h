/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDy_DNAHybrid.Nana\main.Nana.h
*
* @brief 
*/

#ifndef main_nana_H
#define main_nana_H

#include "SetupPage.h"
#include "SeqExpl.h"
#include "FindSondenPage.h"
//#include "MplexPCR.h"
#include "uArray.h"
#include "TmCalcPage.h"

#include "thdy_programs\init_thdy_prog_param.h"

#include <nana/gui/widgets/tabbar.hpp>
//#include <nana/gui/tooltip.hpp>
//#include <nana/gui/widgets/toolbar.hpp>
//#include <nana/gui/widgets/progress.hpp>
//#include <nana/gui/widgets/group.hpp>
extern std::string e_mail_firma;

#include <nana/gui/wvl.hpp>

class ThDyNanaForm : public nana::form, public EditableForm , public ThDyProject
{
    using tabbar = nana::tabbar<std::string> ;
	tabbar                     tabbar_     {*this};
    SetupPage                  setup_      {*this};
    FindSondenPage             findSond_   {*this};
    MplexPCR                   mPCR_       {*this};
    uArray                     uArr_       {*this}; 
    TmCalcPage                 tmCalc_     {*this}; 
    nana::label                _firma     {*this, e_mail_firma};

  public:    
    std::vector<std::unique_ptr<nana::form>> _results;
    SeqExpl                         mExpl_      {*this};

    ThDyNanaForm (int argc, char *argv[])  ;
    //~ThDyNanaForm();

    void SetDefLayout       () override;
    void AsignWidgetToFields() override;
    void add_page           (widget& w);
    void ShowExpl           (){tabbar_.activated(1);}
};


#endif