/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
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

#include <nana/gui.hpp>
#include <nana/gui/widgets/tabbar.hpp>

#include "ThDy_programs/init_ThDy_prog_param.h"

#include "SetupPage.h"
#include "SeqExpl.h"
#include "FindSondenPage.h"
#include "uArray.h"
#include "TmCalcPage.h"

class ThDyNanaForm : public nana::form, public EditableForm , public ThDyProject
{
  public:
    const std::string          e_mail_firma;

  private:
    using tabbar = nana::tabbar<std::string> ;
    tabbar                     tabbar_     {*this};
    SetupPage                  setup_      {*this};
    FindSondenPage             findSond_   {*this};
    MplexPCR                   mPCR_       {*this};
    uArray                     uArr_       {*this}; 
    TmCalcPage                 tmCalc_     {*this}; 
    nana::label                _firma     {*this, this->e_mail_firma};

  public:    
    std::vector<std::unique_ptr<nana::form>> _results;
    SeqExpl                         mExpl_      {*this};

    ThDyNanaForm (int argc, char *argv[])  ;

    void SetDefLayout       () override;
    void AsignWidgetToFields() override;
    void add_page           (widget& w);
    void ShowExpl           (){tabbar_.activated(1);}
};


#endif