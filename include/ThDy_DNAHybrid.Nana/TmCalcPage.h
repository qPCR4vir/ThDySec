/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
*
* @file  ThDySec\include\ThDy_DNAHybrid.Nana\TmCalcPage.h
*
* @brief 
*
*/

#ifndef TmCalcPage_H
#define TmCalcPage_H

#include <nana/gui/widgets/group.hpp>
#include <nana/gui/widgets/checkbox.hpp>
//#include <nana/gui/widgets/progress.hpp>

#include <nanaBind.hpp>
#include <EditableForm.hpp>
#include <number.hpp>

#include "ThDy_programs/init_ThDy_prog_param.h"

class ThDyNanaForm ;
extern std::string e_mail_firma;

class TmCalcPage : public CompoWidget
{
    ThDyNanaForm            &_Pr;
    nana::group             primers             {*this, ("<bold=true> Primers: </>" ), true}, 
                            interaction         {*this, ("<bold=true> Interaction: </>" ), true},
                            align               {*this, ("<bold=true> Alignment: </>" ), true}; 
    nana::textbox           sec_                {primers},  
                            sec2align_          {primers},  
                            txtBx_ResultSec     {align},  
                            txtBx_ResultSec2Align{align};
    nana::checkbox          chkBx_Tm_save_asPCR {*this, ("save")},   
                            chkBx_align         {*this, ("align")},
                            chkBx_copy_rev      {primers, ("rev")},    
                            chkBx_copy_compl    {primers, ("cpl")};
    nana::button            run_                {*this, ("Tm !")},
                            copy_f_s_2          {primers, ("copy")},   
                            copy_s              {primers, ("c")},
                            copy_s_a            {primers, ("c")};      
    nana::label             error_              {primers, ("no error")},
                            _firma;
    ;
    nana::NumberBox         Tm_min_Up{interaction}, Tm_Up{interaction}, Tm_max_Up{interaction} ,
                            Tm_min_Dw{interaction}, Tm_Dw{interaction}, Tm_max_Dw{interaction} ,
                            Tm_min_In{interaction}, Tm_In{interaction}, Tm_max_In{interaction} ,
                            G_min_Up {interaction},  G_Up{interaction},  G_max_Up{interaction} ,
                            G_min_Dw {interaction},  G_Dw{interaction},  G_max_Dw{interaction} ,
                            G_min_In {interaction},  G_In{interaction},  G_max_In{interaction} ;

    ParamGUIBind::BindGroup              _TmCalc;
public:     
    TmCalcPage (ThDyNanaForm& tdForm);

    void SetDefLayout   () override
    {
        _DefLayout= 	
	"vertical      gap=8  min=150    margin=5                              		\n\t"
	"		       < weight=95  primers  >                                               		\n\t"
	"		       < weight=95  gap=2  <weight=120 vertical gap=2  margin=[15,45,0,15]   Left  >              	\n\t"
	"                                                                       <weight=320   Table     >  >           		\n\t"
	"		        < weight=70 vertical  ResAlign>    		\n\t"
		           "<>"
		           "< weight = 21 < > <Firma weight = 180> <weight = 3 > >       \n\t  "
	"		\n\t"
            ;

      primers.div("vert <weight=50  margin=[0,5,0,5] <min=100 vertical gap=2 InputSec>                            "  
                        "                            <weight=50 gap=1 CopyBut grid=[2,2]  collapse(0,0,1,2)> > \n\t "
                        "<weight=23   <weight=20>"
                        "                       <min=50    error     > "
                        "                       <weight=80 rev_compl >     >         \n\t  "
		  );


         interaction.div("vert <min=280    margin=[0,5,5,5] Table    grid=[7,4]    >                "  );

         align.div("vert < weight=50 vertical margin=[0,5,0,5]  ResAlign  >                            " );
    }

    void AsignWidgetToFields() override
    {
	    _place.field("primers"  )<< primers    ;
	    _place.field("Left"     )<< run_        << chkBx_align;
	    _place.field("Table"    )<< interaction;
	    _place.field("ResAlign" )<< align      ;

	    primers["InputSec" ]<< sec_          << sec2align_ ;
	    primers["CopyBut"  ]<<  copy_f_s_2   << copy_s      << copy_s_a ;
	    primers["error"    ]<< error_        ;
	    primers["rev_compl"]<< chkBx_copy_rev << chkBx_copy_compl ;
		_place.field("Firma") << _firma;


	    interaction["Table" ]<< ""          << "   min-" << u8"Tm(�C)"   << "-max"  << "   min-"  << "G(kJ)"    << "-max   "
	                         << "Up"        << Tm_min_Up << Tm_Up        << Tm_max_Up<<G_min_Up   <<  G_Up      <<  G_max_Up   
	                         << "Down"      << Tm_min_Dw << Tm_Dw        << Tm_max_Dw<<G_min_Dw   <<  G_Dw      <<  G_max_Dw   
	                         << "Interact"  << Tm_min_In << Tm_In        << Tm_max_In<<G_min_In   <<  G_In      <<  G_max_In  ;

        align["ResAlign" ]  << txtBx_ResultSec << txtBx_ResultSec2Align ;
    }

    void Run();
    void Copy();
    void Self();
    void Rev();
};

#endif