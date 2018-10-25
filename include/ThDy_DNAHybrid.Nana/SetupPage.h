/**
* Copyright (C) 2009-2018, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2018
*
* @file  ThDySec\include\ThDy_DNAHybrid.Nana\SetupPage.h
*
* @brief 
*
*/

#ifndef SetupPage_H
#define SetupPage_H

#include <EditableForm.hpp>    // todo: <nana_ext/EditableForm.hpp>
#include <number.hpp>
#include <nanaBind.hpp>
#include <Units.hpp>

#include <nana/gui/widgets/group.hpp>
#include <nana/gui/widgets/checkbox.hpp>

//#include <nana/gui/widgets/progress.hpp>
//#include <nana/gui/widgets/tabbar.hpp>
//#include <nana/gui/widgets/treebox.hpp>
//#include <nana/gui/widgets/listbox.hpp>
//#include <nana/gui/widgets/toolbar.hpp>
//#include <nana/gui/tooltip.hpp>
//#include <nana/gui/wvl.hpp>

//#include <iostream>    // temp, for debugging
//#include <fstream>     // temp, for debugging
//#include <filesystem>


class ThDyNanaForm ;
 

class SetupPage : public CompoWidget
{
    ThDyNanaForm       &_Pr;
    nana::group         _gr_dir			{*this,    (" <bold=true> Directories: </>"               ), true};
    FilePickBox         _results            { _gr_dir, ("Results: ") } ;

    nana::group         _gr_seq         {_gr_dir, (" <bold=true> Sequences </>File or Directory"      ), true};
    nana::group         _gr_targ            {_gr_seq, (""                                                 ), true};
    FilePickBox         _targets            { _gr_targ, ("Targets: ") }  ;
    nana::checkbox      _chkTargRecDir      { _gr_targ, ("Targets - Recur Dir") },
                        _chkTargOnlyStruct  { _gr_targ, ("Only reproduce Dir Structure") };

    nana::group         _gr_ntarg       {_gr_seq, (""                                                 ), true};
    FilePickBox         _nTsec              { _gr_ntarg, ("Non targets: "),("FindSonden-OSB.NonTarg.lay.txt")};
    nana::checkbox      _chk_nTgRecDir      { _gr_ntarg, ("Non Targets - Recur Dir") },
                        _chk_nTgOnlyStruct  { _gr_ntarg, ("Only reproduce Dir Structure") };

    nana::group         _gr_PCRfiltre   {_gr_seq, (""                                             ), true};
    FilePickBox         _PCRfiltre          { _gr_PCRfiltre, ("PCR-filtre: ")};

    nana::group         _gr_PrimersFilePCR{_gr_seq, (""                                                 ), true};
    FilePickBox         _PrimersFilePCR     { _gr_PrimersFilePCR, ("Primers: ") };
    nana::checkbox      _chkPrimRecDir      { _gr_PrimersFilePCR, ("Primers - Recur Dir") },
                        _chkPrOnlyStruct    { _gr_PrimersFilePCR, ("Only reproduce Dir Structure") };

    nana::group         _gr_uArr        {_gr_seq, (""                                                 ), true};
    FilePickBox         _Prob_uArr          { _gr_uArr, ("Probes: ") };
    nana::checkbox      _chkProbRecDir      { _gr_uArr, ("Probes - Recur Dir") },
                        _chkProbOnlyStruct  { _gr_uArr, ("Only reproduce Dir Structure") };

    OpenSaveBox         _NNParamFile        { _gr_dir, ("NN param: ")};

    nana::group         _gr_salt        {*this, (" <bold=true> Input & analisis parameters: </>"                 ), true};
    nana::combox        comBoxSalMeth       { _gr_salt}, 
                        comBoxTAMeth        { _gr_salt};
    nana::NumUnitUpDown numUpDowTgConc      { _gr_salt, ("Target Conctr:"      ), 50, 0.1    , 1000000,  "nM"}, 
                        numUpDowSalConc     { _gr_salt, ("Salt Conc [Cations]:"), 50, 0.0001 , 10000,    "mM"} , 
                        numUpDowTa          { _gr_salt, ("Temp. Anneling:"     ), 55, 40     , 75,       RTunits::grC},
                        numUpDowSdConc      { _gr_salt, ("Sonde Conctr:"       ), 0.8, 0.001 , 1000,     RTunits::uM}  ;

    nana::button        _set_def_proj       { *this,("Set as Def. project") },
                        _load_def_proj      { *this,("ReLoad Def. project") };

    nana::group         _gr_checks {*this, (" <bold=true> Save in results: </>"                 ), true};
    nana::checkbox      ckBx_savTm          { _gr_checks, ("Tm"    ) },
                        ckBx_savPos         { _gr_checks, ("Pos"   ) },
                        ckBx_savG           { _gr_checks, ("G"     ) },
                        ckBx_savAlign       { _gr_checks, ("Align" ) },
                        ckBx_savProj        { _gr_checks, ("Proj"  ) },
                        ckBx_savG_Plasm     { _gr_checks, ("G->Plasmid") },
                        ckBx_savTm_Plasm    { _gr_checks, ("Tm->Plasmid") },
                        ckBx_savLog         { _gr_checks, ("log"     ) },
                        ckBx_savExportSond  { _gr_checks, ("Exp. probes" ) },
                        ckBx_savExportTarg  { _gr_checks, ("Exp. targets") },
                        ckBx_loadNNParam    { _gr_dir   , ("load at start") },
                        ckBx_savNNParam     { _gr_checks, ("save NNparam") }/*,*/
                        ;

    ParamGUIBind::BindGroup  _setup;

    void  SetDefLayout   () override ;
    virtual void  AsignWidgetToFields () final  override;
    void  MakeResponive();
    void  SaveProj();
    void  setAsDefProject();
    void  RestDefPr	 ( )	;	// Restore (USE) Deff  Proj File

public:     
    OpenSaveBox         _proj       { *this, ("Project:") };

    SetupPage (ThDyNanaForm& tdForm);

    static FilePickBox& AddFastaFiltre(FilePickBox &fpb)
    {
        return fpb.add_filter(FastaFiltre( ));
    }
    static FilePickBox::filtres FastaFiltre( );


    void AddMenuItems(nana::menu& menu);
    void LoadProjectAndReg(std::string file);
    void LoadProject(std::string file);
};

#endif