/**
* Copyright (C) 2009-2018, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\SetupPage.cpp
*
* @brief 
*
*/

#include "ThDy_DNAHybrid.Nana\SetupPage.h"
#include "ThDy_DNAHybrid.Nana\main.Nana.h"


using namespace ParamGUIBind;

void   SetupPage::SetDefLayout()
{
	_DefLayout =
R"(
	vertical      gap=3    margin=10    			
	    <min=465  gap=5    <weight=5>   <min=450   vertical   gap=5 	    	 
	                                                   <weight=26  Project       >       	 
	    	                                           <weight=400 dir>       	 
	    	                                           <weight=30   <weight=20>   <min=280 max=700  gap=5 buttons>   <>>	  
	                                      >   <weight=5>   <weight=120 vertical   
                                                                         <weight=235  checks > 
                                                                         <>  
                                                           >    		 
	    
	    >			                                 
	    < weight=70  salt    >   
	    <>				
   		< weight=21 <><Firma weight=180> <weight=3 > >                   		 

 )";

	_results.ResetLayout(70);
	_targets.ResetLayout(70);
	_nTsec.ResetLayout(70);
	_PCRfiltre.ResetLayout(70);
	_PrimersFilePCR.ResetLayout(70);
	_Prob_uArr.ResetLayout(70);
	_NNParamFile.ResetLayout(70);

	numUpDowSdConc.ResetLayout(80);
	numUpDowTa.ResetLayout(90);
	numUpDowTgConc.ResetLayout(80);
	numUpDowSalConc.ResetLayout(110);

	_gr_dir.div(" vert "
		"  <weight=26  Results  >       \n\t"
		"  <min=280 margin=5 seq>    	    \n\t"
		"  <weight=26  <NN_param><weight=80 ckBx_loadNNParam> >       \n\t");

	_gr_seq.div(" vert "
		"  <weight=62  margin=[0,5,0,5] _targets       >    		         \n\t"
		"  <weight=62  margin=[0,5,0,5] _nTsec         >    		         \n\t"
		"  <weight=50  margin=[0,5,0,5] _PCRfiltre     >    		         \n\t"
		"  <weight=62  margin=[0,5,0,5] _PrimersFilePCR>    		         \n\t"
		"  <weight=67  margin=[0,5,5,5] _Prob_uArr     >    		         \n\t");

	_gr_targ.div(" vert "
		"<weight=22 dir>                                           \n\t"
		"<weight=23 gap=10 <weight=10%>< Opt ><weight=10%>    >   \n\t"
		"<>                                                        \n\t");

	_gr_ntarg.div(" vert "
		"<weight=22 dir>                                           \n\t"
		"<weight=23 gap=10 <weight=10%>< Opt ><weight=10%>    >   \n\t"
		"<>                                                        \n\t");

	_gr_PCRfiltre.div(" vert "
		"<weight=22 dir>                                           \n\t"
		"<weight=23 gap=10 <weight=10%>< Opt ><weight=10%>    >   \n\t"
		"<>                                                        \n\t");

	_gr_PrimersFilePCR.div(" vert "
		"<weight=22 dir>                                           \n\t"
		"<weight=23 gap=10 <weight=10%>< Opt ><weight=10%>    >   \n\t"
		"<>                                                        \n\t");

	_gr_uArr.div(" vert "
		"<weight=22 dir>                                           \n\t"
		"<weight=23 gap=10 <weight=10%>< Opt ><weight=10%>    >   \n\t"
		"<>                                                        \n\t");

	_gr_checks.div(" <vertical weight=210 checks>     				\n\t");

	_gr_salt.div(
		"     < <> 		                                    \n\t"
		"                 <weight=200 vertical ConcST        gap=2>     \n\t"
		"                 <> 		                                    \n\t"
		"		          <weight=230 vertical ConcSaltTa    gap=2>     \n\t"
		"                 <> 		                                    \n\t"
		"		          <weight=250 vertical gap=5  <weight=23 SMeth >        	\n\t"
		"				                              <weight=23 AMeth >  >       \n\t"
		"                 <>  		 >                                  \n\t");


}

void  SetupPage::AsignWidgetToFields ()  
    {
      _setup<< link( _Pr._cp._OutputFile      ,       _results  )
            << link( _Pr._cp._InputTargetFile ,       _targets  )
            << link( _Pr._cp._TRecurDir         ,     _chkTargRecDir)
            << link( _Pr._cp._TDirStrOnly      ,     _chkTargOnlyStruct)
            << link( _Pr._cp._NonTargetFile     ,       _nTsec  )
            << link( _Pr._cp._nTRecurDir      ,     _chk_nTgRecDir)
            << link( _Pr._cp._nTDirStrOnly      ,     _chk_nTgOnlyStruct)
            << link( _Pr._cp._PCRfiltrPrFile  ,       _PCRfiltre)
            << link( _Pr._mPCR._InputSondeFile , _PrimersFilePCR)            
            << link( _Pr._mPCR._PrRecurDir      ,     _chkPrimRecDir)
            << link( _Pr._mPCR._PrDirStrOnly      ,     _chkPrOnlyStruct)
            << link( _Pr._uArr._InputSondeFile , _Prob_uArr)            
            << link( _Pr._uArr._PrRecurDir      ,     _chkProbRecDir)
            << link( _Pr._uArr._PrDirStrOnly      ,     _chkProbOnlyStruct)
            << link( _Pr._cp._InputNNFile       , _NNParamFile  )
            << link( _Pr._cp.ConcSd	    ,       numUpDowSdConc  )
            << link( _Pr._cp.ConcSalt	     , numUpDowSalConc  )
            << link( _Pr._cp.ConcTg	    ,       numUpDowTgConc  )
            << link( _Pr._cp.Ta	            ,       numUpDowTa  )        
            << link( _Pr._cp.SaltCorr	  ,      comBoxSalMeth  )        
            << link( _Pr._cp.TAMeth       ,       comBoxTAMeth  )        
            << link( _Pr._cp.st_savTm	  ,       ckBx_savTm    )        
            << link( _Pr._cp.st_savPos	  ,       ckBx_savPos   )        
            << link( _Pr._cp.st_savG 	  ,         ckBx_savG   )        
            << link( _Pr._cp.st_savAlign  ,     ckBx_savAlign   )        
            << link( _Pr._cp.st_savProj   ,     ckBx_savProj    )        
            << link( _Pr._cp.st_savG_Plasm , ckBx_savG_Plasm    )        
            << link( _Pr._cp.st_savTm_Plasm , ckBx_savTm_Plasm  )     
            << link( _Pr._cp.st_savLog 	  ,       ckBx_savLog   )        
            << link( _Pr._cp.st_Exp_sond  ,ckBx_savExportSond   )        
            << link( _Pr._cp.st_ExpTarg	 ,ckBx_savExportTarg    )        
            << link( _Pr._cp.saveNNPar     , ckBx_savNNParam    )        
            << link( _Pr._cp.loadNNPar    , ckBx_loadNNParam    )        
          ;
            
        _place.field("Project"  )    <<  _proj      ;
        _place.field("dir"      )    <<  _gr_dir    ;
	    _place.field("buttons"  )    <<  _set_def_proj << _load_def_proj;
	    _place.field("checks"   )    <<  _gr_checks ;
	    _place.field("salt"     )    <<  _gr_salt   ;

        _gr_dir ["Results"  ]    <<  _results   ;
        _gr_dir ["seq"      ]    <<  _gr_seq    ;
	    _gr_dir ["NN_param" ]    << _NNParamFile  ;
	    _gr_dir ["ckBx_loadNNParam"]    <<   ckBx_loadNNParam ;
        _gr_seq ["_targets"        ]    <<   _gr_targ  ;
        _gr_seq ["_nTsec"          ]    <<   _gr_ntarg  ;
        _gr_seq ["_PCRfiltre"      ]    <<   _gr_PCRfiltre  ;
        _gr_seq ["_PrimersFilePCR" ]    <<   _gr_PrimersFilePCR  ;
        _gr_seq ["_Prob_uArr"      ]    <<   _gr_uArr  ;

        _gr_targ ["dir"  ]    <<   _targets  ;
        _gr_targ ["Opt"  ]    <<   _chkTargRecDir  << _chkTargOnlyStruct ;

        _gr_ntarg ["dir"  ]    <<   _nTsec  ;
        _gr_ntarg ["Opt"  ]    <<   _chk_nTgRecDir <<  _chk_nTgOnlyStruct          ;

        _gr_PCRfiltre ["dir"  ]    <<   _PCRfiltre  ;
        //_gr_PCRfiltre .plc["Opt"  ]    << _chkTargRecDir  << _chkTargOnlyStruct ;

        _gr_PrimersFilePCR ["dir"  ]    <<   _PrimersFilePCR  ;
        _gr_PrimersFilePCR ["Opt"  ]    <<   _chkPrimRecDir  << _chkPrOnlyStruct ;

        _gr_uArr ["dir"  ]    <<   _Prob_uArr  ;
        _gr_uArr ["Opt"  ]    <<   _chkProbRecDir  << _chkProbOnlyStruct ;

	    _gr_checks["checks"  ]   << ckBx_savTm    << ckBx_savPos     <<ckBx_savG         << ckBx_savAlign 
                                     << ckBx_savProj  << ckBx_savG_Plasm << ckBx_savTm_Plasm << ckBx_savLog
                                     << ckBx_savExportSond << ckBx_savExportTarg<< ckBx_savNNParam;

	    _gr_salt["ConcST"     ]   << numUpDowSdConc       << numUpDowTgConc ;
	    _gr_salt["ConcSaltTa" ]   << numUpDowSalConc      << numUpDowTa ;
	    _gr_salt["SMeth"      ]   << " Salt Correct. Method:"	   <<  comBoxSalMeth;
	    _gr_salt["AMeth"      ]   << " ThDy Align. Method"       <<  comBoxTAMeth ;

		_place.field("Firma") << _firma;

    }
void  SetupPage::MakeResponive()
    {
        _proj.add_filter(("ThDy project"),("*.ThDy.txt"));
        _proj.onOpenAndSelectFile ([this](const std::string &file)
        { 
            this->LoadProject (  file  );
                std::cerr << "onOpenAndSelectFile: Loaded Project file: " << file << std::endl;
        } );
        _proj.onSaveFile          ([this](const std::string &file)
        { 
            this->_Pr.save (nana::charset ( (  file  ))); 
                std::cerr << "onSaveFile: Saved Project file: " << file << std::endl;
        } );

        AddFastaFiltre(_targets );
        AddFastaFiltre(_nTsec );
        AddFastaFiltre(_PCRfiltre );
        AddFastaFiltre(_PrimersFilePCR );

        _NNParamFile.add_filter      (("Nearest Neibrhud Parametr"),("*.NN.csv"));
        _NNParamFile.onOpenAndSelect ([this]()
            { 

            //assert((  std::cerr<< "\nBefore loading NNfile, SetupPage: " , true  ));;
            //assert((  std::wcerr<< caption() << std::endl , true  ));;
            
            std::ifstream nn(_Pr._cp._InputNNFile.get());
                _Pr._cp._pSaltCorrNNp->LoadNNParam(nn ) ;

            //assert((  std::cerr << "onOpenAndSelect: Opened NNp file: " << _Pr._cp._InputNNFile.get() << std::endl, true  ));;
                return true;

            } );
        _NNParamFile.onSave          ([this]()
        {

            //assert((  std::cerr<< "\nBefore saving NNfile, SetupPage: " , true  ));;
            //assert((  std::wcerr<< caption() << std::endl , true  ));;
            
            std::ofstream{ _Pr._cp._InputNNFile.get() /* + (".NN.csv") */} << *_Pr._cp._pSaltCorrNNp;        

            //assert((  std::cerr << "onSave: Saved NNp file: " << _Pr._cp._InputNNFile.get() << std::endl, true  ));;
         } );

        _set_def_proj .events().click ([&](){ setAsDefProject() ;} );
        _load_def_proj.events().click ([&](){ RestDefPr      () ;} );
    }
void  SetupPage::SaveProj()
	{	 
        if(  _proj.Canceled () )  return;
        _Pr.save ( _proj.FileName());
	}
void  SetupPage::setAsDefProject()
    {
		std::string caption = "Set current setting as Default project";
		std::string message = std::string("This will overwrite the current Default Project." )        + "\n\n"
						+  "Are you sure?"   + "\n\n"
						+ "\tYes:  The default project will be overwrited. " + "\n"
						+ "\tNo:  No action will be taken. " + "\n"
						;
		switch ( (nana::msgbox(  *this, nana::charset (caption) , nana::msgbox::yes_no )
                        <<  message
                    ).icon(nana::msgbox::icon_question ) .show (  ))
		{
			case  nana::msgbox::pick_yes :  
                                    _Pr.save_asDefPr() ; 					 // crea el Def Project.
				return;

			case  nana::msgbox::pick_no:    
            default:;
        }
    }
void  SetupPage::RestDefPr	 ( )		// Restore (USE) Deff  Proj File
{		 
    try{
		    _Pr.load_defPr() ;			// cuando no existe Def Project: 1ra vez que se usa el prog??
		}
	catch ( std::exception& e)
	{ 
		(nana::msgbox ( ("Error loading Def Project" ) )<< e.what()).show() ;
 	}		 
}
void  SetupPage::AddMenuItems(nana::menu& menu)
{
    menu.append(("New"    )  , [&](nana::menu::item_proxy& ip)  {  ;  } );// ??
    menu.append(("Open...")  , [&](nana::menu::item_proxy& ip)  
    { 
        _proj.open(_proj.FileName());
        if (!_proj.Canceled())
            LoadProject( _proj.FileName());
    } );
    menu.append(("Save...")  , [&](nana::menu::item_proxy& ip)  
    { 
        _proj.save(_proj.FileName()); 
        if( ! _proj.Canceled () )   
            _Pr.save ( _proj.FileName());
    } );
    menu.append_splitter();
    menu.append(("Set as default") , [&](nana::menu::item_proxy& ip)  {;  });
    menu.append(("Restore default"), [&](nana::menu::item_proxy& ip)  {;  });
    menu.append_splitter();
    menu.append(("Exit"    )  , [&](nana::menu::item_proxy& ip)  {  _Pr.close();  } );

}
void  SetupPage::LoadProjectAndReg(std::string file)
    {
        LoadProject(file);
 		_proj.FileName(file  );
   }
void  SetupPage::LoadProject(std::string file)
	{
		try
		{
			_Pr.load( file );
 		}
		catch (std::exception& e)
		{
			std::string caption = "Error trying to load the project file:";
			std::string message =   file          + "\n\n"
							+  e.what()   + "\n\n"
							+  "Use the Default project?"   + "\n\n"
							+ "\tYes:  The default project file will be loaded. " + "\n"
							+ "\tNo:  Select a project file to be loaded. " + "\n"
							+ "\tCancel: Use the values correctly loaded mixed with the\t\t\t previous existing. "
							;
			switch ( (nana::msgbox(  *this, nana::charset (caption) , nana::msgbox::yes_no_cancel )
                            <<  message
                        ).icon(nana::msgbox::icon_error) .show (  ))
			{
				case  nana::msgbox::pick_yes :  
					    _Pr.load_defPr();
                        _proj.FileName(_Pr.ProjectFile ()  );
					return;

				case  nana::msgbox::pick_no:    
                        _proj.open (nana::charset (file));
                        if ( !  _proj.Canceled() )
                                LoadProject( _proj.FileName());
                        return;
			}
		}
	}
   
SetupPage::SetupPage          (ThDyNanaForm& tdForm) try
        : _Pr           (tdForm), 
          CompoWidget   (tdForm, ("Setup"), ("Setup.lay.txt")),
		  _firma{*this, _Pr .e_mail_firma}
    {
		_results.folder = true;
		_targets.folder = true;
		_nTsec.folder   = true;

	    InitMyLayout();
        SelectClickableWidget( _set_def_proj);
        SelectClickableWidget( *this);

        MakeResponive();
    }
   catch (std::exception & e)
   {
	   throw std::runtime_error(std::string("An error occurred during initialization of the Setup page window:\n") + e.what());
   }
   catch (...)
   {
	   throw std::runtime_error(std::string("An unknown error occurred during initialization of the Setup page window"));
   }



FilePickBox::filtres SetupPage::FastaFiltre( )
    {
        return FilePickBox::filtres       { {("fasta")       , ("*.fas;*.fasta"     ) },
                                            {("NCBI BLAST")  , ("*-Alignment.xml"   ) },
                                            {("GB"        )  , ("*.gb;*-sequence.xml")},
                                            {("Text"      )  , ("*.txt"             ) },
                                            {("All sequences"), ("*.fas;*.fasta;*.txt;*-Alignment.xml;*.gb;*-sequence.xml")},
                                            {("All files" )  , ("*.*"               ) }
                                          } ;
    }


