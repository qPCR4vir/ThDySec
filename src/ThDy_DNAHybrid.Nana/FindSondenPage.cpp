/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\FindSondenPage.cpp
*
* @brief 
*/

#include "ThDy_DNAHybrid.Nana\FindSondenPage.h"
#include "ThDy_DNAHybrid.Nana\main.Nana.h"

int      n_len{ 6 }, n_dec{ 1 };

std::string str(double val, int ln=n_len, int dc=n_dec) 
{
	std::string s(n_len, 0);

	auto l = snprintf(&s[0], n_len + 1, ("% *.*f"), n_len, n_dec, val );

	return s;
}

std::string str(long val, int ln=n_len) 
{
	std::string s(n_len, 0);

	auto l = snprintf(&s[0], n_len + 1, ("% *d"), n_len,  val );

	return s;
}


class TableCandRes : public nana::form, public EditableForm
{
	//CProgParam_SondeDesign::targets_comp  targets_comparitions;
	nana::listbox                         list{*this};
  public:

	TableCandRes ( std::vector < CProgParam_SondeDesign::targets_comp > &&targets_comparitions,
		           std::string name  ) 
		:   /*targets_comparitions ( std::move(targets_comparitions) ),*/
			nana::form(nana::rectangle(nana::point(50, 5), nana::size(1000, 650))),
			EditableForm(nullptr, *this, std::string("Target comparison & candidate probes: ") + name, "TableCand.lay.txt" )
	{ 
		auto& sch = list.scheme();
		sch.header_height = 20;
		sch.text_margin = 2;
		sch.item_height_ex = 1;  ///< Set !=0 !!!!  def=6. item_height = text_height + item_height_ex
		//sch.item_height = sch.text_height + sch.item_height_ex;
		sch.header_splitter_area_before     = 4;
		sch.header_splitter_area_after      = 4;

		int i = 0;
		list.append_header(   "T Pos",        60);
		list.append_header(   "T Cand",       60);

		list.append_header(   "N Pos",        60);
		list.append_header(   "N Cand",       60);
		list.append_header(   "Targ #",       50);
		list.append_header(   "Targ name",        100);
		list.append_header(   "N Cand",       60);
		list.append_header(   "N Pos",        60);

		list.append_header(   "Num T Hits",       80);
		list.append_header(   "Num Hits OK",      80);

		list.append_header(   "N Pos",        60);
		list.append_header(   "N Cand",       60);
		list.append_header(   "Iterat#",           80);
		list.append_header(   "Targ name",         100);
		list.append_header(   "N Cand",       60);
		list.append_header(   "N Pos",        60);

		list.append_header(   "T Pos",        60);
		list.append_header(   "T Cand",       60);

		auto value_translator = [](const std::vector<nana::listbox::cell>& cells)
		{
			CProgParam_SondeDesign::targets_comp p;
			int i = 0;

			p.before.t_n_pos				 = std::stol(cells[i++].text);      //   "Num T Pos", 
			p.before.t_n_cand				 = std::stol(cells[i++].text);		//   "Num T Cand", 
							     
			p.before.target_1_n_cand_pos	 = std::stol(cells[i++].text);		//   "Num Pos", 
			p.before.target_1_n_cand		 = std::stol(cells[i++].text); 		//   "Num Cand", 
			p.target_num					 = std::stol(cells[i++].text);		//   "Targ Num", 
			p.target_1_name                  = cells[i++].text ;			    //   "Targ name", 
			p.after.target_1_n_cand          = std::stol(cells[i++].text); 	    //   "Num Cand", 
			p.after.target_1_n_cand_pos      = std::stol(cells[i++].text); 		//   "Num Pos", 
			 			     
			p.THits                          = std::stol(cells[i++].text);					//   "Num T Hits",      ????????????????????
			p.HitsOK                         = std::stol(cells[i++].text);					//   "Num Hits OK",     ????????????????????
			 		     
			p.before.target_2_n_cand_pos     = std::stol(cells[i++].text);		//   "Num Pos", 
			p.before.target_2_n_cand         = std::stol(cells[i++].text); 		//   "Num Cand", 
			p.iteration_num                  = std::stol(cells[i++].text);		//   "Iterat#", 
			p.target_2_name                  = cells[i++].text;			        //   "Targ name", 
			p.after.target_2_n_cand          = std::stol(cells[i++].text); 	    //   "Num Cand", 
			p.after.target_2_n_cand_pos      = std::stol(cells[i++].text); 		//   "Num Pos", 

			p.after.t_n_pos                  = std::stol(cells[i++].text);      //   "Num T Pos", 
			p.after.t_n_cand                 = std::stol(cells[i++].text);		//   "Num T Cand", 
			
			return p;
		};

		auto cell_translator = [](const CProgParam_SondeDesign::targets_comp& p)
		{
			std::vector<nana::listbox::cell> cells;

			cells.emplace_back( str (p.before.t_n_pos				) );
			cells.emplace_back( str (p.before.t_n_cand				) );
							        
			cells.emplace_back( str (p.before.target_1_n_cand_pos	) );
			cells.emplace_back( str (p.before.target_1_n_cand		) );
			cells.emplace_back( str (p.target_num					) );
			cells.emplace_back(     (p.target_1_name               )   );
			cells.emplace_back( str (p.after.target_1_n_cand       )   );
			cells.emplace_back( str (p.after.target_1_n_cand_pos   )   );
							        
			cells.emplace_back( str (p.THits                       )   );
			cells.emplace_back( str (p.HitsOK                      )   );
								    
			cells.emplace_back( str (p.before.target_2_n_cand_pos  )   );
			cells.emplace_back( str (p.before.target_2_n_cand      )   );
			cells.emplace_back( str (p.iteration_num               )   );
			cells.emplace_back(     (p.target_2_name               )   );
			cells.emplace_back( str (p.after.target_2_n_cand       )   );
			cells.emplace_back( str (p.after.target_2_n_cand_pos   )   );
							      
			cells.emplace_back( str (p.after.t_n_pos               )   );
			cells.emplace_back( str (p.after.t_n_cand              )   );

			return cells;
		};

		list.at(0).model<std::recursive_mutex>(std::move(targets_comparitions), value_translator, cell_translator);

		InitMyLayout();
		SelectClickableWidget(list);
		SelectClickableWidget(*this);


	}



	void SetDefLayout() override
	{
		_DefLayout =
			"table"
			;
	}

	void AsignWidgetToFields() override
	{
		_place["table"] << list;
	}

};

FindSondenPage::FindSondenPage(ThDyNanaForm& tdForm)    try
        : _Pr        (tdForm), 
          CompoWidget(tdForm, "Find probes", "FindSonden.lay.txt")
	{
		bgcolor (static_cast<nana::color_rgb>(0xAAAAAA));  ///\todo: use code

		chkBx_showFindedProbes.check(true);
		InitMyLayout();
		SelectClickableWidget( *this);

		_design .events().click([&]() 
		{
			Run_Design(true );  
		});    

		_compare.events().click([&]() 
		{
			Run_Design(false);  
		});  

		//_Gmin.tooltip(("Only probes with stronger interaction with target (smaller G by selected Ta) will be \"include\""));
	}
catch (std::exception & e)
{
	throw std::runtime_error(std::string("An error occurred during initialization of the Find probes page window:\n") + e.what());
}
catch (...)
{
	throw std::runtime_error(std::string("An unknown error occurred during initialization of the Find probes page window"));
}

void FindSondenPage::SetDefLayout()
{
	_DefLayout = 

R"( 
	vertical   gap=2    margin=5               
			<weight=10     >       			 
	 	    <weight=260 gap=8 <weight=5> <weight=388 vertical   		                       
		                                                   <weight=115 <weight=388 Sonde  > >		  
	 		                                               <weight=10>			 
	 			                                           <weight=72 TargCov        >    		 
	 		                                               <weight=10> 			 
	 			                                           <weight=40 <   <> <weight=300   gap=20 Run>  <>     > >    			 
	 			                                           <weight=10>                                 			 
	 			                         >   <> <weight=230 gap=1 vertical  options> <weight=5>      
	 	    >   			 
	 		<weight=23   <weight=140> <Output>   <> >       		          
	        <>				
   		    < weight=21 <><Firma weight=180> <weight=3 > >                   		 
 )";


	_gr_probes.div("vert < Sonde  margin=2 gap= 2 grid=[2,4]  	    \n\t"
		"					                                    		\n\t"
		"						                                  >	\n\t");

	_gr_prob_tg.div("<  margin=2 gap= 2 vertical   options>");
	_gr_prob_ntg.div("<  margin=2 gap= 2 vertical   options>");
	_gr_probself.div("<  margin=2 gap= 2 vertical   options>");
	_gr_find_prb.div("<  margin=5 gap= 2 TargCov grid=[2,2]>");

	_Gmin.ResetLayout(45, 40, 55);   _Gmax.ResetLayout(1, 40, 75);
	_Tmmin.ResetLayout(45, 40, 55);  _Tmmax.ResetLayout(1, 40, 75);
	_Lengthmin.ResetLayout(45, 40, 55);  _Lengthmax.ResetLayout(1, 40, 75);

	_MaxG.ResetLayout(110, 45, 50);
	_MinTm.ResetLayout(110, 45, 50);

	_MinG.ResetLayout(110, 45, 50);
	_MaxTm.ResetLayout(110, 45, 50);

	_MinSelfG.ResetLayout(110, 45, 50);
	_MaxSelfTm.ResetLayout(110, 45, 50);

	numUpDw_MinTargCov.ResetLayout(30, 40, 40);
	numUpDw_MaxTargCov.ResetLayout(30, 40, 40);
}





void FindSondenPage::AsignWidgetToFields()
{
	using ParamGUIBind::link;

	_findSond << link(_Pr._SdDes.G_sig, _MaxG)
		<< link(_Pr._SdDes.Tm_sig, _MinTm)
		<< link(_Pr._SdDes.MinSd_nTgG, _MinG)
		<< link(_Pr._SdDes.MaxSd_nTgTm, _MaxTm)
		<< link(_Pr._SdDes.MinSelfG, _MinSelfG)
		<< link(_Pr._SdDes.MaxSelfTm, _MaxSelfTm)
		<< link(_Pr._SdDes.sL.G, _Gmin, _Gmax)
		<< link(_Pr._SdDes.sL.T, _Tmmin, _Tmmax)
		<< link(_Pr._SdDes.sL.L, _Lengthmin, _Lengthmax)
		<< link(_Pr._SdDes.common, chkBx_common)
		<< link(_Pr._SdDes.unique, chkBx_unique)
		<< link(_Pr._SdDes.Coverage, numUpDw_MinTargCov, numUpDw_MaxTargCov)

		;

	/// Use room (wd,w,h) in combination with a <Table grid=[W,H]>
	_place["Sonde"] << _gr_probes;
	_place["TargCov"] << _gr_find_prb;
	_place.field("Run") << _design << _compare;
	_place.field("options") << _gr_prob_tg << _gr_prob_ntg << _gr_probself;
	_place.field("Output") << chkBx_showFindedProbes;

	_gr_probes["Sonde"] << ("                               Min.") << ("           Max.")
		<< _Gmin << _Gmax
		<< _Tmmin << _Tmmax
		<< _Lengthmin << _Lengthmax;

	_gr_find_prb["TargCov"] << chkBx_unique << numUpDw_MinTargCov
		<< chkBx_common << numUpDw_MaxTargCov;

	_gr_prob_tg["options"] << _MaxG << _MinTm;
	_gr_prob_ntg["options"] << _MinG << _MaxTm;
	_gr_probself["options"] << _MinSelfG << _MaxSelfTm;

	_place.field("Firma") << e_mail_firma;

}

void FindSondenPage::Run_Design(bool design)
{
    _Pr._SdDes._design	 = design ;		
		 
	try{                                   
			_Pr._SdDes.probes=  _Pr._mPCR._probesMS.get();  /// Use _Pr._SdDes.probes to attach the results

            _Pr.Run(_Pr._SdDes);	                        ///  Why not   _Pr._SdDes.Run ();? \todo make paralel with nana::progres

                if (chkBx_showFindedProbes.checked()) 
                ( dynamic_cast<ThDyNanaForm&>(_Pr)).mExpl_.ShowFindedProbes_in_mPCR();

				for (auto & p : _Pr._SdDes.targets_comparitions)
					std::cout << "\n" << p.iteration_num << ". " << p.target_1_name << " vs. " << p.target_2_name;

				_Pr._results.emplace_back(new TableCandRes( std::move(_Pr._SdDes.targets_comparitions ), _Pr._cp._OutputFile.get() ) );
				_Pr._results.back()->show();

 	}
	catch ( std::exception& e)
	{ 
        (nana::msgbox(*this,"Error during probe Design !", nana::msgbox::button_t::ok)<<e.what()) (  ) ;
		return;
	}	 	        		 
}   

