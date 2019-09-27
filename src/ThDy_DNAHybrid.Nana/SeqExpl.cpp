/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2019
*
* @file  ThDySec/src/ThDy_DNAHybrid.Nana/SeqExpl.cpp
*
* @brief 
*
*/

#include <algorithm>

#include "ThDy_DNAHybrid.Nana/SeqExpl.h"
#include "ThDy_DNAHybrid.Nana/main.Nana.h"

namespace fs = std::filesystem ;

SeqExpl::SeqExpl  (ThDyNanaForm& tdForm)
        : _Pr             (tdForm), 
          CompoWidget     (tdForm, ("Seq Explorer"), ("SeqExpl.lay.txt")),
          _firma{*this, _Pr .e_mail_firma}
    {
        auto& sch = _list.scheme();
		sch.header_padding_bottom = sch.header_padding_top = 1;
        sch.text_margin = 2;
        sch.item_height_ex = 1;
        sch.header_splitter_area_before = 4;
        sch.header_splitter_area_after = 4 ;

        _statusbar.format(true).bgcolor(nana::colors::turquoise );
        
        InitMyLayout();
        SelectClickableWidget( _tree);
        SelectClickableWidget( *this);
        _tree.checkable(true);
        _list.checkable(true);
        _list.append_header(("Name")  , 120);     // col 0: name  
        _list.append_header(("Length"), 50);      // col 1: len
        _list.append_header(std::string("Tm ")+RTunits::grC , 60);      //case 2: Tm 
        _list.append_header(("Deg")   , 50);      // case 3: deg   
        _list.append_header(("Description")   , 220);   // case 4: descr  
        _list.append_header(("Beg"), 50);         // case 5: beg in aln 
        _list.append_header(("End"), 50);         // case 6: end in aln    
        _list.append_header(("Seq")   , 420);

        AddMenuItems(_menuProgram);
        MakeResponive();

        _tree.auto_draw(true);
        _list.auto_draw(true);
    }

SeqExpl::Node SeqExpl::AddNewSeqGr  (Node node)
{
    try{
            return appendNewNode(node, _Pr._cp.AddSeqGroup(**node.value<MSecIt>(),"New group")).expand(true);
        }
    catch ( std::exception& e)
        {
          (nana::msgbox ( ("Error adding new group" ) )<< e.what()).show() ;
          return node;
        }
}

SeqExpl::Node SeqExpl::AddMSeqFiles (const std::filesystem::path &file, bool  all_in_dir)
    {     
    try{ 
            Node   tn = _tree.selected();
			if (tn.empty()) return tn; // trhow

            CMultSec::pMSec ms = *tn.value<MSecIt>();
            _Pr._cp.AddSeqFromFile(*ms, file, all_in_dir);
            return Refresh(tn);
        }
        catch ( std::exception& e)
        { 
            (nana::msgbox ( ("Error adding sequences" ) )<< e.what()).show() ;
			Node   tn = _tree.selected();
			if (tn.empty()) return tn; //?
            return tn;
         }         
    }

void SeqExpl::AddMenuItems(nana::menu& menu)
    {
        menu.append_splitter();
        using Mitem = nana::menu::item_proxy;
        menu.append(("Add a new, empty, group for sequences"      ),[&](Mitem& ip) 
			{	Node   tn = _tree.selected();
				if (tn.empty()) return; //?
				AddNewSeqGr(tn); 
			});
        menu.append(("Add a group of sequences from a file..."    ),[&](Mitem& ip) { Click(_loadFile);              });
        menu.append(("Add groups of sequences from a directory..."),[&](Mitem& ip) { Click(_loadDir);               });

        menu.append_splitter();

        menu.append(("Scan the structure of directory..."),[&](Mitem& ip) { Click(_scanDir);              });
        menu.append(("Re-scan the original directory"   ), [&](Mitem& ip) { Click(_re_scanDir);           });
        menu.append(("Reload the original file"         ), [&](Mitem& ip) 
			{	Node   tn = _tree.selected();
				if (tn.empty()) return; //?
				ReloadFile(tn); 
			});
        menu.append(("Reload the original directory"    ), [&](Mitem& ip) 
			{	Node   tn = _tree.selected();
				if (tn.empty()) return; //?
				ReloadDir (tn); 
			});
        menu.append(("Replace from a file . . ."        ), [&](Mitem& ip)
        {
            auto tn= _tree.selected();
			if (tn.empty()) return; //?
			MSecIt ms = tn.value<MSecIt>();
            fs::path pt{(*ms)->_orig_file};
            pt= fs::canonical(pt).make_preferred();
            nana::filebox  fb{ *this, true };
            fb .add_filter ( SetupPage::FastaFiltre( )   )
                    .init_file(pt.string())
                    .title(("Replace a group of sequences from a file"));
            auto path=fb.show();
            if ( path.empty() ) return;

            Replace(tn, ms, path[0], false);
        });

        menu.append(("Replace from directory . . ."), [&](Mitem& ip)
        {
            Node   tn = _tree.selected();
			if (tn.empty()) return; //?
			MSecIt ms = tn.value<MSecIt>();

            nana::folderbox  fb{ *this, (*ms)->_orig_file.parent_path(),
                                 "Replace/reload a group of sequences from a directory" };
            auto p=fb.show();
            if (p.empty()) return;

            Replace(tn, ms, p[0], true);
        });

        menu.append_splitter();

        menu.append     ( ("Show Only local sequences"),[&](Mitem& ip) { ShowLocals( menu.checked(ip.index())); })
            .check_style( nana::menu::checks::option)
            .checked    ( false );

        menu.append     ( ("Show filtered sequences"  ),[&](Mitem& ip) { ShowFiltered( menu.checked(ip.index())); })
            .check_style( nana::menu::checks::highlight )
            .checked    ( true );
 
        menu.append_splitter();
        menu.append(("Cut selected sequences from list"          ),[&](Mitem& ip)  {  Click(_cutSec);  });
        menu.append(("Cut selected groups of sequences from tree"),[&](Mitem& ip)  {  Click(_cut   );  });
        menu.append(("Paste the sequences"                       ),[&](Mitem& ip)  {  Click(_paste );  });

        menu.append_splitter();
        menu.append(("Del selected sequences from list"          ),[&](Mitem& ip)  {  Click(_delSec);  });
        menu.append(("Del selected groups of sequences from tree"),[&](Mitem& ip)  {  Click(_del   );  });
        menu.append(("Rename the selected group of sequences"    ),[&](Mitem& ip)
        {
            
			Node   tn = _tree.selected();
			if (tn.empty()) return; //?
			RenameFrom rnm(_tree, tn.text());
            nana::API::modal_window( rnm );
            tn.text(rnm.Name());
            tn.value<CMultSec*>()->_name = rnm.Name() ;

        }).enabled(true);

        menu.append_splitter();
        auto  indxFASTA = menu.append(("Export FASTA . . ."),[&](Mitem& ip) {  /*_tree.selected().value<CMultSec*>()->ExportFASTA();*/  }).index();
        auto& menuFASTA = *menu.create_sub_menu(indxFASTA);
        menuFASTA.append(("Only current sequences"     ),[&](Mitem& ip){  ;  });
        menuFASTA.append(("Selected sequences in group"),[&](Mitem& ip)
			{	Node   tn = _tree.selected();
				if (tn.empty()) return; //?
				_Pr.ExportFASTA(tn.value<MSecIt>()->get(), true );  
			});
        menuFASTA.append(("All sequences in group"     ),[&](Mitem& ip)
			{	Node   tn = _tree.selected();
				if (tn.empty()) return; //?
				_Pr.ExportFASTA(tn.value<MSecIt>()->get(), false);  
			});
        menuFASTA.append(("All selected sequences"     ),[&](Mitem& ip){ _Pr._cp._pSeqTargets->Export_as("export.fasta", true )  ;  });
        menuFASTA.append(("All sequences"              ),[&](Mitem& ip){ _Pr._cp._pSeqTargets->Export_as("export.fasta", false)  ;  });

    }

void SeqExpl::MakeResponive()
    {
        // the "selected" feature in the GUI have no efect in the data, it is a pure GUI feature,
        // but the check status go to the data with the name Selected() !!!

        _tree.events().selected ( [&]( const nana::arg_treebox &tbox_arg_info ) 
                                 { if (tbox_arg_info.operated) RefreshList(tbox_arg_info.item); });

        _tree.events().checked  ( [&]( const nana::arg_treebox &tbox_arg_info )
        {                                              
            (*(tbox_arg_info.item.value<MSecIt>()))->Selected(tbox_arg_info.operated);
            if (tbox_arg_info.item== _tree.selected())
                RefreshList(tbox_arg_info.item);                //  ??????? Only RefreschList
        });

        _list.events().checked  ( [&](  const nana::arg_listbox &lbox_arg_info )
        {                                               
            (*(lbox_arg_info.item.value<SecIt>()))->Selected(lbox_arg_info.item.checked());
        });
 
        _loadFile   .events().click([this]()                                  //  ------------  LOAD--------------
                        {
                            auto      tn    = _tree.selected();
							if (tn.empty()) return;
							MSecIt ms    = tn.value<MSecIt>();
                            fs::path pt= fs::canonical((*ms)->_orig_file).make_preferred();

                            nana::filebox  fb{ *this, true };
                            fb .add_filter ( SetupPage::FastaFiltre( )                   )
                                    .init_file(pt.string())
                                    .title      ( ("File load: Add a group of sequences from a file") );

                            if (auto path=fb(); ! path.empty())
                               AddMSeqFiles(path[0].string(), false);
                        });
        //_loadFileTT.set(_loadFile,("File load: Add a group of sequences from a file"));

                                                                              //  ------------  RE-LOAD-------------- ?
        _re_loadFile.events().click([this]()  
			{  
				auto      tn = _tree.selected();
				if (tn.empty()) return;
				ReloadFile(tn);    
			});

                                                                             //  ------------  LOAD-DIR --------- ?
        _loadDir    .tooltip(("Directory load: Add a tree of groups of sequences from a directory."))
                    .events().click([this]()
                        {
							auto      tn = _tree.selected();
							if (tn.empty()) return;
							MSecIt ms = tn.value<MSecIt>();
                            nana::folderbox  fb{ *this, (*ms)->_orig_file, "Directory load: Add a tree of groups of sequences from a directory" };
                            auto path=fb();
                            if (path.empty()) return;
                            AddMSeqFiles(path[0].string(), true);
                        });
                                                                             //  ------------  RE-LOAD-DIR --------- ?
        _re_loadDir .tooltip(("Directory reload: Reload a tree of groups of sequences from a directory,\npossibly using new filters."))
                    . events().click([this]()  
						{  auto      tn = _tree.selected();
		                   if (tn.empty()) return;
		                   ReloadDir (tn);    
						});

                                                                             //  ------------  RE_SCAN-DIR ---------
        _re_scanDir .tooltip(("Re-Scan the original directory: reproduce the structure."))
                .events().click([this]() 
					{ auto      tn = _tree.selected();
					  if (tn.empty()) return;
				      ReloadDir (tn, true); 
					});

                                                                             //  ------------  SCAN-DIR ---------
        _scanDir    .tooltip(("Directory scan: Reproduce the structure of directory..."))
                    .events().click([this]()
                        {
							auto      tn = _tree.selected();
							if (tn.empty()) return;
							MSecIt ms = tn.value<MSecIt>();
                            nana::folderbox  fb{ *this, (*ms)->_orig_file, "Directory scan: Reproduce the structure of directory..." };
                            auto path=fb();
                            if (path.empty()) return;
                            Replace(tn, ms, path[0], true, true);
                        });
                                                                                //  ------------  PASTE ---------
        _paste      .tooltip(("Paste sequences"))
                    .events().click([this]()
        {
            Node    tn = _tree.selected();
            if (tn.empty()) return;
			CMultSec *pms = tn.value<MSecIt>()->get();

            for (MSecIt ms : _dragMSec)
                pms->MoveMSec(ms);     /// \todo use MoveMSec   ?????!!!!!!
            for (SecIt s : _dragSec)
                pms->MoveSec(s);

            _dragMSec.clear();
            _dragSec .clear();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            populate(tn);
            populate(_tree.find(("Don t use") )); ///\todo better - this is error-prone
            _list.clear();
            populate_list_recur(tn);
            tn.select(true).expand(true);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });

        _cut        .tooltip(("Cut a group of sequences"))                    //  ------------  CUT Gr--------------
                    .events().click([this]()
        {
            Node tn= _tree.selected();
			if (tn.empty()) return;
			if (tn->owner().empty() || tn->owner()->owner().empty())    //   ???  if( tn->level() < 2 );  ///\todo use isRoot() ?
            {
                (nana::msgbox ( _tree , ("Cut a group of sequences " + tn->text()) )
                          << ("Sorry, you can�t cut the group: ") + tn->text() )
                          .icon(nana::msgbox::icon_error )
                          .show() ;
                return;
            }
            MSecIt ms = tn.value<MSecIt>();
            CMultSec *pms = (*ms)->_parentMS;
            assert(pms);

            _dragMSec.push_back(_Pr._cp._pSeqNoUsed->MoveMSec(ms));
            auto own = tn->owner();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            _tree.erase(tn);
            populate(appendNewNode (_tree.find(("Don t use") ), ms ));
            own.select(true).expand(true);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });

        _del        .tooltip(("Delete a group of sequences "))               //  ------------  DEL --------------
                    .events().click([this]()
        {
            auto tn= _tree.selected();
			if (tn.empty()) return;
            if (tn->owner().empty() || tn->owner()->owner().empty())
            {
                (nana::msgbox ( _tree , ("Deleting a group of sequences " + tn->text()) )
                          << ("Sorry, you can�t delete the group: ") + tn->text() )
                          .icon(nana::msgbox::icon_error )
                          .show() ;
                return;
            }
            MSecIt ms = tn.value<MSecIt>();
            CMultSec *pms = (*ms)->_parentMS;
            assert(pms);

            _Pr._cp._pSeqNoUsed->MoveMSec(ms);
            auto own = tn->owner();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            _tree.erase(tn);
            populate(appendNewNode (_tree.find(("Don t use") ), ms ));
            own.select(true).expand(true);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });

        _cutSec     .tooltip(("Cut selected sequences from list"))            //  ------------  CUT sec --------------
                    .events().click([this]()
        {
            auto sel =    _list.selected() ;
            // std::cout << "\nGoing to Cut selected sequences from list - items: " << sel.size();
            if (sel.empty()) return;
            _list.auto_draw(false);
            for (auto i : sel)
            {
                SecIt s=_list.at(i ).value<SecIt>();
                _dragSec.push_back(_Pr._cp._pSeqNoUsed->MoveSec(s));
            }
            _list.erase(sel);
            RefreshStatusInfo();
            _list.auto_draw(true);
        });

        _delSec     .tooltip(("Delete selected sequences from list"))          //  ------------  DEL sec --------------
                    .events().click([this]()
        {
            auto sel =    _list.selected() ;
            // std::cout << "\nGoing to delete selected sequences from list - items: " << sel.size();
            if (sel.empty()) return;
            _list.auto_draw(false);
            for (auto i : sel)
            {
                SecIt s=_list.at(i ).value<SecIt>();
                _Pr._cp._pSeqNoUsed->MoveSec(s);
            }
            _list.erase(sel);
            RefreshStatusInfo();
            _list.auto_draw(true);
        });

        _show_locals_s.enable_pushed(true)
                      .pushed(false)
                      .tooltip(("Show only local sequences, and not the sequences in internal trees"))   
                      .events().click([this]() { ShowLocals( _show_locals_s.pushed());  });

        _show_filt_s.enable_pushed(true)  
                    .pushed(true) 
                    .tooltip(("Show filtered sequences too"))   
                    .events().click([this]() { ShowFiltered( _show_filt_s.pushed());  });
    }

/// The 'old' node tn will be eliminated, and in his previous parent a new node will be created and returned,
/// The sequences in ms will be moved to DontUse and in his parent a new CMultiSec will be loaded from file.
SeqExpl::Node SeqExpl::Replace (Node tn, MSecIt ms, const fs::path& Path, bool all_in_dir, bool scan_only/*=false*/)
{        
try{ 
        if (isRoot(tn) || isRoot(tn->owner()) )
        {
            nana::msgbox ( ("Sorry, you can't replace group " + tn.text()) ).show() ;
            return tn;
        }
        Node      own = tn->owner();
        CMultSec *pms = (*ms)->_parentMS;
        assert(pms);
        if (pms != (*own.value<MSecIt>()).get())
            throw std::logic_error("The group of sequences to be replaced was not in the parent tree node.");

        MSecIt newms = _Pr._cp.AddSeqFromFile(*pms, Path, all_in_dir, scan_only);

        _Pr._cp._pSeqNoUsed->MoveMSec(ms);
        _tree.auto_draw(false);
        populate(_tree.find( _Pr._cp._pSeqNoUsed->_name));
        return Refresh(own) ;
    }
    catch ( std::exception& e)
    { 
        (nana::msgbox ( ("Error replacing sequences: " ) ).icon(nana::msgbox::icon_error)
            << "into group:    "  << tn.key()                                 
            << "\n from " << (all_in_dir?"directory: " : "file: ") << Path.string()     <<"\n"<< e.what()
        ).show() ;
     }        
    catch(...)
    {
            (nana::msgbox(("An uncaptured exception during replacing sequences: "))
                .icon(nana::msgbox::icon_error) 
            << "into "<< tn.key()                                
            << "from "<< (all_in_dir?"directory ":"file ") << Path.string()
            ).show();
    }
    return tn;
}

void SeqExpl::ShowFindedProbes_in_mPCR(bool show_/*=true*/)
{
    //auto idp = _Pr._cp.MaxTgId.get();
    //_Pr._cp.MaxTgId.set(100);
    //CMultSec *ms= _Pr._cp.AddSeqFromFile    (_Pr._mPCR._probesMS.get() , _Pr._cp._OutputFile.get() + ".sonden.fasta", false    );
    //_Pr._cp.MaxTgId.set(idp);

    RefreshProbes_mPCR( show_ );
}

void SeqExpl::RefreshProbes_mPCR(bool show_/*=true*/)
{
    Node probNode = _tree.find(_Pr._mPCR._probesMS->_name);
    Refresh(probNode).expand(true).select(true);
    if (show_) 
        _Pr.ShowExpl();
}

void SeqExpl::SetDefLayout()
{
    _DefLayout =
            "vertical                                                 	\n\t"
            "	          <height=23 <toolbar width=680 margin=2 ><>>    	\n\t"
            "	          <      <Tree  > |75% <List >   >                	\n\t"
            "	          <height=22 margin=[1,5,1,5] <width=80 tflab> <min=640 max=1000  TargetsOptions   > >  	\n\t"
            "	          <height=22 <statusbar margin=2  min=700 > <Firma max=180> <width=3 > >   	\n\t"
            "		\n\t"
        ;

    numUpDwMaxTgId.ResetLayout(60, 40, 30);
    numUpDw_TgBeg.ResetLayout(35, 40, 30);
    numUpDw_TgEnd.ResetLayout(35, 40, 30);
    numUpDw_SLenMin.ResetLayout(60, 40, 30);
    numUpDw_SLenMax.ResetLayout(60, 70, 30);
}

void SeqExpl::AsignWidgetToFields()
{
    using ParamGUIBind::link;

    _commPP << link(_Pr._cp.MaxTgId, numUpDwMaxTgId)
        << link(_Pr._cp.SecLim, numUpDw_TgBeg, numUpDw_TgEnd)
        << link(_Pr._cp.SecLenLim, numUpDw_SLenMin, numUpDw_SLenMax)
        ;


    _place["toolbar"] << "   Files:" << _loadFile << _re_loadFile << _paste
        << "      Dir:" << _loadDir << _re_loadDir << _scanDir << _re_scanDir << _cut << _del
        << "      Seq:" << _show_locals_s << _show_filt_s << _cutSec << _delSec
        ;
    _place["Tree"  ] << _tree;
    _place["List"  ] << _list;
    _place["tflab"] << "Targets filters: ";
    _place["TargetsOptions"]    << numUpDwMaxTgId
                                << numUpDw_TgBeg << numUpDw_TgEnd
                                << numUpDw_SLenMin << numUpDw_SLenMax;

    _place["statusbar"] << _statusbar;
    _place.field("Firma") << _firma;
}

std::string K2s(Temperature t) { return temperature_to_string(KtoC(t)); }

void SeqExpl::RefreshStatusInfo(const CMultSec &ms)
{
	//std::string local = 
    _statusbar.caption("  <bold>Local</> (sq: " + std::to_string(ms._Local._NSec)
        + ", gr: " + std::to_string(ms._Local._NMSec) + "),"
        + (ms._Local._NSec ? "    Length[" + std::to_string(ms._Local._Len.Min())
            + "-" + std::to_string(ms._Local._Len.Max()) + "]"
            + ", Tm[ " + K2s(ms._Local._Tm.Min()) + "-"
            + K2s(ms._Local._Tm.Max()) + "]  "
            : "")
        + ".        <bold>Global</> (sq: " + std::to_string(ms._Global._NSec)
        + ", gr: " + std::to_string(ms._Global._NMSec) + "),"
        + (ms._Global._NSec ? "    Length[" + std::to_string(ms._Global._Len.Min())
            + "-" + std::to_string(ms._Global._Len.Max()) + "]"
            + ", Tm[ " + K2s(ms._Global._Tm.Min()) + "-"
            + K2s(ms._Global._Tm.Max()) + " ]  "
            : "")
    );
}

void SeqExpl::InitTree()
    {
        _list.auto_draw(false);
        _tree.auto_draw(false);

        auto& ms = _Pr._cp._pSeqTree->MSecL();
        for ( MSecIt CurMSec = ms.begin(); CurMSec != ms.end(); CurMSec++)
            populate( AddRoot(CurMSec)) ;

        _tree.find(("Target seq")).select(true);
        populate_list_recur(*_Pr._cp._pSeqTree);

		_list.auto_draw(true);
		_tree.auto_draw(true);

    }

nana::listbox::oresolver& operator<<(nana::listbox::oresolver & ores, CMultSec::SecIt const sec_it )
{
    static const long    blen{ 50 }, slen{ 20000 };
    char val[blen];

    CSec *sec = sec_it->get();

    snprintf(val,blen,     ("%*d")  , 6,           sec->Len()       );

    ores <<  sec->Name()                                         // col 0: name  
         <<  val  ;                                              // col 1: len

    Temperature t=KtoC( sec->NonDegSet() ? sec->NonDegSet()->_Local._Tm.Ave() : sec->_Tm.Ave());
    //snprintf(val,blen, (u8"% *.*f �C"), 6, 1,   t );
    Temperature min=57.0, max=63.0;   // def
    if (sec->Len() > 50)
    {
        min = 76.0;
        max = 82.0;
    }
    else
    {
        min = 45.0;
        max = 65.0;
    }

    double fade_rate=  t<min? 0.0 : t>max? 1.0 : (t-min)/(max-min);
    nana::color tc{static_cast<nana::color_rgb>(0xFFFFFFFF)} , 
                bc = nana::color(nana::colors::blue).blend( nana::colors::red, fade_rate); 

    ores << nana::listbox::cell{ temperature_to_string(t), {bc , tc} };                             //case 2: Tm

    snprintf(val,blen,     ("%*d")  , 5,           sec->Degeneracy());


    std::string desc = sec->Description();
    if (nana::review_utf8(desc))
        sec->Description(desc);

    ores <<  val                                                     // case 3: deg    
         <<  desc   ;                                                // case 4: descr  

    if( sec->_aln_fragment && sec->_aln_fragment->aln.lenght())
    {
        snprintf(val,blen,     ("%*d")  , 6,           sec->_aln_fragment ->aln.Min()       );
        ores <<  val                       ;                                                    // case 5: beg in aln   
        snprintf(val,blen,     ("%*d")  , 6,           sec->_aln_fragment ->aln.Max()       );
        ores <<  val                       ;                                                    // case 6: end in aln    
    }
    else if (sec->_aln_fragment && sec->_aln_fragment->sq.lenght())
    {
        snprintf(val, blen, ("%*d"), 6, sec->_aln_fragment->sq.Min());
        ores << val;                                                                            // case 5: beg in aln   
        snprintf(val, blen, ("%*d"), 6, sec->_aln_fragment->sq.Max());
        ores << val;                                                                            // case 6: end in aln    
    }
    else

        ores  << " - "      << " - "        ;                                                   // no pos in aln 

    ores <<  (char *)(  sec->Sequence().substr (1, std::min( sec->Len(), slen)).c_str()    );   // sec                                                     

    return ores;
}


/*
// a "trick" to get an iterator to the CMultiSec *   !!!!!!!!!
auto it=std::find_if(parentMs->MSecL().begin(), parentMs->MSecL().end(),
                     [&ms](auto & sp_ms) {return ms == sp_ms.get(); });

if (it == parentMs->MSecL().end()) return;
*/
