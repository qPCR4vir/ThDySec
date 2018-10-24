/**
* Copyright (C) 2009-2017, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\SeqExpl.cpp
*
* @brief 
*
*/

#include "ThDy_DNAHybrid.Nana\SeqExpl.h"
#include "ThDy_DNAHybrid.Nana\main.Nana.h"

#include <algorithm>

namespace fs = std::experimental::filesystem ;

SeqExpl::SeqExpl              (ThDyNanaForm& tdForm)
        : _Pr             (tdForm), 
          CompoWidget     (tdForm, ("Seq Explorer"), ("SeqExpl.lay.txt"))
    {
        auto& sch = _list.scheme();
		sch.header_padding_bottom = sch.header_padding_top = 1;//sch.header_height = 20;
        sch.text_margin   = 2;
        sch.item_height_ex= 1;  ///< Set !=0 !!!!  def=6. item_height = text_height + item_height_ex
        //sch.item_height   = sch.text_height + sch.item_height_ex;
        sch.header_splitter_area_before = 4;
        sch.header_splitter_area_after  = 4 ; 

        //auto& tree_sch = _list.scheme();
        //tree_sch.item_height_ex = 1;  ///< Set !=0 !!!!  def=6. item_height = text_height + item_height_ex
        //tree_sch.item_height = tree_sch.text_height + tree_sch.item_height_ex;

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
        //_list.resolver(ListSeqMaker());


        AddMenuItems(_menuProgram);
        MakeResponive();

        _tree.auto_draw(true);
        _list.auto_draw(true);
    }

SeqExpl::Node SeqExpl::AddNewSeqGr  (Tree::item_proxy node)
        {    try{    
                    return appendNewNode(node, _Pr._cp.AddSeqGroup(node.value<CMultSec*>(),"New group")).expand(true);
                }
                catch ( std::exception& e)
                { 
                  (nana::msgbox ( ("Error adding new group" ) )<< e.what()).show() ;
                  return node;
                }        
        }

SeqExpl::Node SeqExpl::AddMSeqFiles (const std::string &file, bool  all_in_dir) 
    {     
    try{ 
            auto      tn    = _tree.selected();
            CMultSec* ms    = tn.value<CMultSec*>();
            CMultSec* newms = _Pr._cp.AddSeqFromFile    (ms , file, all_in_dir    );
            return Refresh(   tn);
        }
        catch ( std::exception& e)
        { 
            (nana::msgbox ( ("Error adding sequences" ) )<< e.what()).show() ;
            return _tree.selected();
         }         
    }

void SeqExpl::AddMenuItems(nana::menu& menu)
    {
        menu.append_splitter();

        menu.append(("Add a new, empty, group for sequences")  , [&](nana::menu::item_proxy& ip) {  AddNewSeqGr(_tree.selected());    });
        menu.append(("Add a group of sequences from a file..."), [&](nana::menu::item_proxy& ip) {  Click(_loadFile);                 });
        menu.append(("Add a tree of groups of sequences from a directory..."),[&](nana::menu::item_proxy& ip) {  Click(_loadDir);     });

        menu.append_splitter();

        menu.append(("Reproduce only the structure of directory..."),[&](nana::menu::item_proxy& ip)  {  Click(_scanDir);     });
        menu.append(("Reload from the original file" )  , [&](nana::menu::item_proxy& ip)   {  ReloadFile(_tree.selected());    });
        menu.append(("Reload from the original directory"), [&](nana::menu::item_proxy& ip) {  ReloadDir(_tree.selected());     });
        menu.append(("Replace from a file . . ." )  , [&](nana::menu::item_proxy& ip) 
        {
            auto tn= _tree.selected();
            if (isRoot(tn))
            {
                nana::msgbox ( ("Sorry, you can't replace group " + tn.text()) ).show() ;
                return;
            }
            CMultSec *ms = tn.value<CMultSec*>();
            CMultSec *pms = ms->_parentMS;
            assert(ms);
            assert(pms);
            fs::path pt{ms->_Path};
            pt= fs::canonical(pt).make_preferred();
            nana::filebox  fb{ *this, true };
            fb .add_filter ( SetupPage::FastaFiltre( )   )
               .init_file(pt.string())
               .title(("Replace/reload a group of sequences from a file"));
            if (!fb()) return;

            _Pr._cp._pSeqNoUsed->AddMultiSec(ms);
            _Pr._cp.AddSeqFromFile    ( pms, fb.file(), false    );
            Refresh(tn->owner());
        });

        menu.append(("Replace from directory . . ."), [&](nana::menu::item_proxy& ip) 
        {
            auto tn= _tree.selected();
            if (tn->owner()->owner().empty())
            {
                nana::msgbox ( ("Sorry, you can't replace group " + tn->text()) ) ;
                return;
            }
            CMultSec *ms = tn.value<CMultSec*>();
            CMultSec *pms = ms->_parentMS;
            assert(ms);
            assert(pms);
            std::cout<<"\n Replacing "<< ms->_Path << " from "  << pms->_Path ;

            fs::path pt{ms->_Path};
            nana::folderbox  fb{ *this, pt.parent_path() , "Replace/reload a group of sequences from a directory" };
            auto p=fb();
            if (!p) return;
            auto it=std::find_if(pms->MSecL().begin(), pms->MSecL().end(),
                                 [&ms](auto & sp_ms) {return ms == sp_ms.get(); });
            if (it == pms->MSecL().end()) return;
            _Pr._cp._pSeqNoUsed->MoveMSec(it);     // hight level MoveMSec !! (actualize globals) ?
            
            auto own = tn->owner();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            CMultSec* newms = _Pr._cp.AddSeqFromFile    ( pms, p->string(), true    );
            _tree.erase(tn);
            populate(appendNewNode  (own, newms) );
            own.expand(true);

            _list.clear();
            populate_list_recur(pms);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });

        menu.append_splitter();

        menu.append     ( ("Show Only local sequences"),[&](nana::menu::item_proxy& ip) { ShowLocals( menu.checked(ip.index())); })
            .check_style( nana::menu::checks::option)
            .checked    ( false );

        menu.append     ( ("Show filtered sequences"  ),[&](nana::menu::item_proxy& ip) { ShowFiltered( menu.checked(ip.index())); })
            .check_style( nana::menu::checks::highlight )
            .checked    ( true );
 
        menu.append_splitter();
        menu.append(("Cut selected sequences from list"          ),[&](nana::menu::item_proxy& ip)  {  Click(_cutSec);  });
        menu.append(("Cut selected groups of sequences from tree"),[&](nana::menu::item_proxy& ip)  {  Click(_cut   );  });
        menu.append(("Paste the sequences"                       ),[&](nana::menu::item_proxy& ip)  {  Click(_paste );  });

        menu.append_splitter();
        menu.append(("Del selected sequences from list"),[&](nana::menu::item_proxy& ip)            {  Click(_delSec);  });
        menu.append(("Del selected groups of sequences from tree"),[&](nana::menu::item_proxy& ip)  {  Click(_del   );  });
        menu.append(("Rename the selected group of sequences"),[&](nana::menu::item_proxy& ip) 
        {
            
            RenameFrom rnm(_tree, _tree.selected().text());
            nana::API::modal_window( rnm );
            _tree.selected().text(rnm.Name());
            _tree.selected().value<CMultSec*>()->_name = rnm.Name() ;

        }).enabled(true);

        menu.append_splitter();
        auto  indxFASTA = menu.append(("Export FASTA . . ."          ),[&](nana::menu::item_proxy& ip)            {  /*_tree.selected().value<CMultSec*>()->ExportFASTA();*/  }).index();
        auto& menuFASTA = *menu.create_sub_menu(indxFASTA);
        menuFASTA.append(("Only current sequences"     ),[&](nana::menu::item_proxy& ip)            {  ;  });
        menuFASTA.append(("Selected sequences in group"),[&](nana::menu::item_proxy& ip)            { _Pr.ExportFASTA(_tree.selected().value<CMultSec*>(), true );  });
        menuFASTA.append(("All sequences in group"     ),[&](nana::menu::item_proxy& ip)            { _Pr.ExportFASTA(_tree.selected().value<CMultSec*>(), false);  });
        menuFASTA.append(("All selected sequences"     ),[&](nana::menu::item_proxy& ip)            { _Pr._cp._pSeqTargets->Export_as("export.fasta", true )  ;  });
        menuFASTA.append(("All sequences"              ),[&](nana::menu::item_proxy& ip)            { _Pr._cp._pSeqTargets->Export_as("export.fasta", false)  ;  });

    }

void SeqExpl::MakeResponive()
    {
        // the "selected" feature in the GUI have no efect in the data, it is a pure GUI feature,
        // but the check status go to the data with the name Selected() !!!

        _tree.events().selected ( [&]( const nana::arg_treebox &tbox_arg_info ) 
                                 { if (tbox_arg_info.operated) RefreshList(tbox_arg_info.item); });
        _tree.events().checked  ( [&]( const nana::arg_treebox &tbox_arg_info )
        {                                              
            tbox_arg_info.item.value<CMultSec*>()->Selected(tbox_arg_info.operated);
            if (tbox_arg_info.item== _tree.selected())  
                RefreshList(tbox_arg_info.item);                //  ??????? Only RefreschList
        });

        _list.events().checked  ( [&](  const nana::arg_listbox &lbox_arg_info )
        {                                               
            lbox_arg_info.item.value<CSec*>()->Selected(lbox_arg_info.item.checked());
        });
 
        _loadFile   .events().click([this]()
                        {
                            auto      tn    = _tree.selected();
                            CMultSec* ms    = tn.value<CMultSec*>();
                            fs::path pt{ms->_Path};
                            pt= fs::canonical(pt).make_preferred();

                            nana::filebox  fb{ *this, true };
                            fb .add_filter ( SetupPage::FastaFiltre( )                   )
                                    .init_file(pt.string())
                                    .title      ( ("File load: Add a group of sequences from a file") );

                            if (fb()) 
                               AddMSeqFiles(fb.file(), false);
                        });
        //_loadFileTT.set(_loadFile,("File load: Add a group of sequences from a file"));

        _re_loadFile.events().click([this]()  {  ReloadFile(_tree.selected());    });
        _re_loadFileTT.set(_re_loadFile,("File reload: Reload a group of sequences from a file, \nposible using new filtres."));

        _loadDir    .tooltip(("Directory load: Add a tree of groups of sequences from a directory."))
                    .events().click([this]()
                        {
                            nana::filebox  fb{ *this, true };
                            fb .add_filter ( SetupPage::FastaFiltre( )                   )
                               .title(("Directory load: Add a tree of groups of sequences from a directory"));
                            if (fb()) 
                                AddMSeqFiles(fb.file(), true);
                        });
        _re_loadDir .tooltip(("Directory reload: Reload a tree of groups of sequences from a directory,\nposible using new filtres."))
                    . events().click([this]()  {  ReloadDir (_tree.selected());    });
        _scanDir    .tooltip(("Directory scan: Reproduce the structure of directory..."))
                    .events().click([this]()
                        {
                            nana::filebox  fb{ *this, true };
                            fb .add_filter ( SetupPage::FastaFiltre( )                   )
                               .title(("Directory scan: Reproduce the structure of directory..."));
                            if (!fb()) return;

                            auto      tn    = _tree.selected();
                            CMultSec* ms    = tn.value<CMultSec*>();
                            CMultSec* newms = _Pr._cp.CopyStructFromDir    ( ms, fb.file()    );
                            _tree.auto_draw(false);
                            populate(  appendNewNode  (tn, newms) );
                            tn.expand(true);
                            _tree.auto_draw(true);
                        });
        _paste      .tooltip(("Paste sequences"))
                    .events().click([this]()
        {
            auto       tn = _tree.selected();
            CMultSec *pms = tn.value<CMultSec*>();

            for (auto ms : _dragMSec)
                pms->AddMultiSec(ms);
            for (auto s : _dragSec)
                pms->AddSec(s);

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
        _cut        .tooltip(("Cut a group of sequences"))
                    .events().click([this]()
        {
            auto tn= _tree.selected(); 
            if (tn->owner()->owner().empty())    //   ???  if( tn->level() < 2 );
            {
                (nana::msgbox ( _tree , ("Cut a group of sequences " + tn->text()) )
                          << ("Sorry, you can�t cut the group: ") + tn->text() )
                          .icon(nana::msgbox::icon_error )
                          .show() ;
                return;
            }
            CMultSec *ms = tn.value<CMultSec*>();
            CMultSec *pms = ms->_parentMS;  
            assert(ms);
            assert(pms);

            _Pr._cp._pSeqNoUsed->AddMultiSec(ms);
            _dragMSec.push_back(ms);
            //ms->MoveBefore(_Pr._cp._pSeqNoUsed->goFirstMSec() );  /// \todo: higth level MoveMSec !! (actualize globals)
            auto own = tn->owner();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            _tree.erase(tn);
            populate(appendNewNode (_tree.find(("Don t use") ), ms ));
            own.select(true).expand(true);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });
        _del        .tooltip(("Delete a group of sequences "))
                    .events().click([this]()
        {
            auto tn= _tree.selected();
            if (tn->owner()->owner().empty())
            {
                (nana::msgbox ( _tree , ("Deleting a group of sequences " + tn->text()) )
                          << ("Sorry, you can�t delete the group: ") + tn->text() )
                          .icon(nana::msgbox::icon_error )
                          .show() ;
                return;
            }
            CMultSec *ms = tn.value<CMultSec*>();
            CMultSec *pms = ms->_parentMS;     
            assert(ms);
            assert(pms);

            _Pr._cp._pSeqNoUsed->AddMultiSec(ms); //ms->MoveBefore(_Pr._cp._pSeqNoUsed->goFirstMSec() );  /// \todo: higth level MoveMSec !! (actualize globals)
            auto own = tn->owner();

            _tree.auto_draw(false);
            _list.auto_draw(false);

            _tree.erase(tn);
            populate(appendNewNode (_tree.find(("Don t use") ), ms ));

            own.select(true).expand(true);

            _tree.auto_draw(true);
            _list.auto_draw(true);
        });

        _cutSec     .tooltip(("Cut selected sequences from list"))
                    .events().click([this]()
        {
            auto sel =    _list.selected() ; 
            for (auto i : sel)
            {
                auto s=_list.at(i ).value<CSec*>();
                _Pr._cp._pSeqNoUsed->AddSec( s );
                _dragSec.push_back(s);
            }
            RefreshList();
        });

        _delSec     .tooltip(("Delete selected sequences from list"))
                    .events().click([this]()
        {
            auto sel =    _list.selected() ; 
            for (auto i : sel)
            {
                auto s=_list.at(i ).value<CSec*>();
                _Pr._cp._pSeqNoUsed->AddSec( s );
            }
            RefreshList();
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

SeqExpl::Node SeqExpl::Replace      (Tree::item_proxy tn, CMultSec *ms, const std::string& Path, bool all_in_dir)
{        
try{ 
        Tree::item_proxy   own = tn->owner();
        CMultSec          *pms = ms->_parentMS;  
        assert(ms);
        assert(pms);


        CMultSec* newms = _Pr._cp.AddSeqFromFile    ( pms, Path, all_in_dir    );

        _Pr._cp._pSeqNoUsed->AddMultiSec(ms); 
        _tree.erase(tn);
        populate(_tree.find( _Pr._cp._pSeqNoUsed->_name));
        return appendNewNode(own, newms).expand(true).select(true) ;
    }
    catch ( std::exception& e)
    { 
        (nana::msgbox ( ("Error replacing sequences: " ) ).icon(nana::msgbox::icon_error)
            << "into group:    "  << tn.key()                                 
            << "\n from " << (all_in_dir?"directory: " : "file: ") << Path     <<"\n"<< e.what()
        ).show() ;
     }        
    catch(...)
    {
            (nana::msgbox(("An uncaptured exception during replacing sequences: "))
                .icon(nana::msgbox::icon_error) 
            << "into "<< tn.key()                                
            << "from "<< (all_in_dir?"directory ":"file ") << Path    
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
    auto probNode = _tree.find(_Pr._mPCR._probesMS->_name);
    Refresh(probNode).expand(true).select(true);
    if (show_) 
        _Pr.ShowExpl();
}


void SeqExpl::SetDefLayout()
{
    _DefLayout =
        "vertical                                                 \n\t"
        "          <weight=23 <toolbar weight=680 margin=2 ><>>    \n\t"
        "          <      <Tree  > |75% <List >   >                \n\t"

        "          <weight=21 <weight=5> <weight=80 tflab> <min=640 max=1000  TargetsOptions   >  <weight=5>>  \n\t"
        "          <weight=20 <statusbar margin=1  min=700 > <Firma max=180> <weight=3 > >   \n\t"
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
            << "      Dir:" << _loadDir << _re_loadDir << _scanDir << _cut << _del
            << "      Seq:" << _show_locals_s << _show_filt_s << _cutSec << _delSec
            ;
        _place["Tree"  ] << _tree;
        _place["List"  ] << _list;
        _place["tflab"] << "Targets filters: ";
        _place["TargetsOptions"]    << numUpDwMaxTgId
                                    << numUpDw_TgBeg << numUpDw_TgEnd
                                    << numUpDw_SLenMin << numUpDw_SLenMax;

        _place["statusbar"] << _statusbar;
        _place.field("Firma") << e_mail_firma;

    }
std::string K2s(Temperature t) { return temperature_to_string(KtoC(t)); }
void SeqExpl::RefreshStatusInfo(CMultSec *ms)
{
	//std::string local = 
    _statusbar.caption("  <bold>Local</> (sq: " + std::to_string(ms->_Local._NSec)
        + ", gr: " + std::to_string(ms->_Local._NMSec) + "),"
        + (ms->_Local._NSec ? "    Length[" + std::to_string(ms->_Local._Len.Min())
            + "-" + std::to_string(ms->_Local._Len.Max()) + "]"
            + ", Tm[ " + K2s(ms->_Local._Tm.Min()) + "-"
            + K2s(ms->_Local._Tm.Max()) + "]  "
            : "")
        + ".        <bold>Global</> (sq: " + std::to_string(ms->_Global._NSec)
        + ", gr: " + std::to_string(ms->_Global._NMSec) + "),"
        + (ms->_Global._NSec ? "    Length[" + std::to_string(ms->_Global._Len.Min())
            + "-" + std::to_string(ms->_Global._Len.Max()) + "]"
            + ", Tm[ " + K2s(ms->_Global._Tm.Min()) + "-"
            + K2s(ms->_Global._Tm.Max()) + " ]  "
            : "")
    );
}


void SeqExpl::InitTree()
    {
        _list.auto_draw(false);
        _tree.auto_draw(false);

        CMultSec *ms=_Pr._cp._pSeqTree.get();
        for ( auto& CurMSec :  ms->MSecL() )
            populate( AddRoot( CurMSec.get() )) ;

        _tree.find(("Target seq")).select(true);
        populate_list_recur(_Pr._cp._pSeqTree.get());

		_list.auto_draw(true);
		_tree.auto_draw(true);

    }

List::oresolver& operator<<(List::oresolver & ores, CSec * const sec )
{
    static const long    blen{ 50 }, slen{ 1000 };
    char val[blen];

    snprintf(val,blen,     ("%*d")  , 6,           sec->Len()       );

    ores <<  sec->Name()                                         // col 0: name  
         <<  val  ;                                              // col 1: len

    Temperature t=KtoC( sec->NonDegSet() ? sec->NonDegSet()->_Local._Tm.Ave() : sec->_Tm.Ave());
    //snprintf(val,blen, (u8"% *.*f �C"), 6, 1,   t );
    Temperature min=57.0, max=63.0;
    double fade_rate=  t<min? 0.0 : t>max? 1.0 : (t-min)/(max-min);
    nana::color tc{static_cast<nana::color_rgb>(0xFFFFFFFF)} , 
                bc = nana::color(nana::colors::red).blend( nana::colors::blue, fade_rate); 

    ores << List::cell{ temperature_to_string(t), {bc , tc} };                             //case 2: Tm 

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