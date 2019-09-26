/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2009-2019
*
* @file  ThDySec/include/ThDy_DNAHybrid.Nana/SeqExpl.h
*
* @brief GUI to explore sequences
*
*/

#ifndef SeqExpl_H
#define SeqExpl_H

#include <nana/gui/widgets/treebox.hpp>
#include <nana/gui/widgets/listbox.hpp>
#include <nana/gui/tooltip.hpp>
#include <nana/gui/widgets/toolbar.hpp>

#include <nanaBind.hpp>
#include <EditableForm.hpp>

#include "ThDy_programs/init_ThDy_prog_param.h"

class ThDyNanaForm ;
 
inline std::string temperature_to_string(Temperature t) // use KtoC to convert from Kelvin to Centigr
{
	int w = 10;  // -250,0 �C
	std::string s(w, 0);
	const static std::string format = std::string("% *.*f ") + RTunits::grC ;  // \u00B0C"
	std::snprintf(&s[0], w + 1, (format.c_str()), w - 4, 1, t);
	return s;
}

class SeqExpl : public CompoWidget
{
    using Tree   = nana::treebox;
    using Node   = Tree::item_proxy;
    using MSecIt = CMultSec::MSecIt;   // the "value type" for Tree Node
    using List   = nana::listbox;
    using SecIt  = CMultSec::SecIt;    // the "value type" for List Item

    ThDyNanaForm       &_Pr;
    Tree                _tree{ *this };
    List                _list{ *this };
	nana::toolbar       _tbar{ *this };
    bool				_showAllseq{true}, _showFiltered{true};
    std::vector<SecIt>  _dragSec;
    std::vector<MSecIt> _dragMSec;

    nana::button    _loadFile     {*this, "Load"  },       //nana::toolbar  _tbar { *this };
                    _re_loadFile  {*this, "reLoad"},
                    _loadDir      {*this, "Load"  },
                    _re_loadDir   {*this, "reLoad"},
                    _scanDir      {*this, "Scan"  },
                    _re_scanDir   {*this, "reScan"},
                    _cut          {*this, "Cut"   },
                    _paste        {*this, "Paste" },
                    _del          {*this, "Del"   },
                    _cutSec       {*this, "Cut"   },
                    _delSec       {*this, "Del"   },
                    _show_locals_s{*this, "local" },
                    _show_filt_s  {*this, "filter"}
                    ; 

	ParamGUIBind::BindGroup    _commPP;
	nana::NumUnitUpDown numUpDwMaxTgId  {*this, "Max. ident.:",99,50, 100 ,  "%" },
                        numUpDw_TgBeg   {*this, "Beg.:"       , 0, 0, 100000,"nt"},    /// rev !!
                        numUpDw_TgEnd   {*this, "End.:"       , 0, 0, 100000,"nt"},    /// rev !!
                        numUpDw_SLenMin {*this, "Min.Len.:"   , 0, 0, 100000,"nt"},
                        numUpDw_SLenMax {*this, "Max.Len.:"   , 0, 0, 100000,"nt"};

    nana::tooltip   _loadFileTT {_loadFile,"File load: Add a group of sequences from a file"},
                    _re_loadFileTT {_re_loadFile,"File reload: Reload a group of sequences from a file, \npossibly using new filters."};

    nana::label     _statusbar    { *this },
                    _firma ;

    void SetDefLayout() override;
    void AsignWidgetToFields() override;
    void MakeResponive();

    Node Refresh(Node node)
    {
            _tree.auto_draw(false);
            try{
                populate(node);
                node.expand(true);
                RefreshList(node);
            }
            catch ( std::exception& e )      //
            {
                (nana::msgbox(*this, "Error during sequence lists displaying !\n\t", nana::msgbox::button_t::ok)
                        .icon(nana::msgbox::icon_information )
                        << e.what()
                ).show (  ) ;
            }
            _tree.auto_draw(true);

            return node;
    }
    void RefreshList(                ) { RefreshList(_tree.selected());      } ///< very High Level
    void RefreshList(const Node &node) { RefreshList( node.value<MSecIt>()); } ///< High Level
    void RefreshList(MSecIt ms)  ///< medium Level
    {
            _list.auto_draw(false);

            _list.clear();
            populate_list_recur(**ms);

            _list.auto_draw(true);
			RefreshStatusInfo(**ms);
    }

    void populate_list_recur(const Node& node) ///< High Level
    {
        populate_list_recur(**node.value<MSecIt>());
    }

    void populate_list_recur(CMultSec &ms)  ///< Low Level
    {
        populate_list(ms);
        if ( _showAllseq )
            for ( auto& CurMSec : ms.MSecL() )
                populate_list_recur(*CurMSec);
    }

    void populate_list(CMultSec &ms)
    {
        for (SecIt CurSec = ms.SecL().begin(); CurSec != ms.SecL().end(); CurSec++)
		  if ( _showFiltered || ! (**CurSec).Filtered() )
              AddSecToList(CurSec);
    }

    List::item_proxy AddSecToList     (SecIt s)
    {
        return _list.at(0).append(s)
                          .value(s)
                          .check  ( (**s).Selected() )
                          .fgcolor( static_cast<nana::color_rgb>(
                                              ((**s).Filtered() ? 0xFF00FF   ///\todo: use coding
                                                             :   0x0   )  ));//nana::color::gray_border );
    }

    Node AddRoot(MSecIt const ms)
    {
        std::string name = (**ms)._name;
        return _tree.insert(name, name).value(ms).check((**ms).Selected());
    }
    bool isRoot(const Node &node)
    {
        return node.level() == 1;
    }
 static Node appendNewNode(Node node, MSecIt const ms) /// Add a new node to the child of node.
    {
        std::string name = (**ms)._name;
        return node.append(name, name).value(ms).check((**ms).Selected());
    }

    Node populate     ( Node node)  ///< create & add to the child of node a new node nuevo for each MSec in ms.
    {
        while(node.size()) 
            _tree.erase(node.child());

        auto &ms = (*node.value<MSecIt>())->MSecL(); //  msec(node);
		for (MSecIt CurMSec = ms.begin();CurMSec != ms.end(); CurMSec++)
			populate( appendNewNode(node, CurMSec)) ;
        return node;
    }

    Node AddNewSeqGr  (Node node) ;

    /// Add seq from file to the selected tree Node
    Node AddMSeqFiles (const std::filesystem::path &file, bool  all_in_dir) ;

    /// The 'old' node tn will be eliminated, and in his previous parent a new node will be created and returned,
    /// The sequences in ms will be moved to DontUse and in his parent a new CMultiSec will be loaded from file.
    Node Replace    (Node tn, MSecIt ms, const std::filesystem::path& Path, bool all_in_dir, bool scan_only=false);

    Node ReloadFile   (Tree::item_proxy tn)
    {
        MSecIt ms = tn.value<MSecIt>();
        if ((*ms)->_orig_file.empty())
            return tn;
        else
            return Refresh(Replace(tn, ms, (*ms)->_orig_file, false));
    }

    Node ReloadDir  (Node tn, bool scan_only=false)
    {
        MSecIt ms = tn.value<MSecIt>();
        if ((*ms)->_orig_file.empty())
        {
            for (Node & ntn : tn)
                ReloadDir(ntn, scan_only);
            return tn;
        }
        else return Replace(tn, ms, (*ms)->_orig_file, true, scan_only).select(true);//true
    }

    void ShowLocals(bool showLocals)
    {        
        if(showLocals != _showAllseq) return ;
        else _showAllseq = ! showLocals;

        _list.auto_draw(false);
        _list.clear();
            populate_list_recur(_tree.selected());
        _list.auto_draw(true);
    }
    void ShowFiltered(bool showFiltered)
    {        
        if(showFiltered == _showFiltered) return ;
        else _showFiltered = showFiltered;

        _list.auto_draw(false);
        _list.clear();
            populate_list_recur(_tree.selected());
        _list.auto_draw(true);
    }
	void RefreshStatusInfo(const CMultSec &ms);
    void RefreshStatusInfo() {RefreshStatusInfo(**_tree.selected().value<MSecIt>() );};

public:
    SeqExpl(ThDyNanaForm& tdForm);
    void ShowFindedProbes_in_mPCR(bool show_=true);
    void RefreshProbes_mPCR(bool show_=true);
    void AddMenuItems(nana::menu& menu);
    void InitTree();
};

nana::listbox::oresolver& operator<<(nana::listbox::oresolver & ores, CMultSec::SecIt const sec );

class RenameFrom : public nana::form, public EditableForm
{
    std::string     _name;
    bool            _renamed    {false};
    nana::textbox   edit        {*this};
    nana:: button   OK          {*this, ("rename")}, 
                    Cancel      {*this, ("abort" )};
    ParamGUIBind::BindGroup       _bind     ;


  public:
    RenameFrom(nana::window owner, std::string name) : 
            _name(name),  
             nana::form  (nana::rectangle( nana::point(150,500), nana::size(300,50) )),
             EditableForm(owner, *this, ("Rename") )     
        {
			edit.multi_lines(false);
			edit.caption(_name);
            InitMyLayout();
            OK.events().click([this]()
            {
                std::string  name=edit.caption(); 
                _renamed = (_name!=name);
                if (_renamed)
                    _name=name;
                this->close(); 
            });
            Cancel.events().click([this](){_renamed = false; this->close(); });

        }
    std::string Name(){return _name;}
    void SetDefLayout   () override
    {
		_DefLayout = R"( vertical gap=2 
			                 <weight = 24 Edit > 
                             <weight = 24     < free_left >     <Buttons gap = 10>      >
                       )";  
	}
    void AsignWidgetToFields() override
    {
	    _place["Edit"    ] << edit  ;
	    _place["Buttons" ] << OK <<   Cancel  ;
    }                                        

};

#endif