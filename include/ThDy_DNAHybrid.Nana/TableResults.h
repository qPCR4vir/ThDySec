/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
*
* @file  ThDySec\include\ThDy_DNAHybrid.Nana\TableResults.h
*
* @brief 
*
*/

#ifndef TableResults_H
#define TableResults_H

#include "thdy_programs\init_thdy_prog_param.h"
#include "matrix.h" 
#include <../../nana.ext/include/EditableForm.hpp>
#include <nana/gui/widgets/listbox.hpp>
//#include <nana/gui/tooltip.hpp>

//#include <nana/gui/widgets/checkbox.hpp>
//#include <nana/gui/widgets/toolbar.hpp>



using List = nana::listbox;


class TableHybRes  : public nana::form, public EditableForm
{   
    using Table = CTable<TmGPos> ;
    using index = Table::index;
    struct value
    {
        Table  *table;
        int      n_len{ 6 }, n_dec{ 1 }  ;

        virtual ~value(        ){}
        value (Table &t) :table {&t}{}

        virtual float val  (index row,  index col)const =0 ;
                float operator()(index row,  index col){return val(row,col);}
        
        virtual std::string str(index row,  index col) const
        {
            std::string s (n_len, 0);

            auto l=snprintf( &s[0], n_len+1 , ("% *.*f"), n_len, n_dec, val(row,col) );
            //s.resize(l);
            return s;
        }

        virtual bool         return_bg(){return false;}
        virtual nana::color bg_color(index row,  index col, List &lst)
		{
			return lst.bgcolor(); 
			//return List::scheme_type().header_bgcolor.get_color();
		}
    }   ;
    struct Tm : value
    {
        float val(index row,  index col) const override
        {
            return table->at(row,col )._Tm;
        }
        Tm(Table &t) :value {t}{};
        bool return_bg() override {return true;}
        nana::color bg_color(index row,  index col, List &lst) override
        {
            Temperature t=val(row,col);

            Temperature min=20.0, max=63.0;
            double fade_rate=  t<min? 0.0 : t>max? 1.0 : (t-min)/(max-min);
            nana::color bgc(nana::colors::red);
            return bgc.blend(nana::colors::blue, fade_rate) ;
        }
    };
    struct G : value
    {
        float val(index row,  index col) const override
        {
            return table->at(row,col )._G;
        }

        G(Table &t) :value {t}{};

		bool return_bg() override { return true; }

		nana::color bg_color(index row, index col, List &lst) override
		{
			Energy e = val(row, col);

			Energy min = -15.0, max = +15.0;
			double fade_rate = e<min ? 0.0 : e>max ? 1.0 : (e - min) / (max - min);
			nana::color bgc(nana::colors::red);
			return bgc.blend(nana::colors::blue, fade_rate);
		}

	
	};
    struct Pos : value
    {
        float val(index row,  index col)const override
        {
            return static_cast<float>( table->at(row,col)._Pos );
        }
        Pos(Table &t) :value {t} {n_dec=0;};
   };

    std::shared_ptr<Table> _table;
    List                   _list { *this };

    nana::button           _bTm  {*this,("Tm" )},       //nana::toolbar     _tbar { *this };
                           _bG   {*this,("G"  )},   
                           _bPos {*this,("Pos")},
                           _mix  {*this, ("Consolidate")}; 

    Tm                     _Tm;
    G                      _G;
    Pos                    _Pos;
    value                  *val { &_Tm} ;
    std::size_t            mTm, mG, mP;
 
    void SetValType(value &val_)
    { 
        _list.auto_draw(false);
        bool freeze{true};
        freeze=_list.freeze_sort(freeze);

        val = &val_;
        Repopulate();
        
        _list.auto_draw(true);
        //_list.unsort();
        _list.freeze_sort(freeze);
    }
    void Repopulate()
    {
        for (auto &i : _list.at(0))
            i.resolve_from(i.value<Index>());
    }

    bool comp(index col, nana::any* row1_, nana::any*row2_, bool reverse)
    {
                float  v1{ val->val( nana::any_cast<Index>( row1_)->row , col-1) },
                       v2{ val->val( nana::any_cast<Index>( row2_)->row , col-1) };
                return reverse?  v2<v1 : v1<v2 ;
    }
    void SetDefLayout   () override
    {
        _DefLayout= 
	        "vertical                  	\n\t"
	        "	  <weight=25 <toolbar weight=200 ><>>       	\n\t"
	        "	   <_list  >       	\n\t"
	        "	 	\n\t"
 
                    ;
    }
    void AsignWidgetToFields() override
    {
 	    _place.field("toolbar"       ) <<_bTm << _bG << _bPos ;
 	    _place.field("_list"         ) <<_list;
		_place.field("Firma") << " ArielVina.Rodriguez@fli.bund.de";

     }
 public:
     TableHybRes    (std::shared_ptr<CTable<TmGPos>> table)  : 
                            _table(table), 
                            _Tm{*table.get()}, _G{*table.get()}, _Pos{*table.get()},  
                            nana::form (nana::rectangle( nana::point(50,5), nana::size(1000,650) )),
                            EditableForm    (nullptr, *this,  table->TitTable(), "TableTm.lay.txt") 
   {
        caption( std::string(("Table Tm: ")) +  _Titel);

		auto& sch = _list.scheme();
		sch.header_height = 20;
		sch.text_margin   = 2;
		sch.item_height_ex= 1;  ///< Set !=0 !!!!  def=6. item_height = text_height + item_height_ex
		sch.item_height   = sch.text_height + sch.item_height_ex;
		sch.header_splitter_area_before = 4;
		sch.header_splitter_area_after = 4 ;

        InitMyLayout();
        SelectClickableWidget( _list);
        SelectClickableWidget( *this);

        _list.auto_draw(false);
                
        _list.append_header(("Seq")  , 120);
        for (index col = 1; col <= table->totalCol(); ++col)
        {    
            _list.append_header(   table->TitColumn(col-1)   , 100);
            _list.set_sort_compare(col,[col,this](const std::string&, nana::any* row1_, const std::string&, nana::any*row2_, bool reverse)
            {
                 return comp(col,row1_,row2_,reverse);
            });
        }

        for (index row = 0; row < table->totalRow(); ++row)
            _list.at(0).append( Index{this,row}, true );

        _list.auto_draw(true);

        //MakeResponive();
        _bTm .events().click([this]()
                        {
                            SetFormat(1);
                            SetValType(_Tm);
                            caption( std::string(("Table Tm: ")) +  _Titel);
                            _bTm .pushed(true);
                            _bG  .pushed(false);
                            _bPos.pushed(false);
                            _menuProgram.checked(mTm, true);
                        });
        _bG  .events().click([this]()
                        {
                            SetFormat(1);
                            SetValType(_G);
                            caption( std::string(("Table G: ")) +  _Titel);
                            _bTm .pushed(false);
                            _bG  .pushed(true);
                            _bPos.pushed(false);
                            _menuProgram.checked(mG, true);
                        });
        _bPos.events().click([this]()
                        {
                            SetFormat(0);
                            SetValType(_Pos);
                            caption( std::string(("Table Pos: ")) +  _Titel);
                            _bTm .pushed(false);
                            _bG  .pushed(false);
                            _bPos.pushed(true );
                            _menuProgram.checked(mP, true);
                        });

        _bTm .enable_pushed(true).pushed(true);
        _bG  .enable_pushed(true).pushed(false);
        _bPos.enable_pushed(true).pushed(false);


        _menuProgram.append_splitter();
        
        mTm=_menuProgram.append     ( ("Show Tm")    , [&](nana::menu::item_proxy& ip)  { Click( _bTm); })
                        .check_style( nana::menu::checks::option)
                        .index();

        mG=_menuProgram.append      ( ("Show delta G"), [&](nana::menu::item_proxy& ip) { Click( _bG);  })
                        .check_style( nana::menu::checks::option)
                        .index();

        mP=_menuProgram.append      ( ("Show Pos")    , [&](nana::menu::item_proxy& ip) { Click( _bPos);})
                        .check_style( nana::menu::checks::option)
                        .index();
    }
    void SetFormat(int dec=1 , int len=6)  // ??
    {  
        _Tm.n_len=_G.n_len=len; _Tm.n_dec=_G.n_dec=dec;
    }

    friend struct Index;

    struct Index
    {
        TableHybRes* table;
        index       row;

        friend List::oresolver& operator<<(List::oresolver& ores, const TableHybRes::Index& i)
        {
            auto &t = *i.table->_table.get();
            auto &v = *i.table->val;
			auto &l = i.table->_list ;
            ores<< t.TitRow(i.row)   ;
                
            if  (v.return_bg() )
                for (int col=0; col< t.totalCol() ; ++col)
                    ores<< List::cell{ v.str     (i.row, col),
                                       v.bg_color(i.row, col, l),
                                       nana::colors::white};
            else 
                for (int col=0; col< t.totalCol() ; ++col)
                    ores<< v.str(i.row, col)  ;

            return ores;
        }
    };
    friend List::oresolver& operator<<(List::oresolver& ores, const TableHybRes::Index& i);

};

#endif