/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2019
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\main.Nana.cpp
*
* @brief Entry point (main) for ThDyHybr with Nana GUI
*
*
*/

#include "ThDy_DNAHybrid.Nana\main.Nana.h"


//if you want to keep the Windows subsystem you can just hint at what your entry point is, 
//because you haven't defined ___tmainCRTStartup. You can do this by adding the following 
//to Properties -> Linker -> Command line:
//
//    /ENTRY:"mainCRTStartup"
//
//This way you get rid of the console window.

 class About : public nana::form, public EditableForm
 {
	 nana::group c_r {*this, "Author:"};
     nana::label  copy_r{ c_r, " Copyright (C) 2009-2020, Ariel Vina-Rodriguez                ( arielvina@yahoo.es, <bold> qPCR4vir</> at <bold blue url=\"https://qpcr4vir.github.io/\"> GitHub</>)" };
	 nana::label  comments{ c_r,R"(This work is mentioned in <bold blue url="https://epub.ub.uni-greifswald.de/frontdoor/deliver/index/docId/2175/file/VinaRodriguez.2018.Dissertation.pdf"> my PhD thesis </>.

 Program distributed under the <bold blue url="https://www.gnu.org/licenses/"> GNU General Public License</>)" };

     nana::group build {*this, "Programm:"};
     nana::label  compiled{ build, R"(Compiled on:   )"  __DATE__   R"( / )"  __TIME__   R"(    Version: v0.02.12    )" };
     nana::label  downloads{ build, R"(Downloads and source code: <bold blue url="https://github.com/qPCR4vir/ThDySec"> github.com/qPCR4vir/ThDySec</>
                                             <bold blue url=" https://github.com/qPCR4vir/ThDySec/wiki "> Wiki:</>
   ________________________________________________________________________________)" };

     nana::label  GUI_lib{ build, R"(Powered by <bold blue url="http://nanapro.org/en-us/ "> Nana C++ GUI library:</>
                                        <bold blue url="https://github.com/qPCR4vir/nana-docs/wiki"> Wiki:</>
                      <bold blue url="https://github.com/qPCR4vir/nana/commit/dbf6a7eebac2e4a6412979adaa2f8a2bd7e9b58f">  Nana Version: 1.7.2 develop</>)" };


	 
 nana::button bclose { *this, "Close" };
	 

 public:
	 About() : nana::form   ( nana::rectangle(nana::point(50, 5), nana::size(510, 460)) ),
	           EditableForm ( nullptr, *this, "About ThDy Hybrid" , "about.lay.txt") 
	 {
         copy_r.format(true);
         comments.format(true);
         compiled.format(true);
         GUI_lib.format(true);
         downloads.format(true);
         // .caption(" Copyright (C) 2009-2018, Ariel Vina-Rodriguez ( arielvina@yahoo.es, <bold> qPCR4vir</> at <bold blue url=\"https://qpcr4vir.github.io/\"> GitHub</>)");
		 InitMyLayout();
		 SelectClickableWidget(c_r);
         SelectClickableWidget(build);
         SelectClickableWidget(copy_r);

         bclose.events().click([&](){this->close();});
	 }

 void SetDefLayout() override
	 {
		 _DefLayout = R"(
< vertical     gap=5 margin=[3,20,3,20]
      <CR     min=125 max=150 >
      <all      min=170 max=200 >
     < min=1  max=180>
      <close  min=25  max=80>
     < min=10  max=20>
 >    )";
		 c_r.div     (R"(vert< CR  vert vfit=290 height=60  gap=5 margin=[5,20,3,20] > )");
         build.div   (R"(vert< all vert                     gap=5 margin=[3,20,3,20] > )");
	 }
 void AsignWidgetToFields()
 {
     _place["CR"]    << c_r  ;
     _place["all"]   << build   ;
     _place["close"] << bclose ;

     c_r   ["CR" ] << copy_r   << comments;
     build ["all"] << compiled << downloads << GUI_lib ;
 }
 };


int main(int argc, char *argv[]) 
{
  try	
  {
    using namespace ParamGUIBind;
    IParBind::SetDef(PriorizeDefault::Parametr );

    ThDyNanaForm tdForm(  argc,  argv);
    tdForm.show();
	nana::exec();
	return 0;
  }
    catch (std::exception& e)
        {
            std::cerr<< std::endl<< e.what();
			//nana::msgbox(e.what()).show();    // --> when we are here the GUI is already "close"
			//nana::exec();
			char c; std::cin >> c ;
            //throw ;
	} 
    catch (...)
        {
            std::cerr<< std::endl<< "exception !!";
			char c; std::cin >> c;
			//throw ;
        }
    //std::cout.rdbuf(coutbuf); //reset to standard output again
} 


    ThDyNanaForm::ThDyNanaForm (int argc, char *argv[]) try
                  :nana::form (nana::rectangle( nana::point(50,5), nana::size(1000,650) )),
                   EditableForm    (nullptr, *this,  "ThDy DNA Hybrid" " (" __DATE__ " / " __TIME__")" , "ThDy.lay.txt"),
				   e_mail_firma (" arielvina@yahoo.es")
{
    //nana::API::zoom_window(*this, true);
    //nana::pixel_rgb_t bk;
    //bk.u.color = background ();
    //bk.u.element.blue =0; 
    //bgcolor (nana::color_rgb( 0xEEEEEE));  ///\todo: use codigo
    //foreground(1);
    //this->scheme().background(col)
    //this->b
    add_page( setup_    );// setup_.ReCollocate(); // 0 
    add_page( mExpl_    ); //mExpl_.ReCollocate();// 1
    add_page( findSond_ );// findSond_.ReCollocate();// 2
    add_page( mPCR_     ); //mPCR_.ReCollocate();// 3
    add_page( uArr_     );// uArr_.ReCollocate();// 4
    add_page( tmCalc_   );// tmCalc_.ReCollocate();// 5

    tabbar_.activated (1);

    setup_._proj.FileNameOnly( ProjectFile()  );
    try{ 
			if ( argc > 1 )
				setup_._proj.FileNameOpen( argv[1])    ;
			else
				load() ;						
		}
    catch ( std::exception& e )      // Por ejemplo cuando no existe Def Project: 1ra vez que se usa el prog.
	{   
        (nana::msgbox(*this, "Error during initial project load !\n\t", nana::msgbox::button_t::ok)
                            .icon(nana::msgbox::icon_information )
                        << e.what()    << "\n\n A new Default Project will be created. "
                        ).show (  ) ;
		save_defPr() ; 					                
    }

	//this->comBoxTAMeth->SelectedIndex  = SMStLucia;    
    try{ 
		    _cp.Actualize_All_NNp();
            LoadSequences();
		}
    catch ( std::exception& e )      //  
	{   
        (nana::msgbox(*this, "Error during sequence or NN parameter load !\n\t", nana::msgbox::button_t::ok)
                            .icon(nana::msgbox::icon_information )
                        << e.what()     
                        ).show (  ) ;
		save_defPr() ; 					                
    }

	mExpl_.InitTree();


	try {
		InitMyLayout();
	}
	catch (std::exception & e)
	{
		throw std::runtime_error(std::string("An error occurred during initialization of the windows layout\n") + e.what());
	}
	catch (...)
	{
		throw std::runtime_error(std::string("An unknown error occurred during initialization of the windows layout\n"));
	}

	//catch (std::exception& e)      //  
	//{
	//	
	//	(nana::msgbox(*this, "Error during Initialization of the windows layout\n\t", nana::msgbox::button_t::ok)
	//		.icon(nana::msgbox::icon_information)
	//		<< e.what()
	//		).show();
	//	save_defPr();
	//}

    setup_.AddMenuItems (_menuBar.push_back("P&roject"));     // 0
    mExpl_.AddMenuItems (_menuBar.push_back("&Sequences"));   // 1 
    AddMenuProgram();                                         // 2
	_menuBar.at(2).append_splitter();
	_menuBar.at(2).append("&About", [&](nana::menu::item_proxy& ip) 
	{
        try {
            About ab;
            ab.modality();
        }
        catch (std::exception e)
        {
            (nana::msgbox(*this, "Error in About !\n\t", nana::msgbox::button_t::ok)
                    .icon(nana::msgbox::icon_information )
                    << e.what()
            ).show (  ) ;
        }
		//ab.show();
	});

    SelectClickableWidget( _menuBar);

	//nana::fn_group<void(tabbar&, value_type&)> active;

	tabbar_.events().activated( [this]()// const tabbar & t, const tabbar::value_type& tab)
	{
		bool enable	= tabbar_.activated( )==1;
		auto &m		= _menuBar. at(1);
		auto sz		= m.size();

		for(int i=0; i< sz; ++i)
			m.enabled(i,enable);
	});


}
catch (std::exception & e)
{ 
	throw std::runtime_error(std::string("An error occurred during initialization of the main window of the application ThDy DNA Hybrid:\n") + e.what() );
}
catch (...)
{ 
	throw std::runtime_error(std::string("An unknown error occurred during initialization of the main window of the application ThDy DNA Hybrid") );
}

    //~ThDyNanaForm();
void     ThDyNanaForm::SetDefLayout   () 
{
    _DefLayout=  
	    "vertical      gap=2                                      	  \n\t"
	    "	 vertical      gap=2                   		              \n\t"
	    "		        <weight=25   reserved_for_menubar >                   		\n\t"
	    "		        <weight=23   PagesTag     >      		      \n\t"
	    "		        <min=255     Pages        >      		      \n\t"
	    "		 		                                              \n\t"
	    "		                                                       \n\t"
        ;
}

void     ThDyNanaForm::AsignWidgetToFields() 
{
	_place.field("PagesTag")        << tabbar_  ;
}   

void     ThDyNanaForm::add_page(widget& w)
{   
    //nana::pixel_rgb_t bk;
    //bk.u.color = background ();
    //bk.u.element.blue =0; 
    //w.background (1);
    tabbar_.push_back (                    w.caption());
    tabbar_.attach    (tabbar_.length()-1, w          );
	_place.field("Pages"   ).fasten( w)  ;
}

//std::ifstream in("in.txt");
//std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
//std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!
//
//std::ofstream out("out.txt");
//std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
//std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
//
//std::string word;
//std::cin >> word;           //input from the file in.txt
//std::cout << word << "  ";  //output to the file out.txt
//
//f(); //call function
//
//
//std::cin.rdbuf(cinbuf);   //reset to standard input again
///\todo What about custom colors in nana GUI?
//nana::color::current_schema[nana::color::schema::list_header_border]=nana::color::Red;
//nana::color::current_schema[nana::color::schema::list_header_bg]=nana::color::Yellow;    // 0xF1F2F4
//nana::color::current_schema[nana::color::schema::list_header_highlighted_bg]=nana::color::Rose;    // 0xFFFFFF
//nana::color::current_schema[nana::color::schema::list_header_pressed_bg]=nana::color::AliceBlue;
//nana::color::current_schema[nana::color::schema::list_header_grabed_bg]=nana::color::Ash_Gray;    // 0x8BD6F6
//nana::color::current_schema[nana::color::schema::list_header_floated_bg]=nana::color::Aztech_Purple;	   // 0xBABBBC
//nana::widget_colors::background ;   ::current_schema[nana::color::schema::list_header_floated_bg]=nana::color::Aztech_Purple;	   // 0xBABBBC
//nana::form::scheme_type().background = nana::colors::light_sky_blue;