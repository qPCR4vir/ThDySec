/**
* Copyright (C) 2009-2017, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\main.Nana.cpp
*
* @brief Entry point (main) for ThDyHybr with Nana GUI
*
*
*/

#include "ThDy_DNAHybrid.Nana\main.Nana.h"

 std::string e_mail_firma = " arielvina@yahoo.es";

//if you want to keep the Windows subsystem you can just hint at what your entry point is, 
//because you haven't defined ___tmainCRTStartup. You can do this by adding the following 
//to Properties -> Linker -> Command line:
//
//    /ENTRY:"mainCRTStartup"
//
//This way you get rid of the console window.

 class About : public nana::form, public EditableForm
 {
	 nana::label  copy_r{ *this, R"(               Copyright (C) 2009-2018, Ariel Vina-Rodriguez (qPCR4vir)
                                  (  arielvina@yahoo.es  )
                                  http://qpcr4vir.github.io/)" };


	 nana::label  comments{ *this,R"(This work is mentioned in my PhD thesis at INNT-FLI.

		 Program distributed under the GNU General Public License, see:
 http://www.gnu.org/licenses/)" };
	 
	 

	 nana::label  compiled{ *this, R"(Compiled on:   )"  __DATE__   R"( / )"  __TIME__   R"(    Version: v0.01.04

   Downloads and source code: https://github.com/qPCR4vir/ThDySec
                       Wiki: https://github.com/qPCR4vir/ThDySec/wiki
   ________________________________________________________________________________
   Powered by Nana C++ GUI library:  http://nanapro.org/en-us/
                                  Wiki:  https://github.com/qPCR4vir/nana-docs/wiki
      Nana Version: 1.6.2 cmake-dev : https://github.com/qPCR4vir/nana/
         )" };

	 
 nana::button close { *this, "Close" };
	 

 public:
	 About() : nana::form   ( nana::rectangle(nana::point(50, 5), nana::size(500, 350)) ),
	           EditableForm ( nullptr, *this, "About ThDy Hybrid" , "about.lay.txt") 
	 {
		 InitMyLayout();
		 SelectClickableWidget(copy_r);

	 }

 void SetDefLayout() override
	 {
		 _DefLayout = 	 "vertical      gap=2 margin=2    	all"	 ;
	 }
 void AsignWidgetToFields()
 {
	 _place.field("all") << copy_r << comments << compiled << close ;
 }
 };


int main(int argc, char *argv[]) 
{

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
  try	
  {
    using namespace ParamGUIBind;

    IParBind::SetDef(PriorizeDefault::Parametr );
    ThDyNanaForm tdForm(  argc,  argv);
    //tdForm.ReCollocate();
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
                   EditableForm    (nullptr, *this,  "ThDy DNA Hybrid" " (" __DATE__ " / " __TIME__")" , "ThDy.lay.txt")  
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

    setup_._proj.FileNameOnly( ProjetFile()  );
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

		About ab;
		ab.show();


		(nana::msgbox(this->handle(), "About ThDy Hybrid", nana::msgbox::button_t::ok) <<
			R"(               Copyright (C) 2009-2018, Ariel Vina-Rodriguez (qPCR4vir)
                                  (  arielvina@yahoo.es  )
                                  http://qpcr4vir.github.io/

   This work is mentioned in my PhD thesis at INNT-FLI.

   Program distributed under the GNU General Public License, see:
           http://www.gnu.org/licenses/

   Compiled on:   )"  __DATE__   R"( / )"  __TIME__   R"(    Version: v0.01.04

   Downloads and source code: https://github.com/qPCR4vir/ThDySec
                       Wiki: https://github.com/qPCR4vir/ThDySec/wiki
   ________________________________________________________________________________
   Powered by Nana C++ GUI library:  http://nanapro.org/en-us/
                                  Wiki:  https://github.com/qPCR4vir/nana-docs/wiki
      Nana Version: 1.6.2 cmake-dev : https://github.com/qPCR4vir/nana/
         )") ();
			;
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
 
