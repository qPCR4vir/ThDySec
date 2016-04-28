/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDy_DNAHybrid.Nana\uArray.cpp
*
* @brief 
*
*/


#include "ThDy_DNAHybrid.Nana\uArray.h"
#include "ThDy_DNAHybrid.Nana\TableResults.h"
#include "ThDy_DNAHybrid.Nana\main.Nana.h"






MplexPCR::MplexPCR            (ThDyNanaForm& tdForm) 
try     : _Pr             (tdForm), 
          CompoWidget     (tdForm, ("MplexPCR"), ("MplexPCR.lay.txt"))
    {

        _do_mPCR      .events().click([&](){buttPCR_Click ();});

        InitMyLayout();
        SelectClickableWidget( *this);
    }
catch (std::exception & e)
{
	throw std::runtime_error(std::string("An error occurred during initialization of the MplexPCR page window:\n") + e.what());
}
catch (...)
{
	throw std::runtime_error(std::string("An unknown error occurred during initialization of the MplexPCR page window"));
}

void MplexPCR::buttPCR_Click()  //	  Run      _IPrgPar_mPCR
{	 			
	try{                                   
		_Pr._mPCR._cp.Actualice_NNp();  
 		_Pr.Run(_Pr._mPCR);	

        _Pr.mExpl_.RefreshProbes_mPCR(/*false*/); 

        _Pr._results.emplace_back(new TableHybRes(_Pr._mPCR._rtbl));
        _Pr._results.back()->show();
        _Pr._results.emplace_back(new TableHybRes(_Pr._mPCR._rtbl_self));
        _Pr._results.back()->show();
	}
catch ( std::exception& e)
	{ 
		std::cerr<< e.what()    ;
        (nana::msgbox(*this,("Error during multiplex PCR analysis !"), 
                                            nana::msgbox::button_t::ok)   <<e.what()) (  ) ;
		return;
	}
    //ShowResTbl(_Pr._mPCR._rtbl );
    //_Pr._uArr._rtbl = nullptr;

    //ShowResTbl(_Pr._mPCR._rtbl_self );
    //_Pr._mPCR._rtbl_self = nullptr;

}

   uArray::uArray            (ThDyNanaForm& tdForm) try
        : _Pr             (tdForm), 
          CompoWidget     (tdForm, ("uArray"), ("uArray.lay.txt"))
    {

        _do_uArray      .events().click([&](){buttuArray_Click ();});

        InitMyLayout();
        SelectClickableWidget( *this);
    }
   catch (std::exception & e)
   {
	   throw std::runtime_error(std::string("An error occurred during initialization of the uArray page window:\n") + e.what());
   }
   catch (...)
   {
	   throw std::runtime_error(std::string("An unknown error occurred during initialization of the uArray page window"));
   }


       void uArray::buttuArray_Click()  
	{	 			
	 try{                                   
		  _Pr._uArr ._cp.Actualice_NNp();  
 		  _Pr.Run(_Pr._uArr);	

          _Pr._results.emplace_back(new TableHybRes(_Pr._uArr._rtbl));
          _Pr._results.back()->show();
		}
	catch ( std::exception& e)
		{ 
          std::cerr<< e.what()    ;
          (nana::msgbox(*this,("Error during uArr analysis !"), 
                                                nana::msgbox::button_t::ok)   <<e.what()) (  ) ;
		  return;
		}
	}
