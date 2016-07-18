/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file ThDySec\src\ThDy_programs\ThDy_MultiplexPCR.cpp
*
* @brief
*
*/

//#include "StdAfx.h"
#pragma unmanaged
//#include "ThDySec\th_dy_align.h"
#include "ThDy_programs/prog_comm_functions.h"

int microArrayProg ( CProgParam_microArray *IPrgPar_uArr, 
                    CMultSec &pr, 
                    CMultSec &tg, 
                    time_t t_0,  
                    int MAxGrDegTg, 
                    const std::string& of_x =""	);

void CreateComplProbes(	CMultSec		&pr	)   /// revise !!!!!!!!!!!!
{
	CMultSec *cms{nullptr};

	for (auto &CurMSec : pr.MSecL())   ///\todo recursive ?! use only selected   ?!!!
	{
		if (CurMSec->_name == "compl")        // if there was already one compl group
		{
			cms = CurMSec.get();
			cms->clear();                ///\todo uncount !!!!!        // destroy to reuse  !!
		}
		else if (CurMSec->Selected())         // recursively CreateComplProbes for all the other groups
			CreateComplProbes(*CurMSec);
	}

	if(!cms) 
        cms=pr.AddMultiSec("compl");

	for (auto &CurSec : pr.SecL()) 			// recorre todos las sondas
	{	CSec &s = *CurSec ;
		if( s.Selected())
		    cms->AddSec ( s.Clone(DNAstrand::rev_compl) ); 
	}
}


		//int t=MultiplexPCRProg ( IPrgPar_Calc, primers		)  ;

/// main function
int MultiplexPCRProg ( CProgParam_MultiplexPCR *IPrgPar_uArr, CMultSec		&pr)  
{
	time_t t_0 = time(nullptr);
    IPrgPar_uArr->_cp.Check_NNp_Targets ();
	CreateComplProbes(	pr	);

	          microArrayProg   ( IPrgPar_uArr, 
                                 pr	, 
                                 pr, 
                                 t_0, 300    , 
                                 "_self"	
                               )  ;

   IPrgPar_uArr->_rtbl_self = IPrgPar_uArr->_rtbl;
   IPrgPar_uArr->_rtbl->TitTable( IPrgPar_uArr->_cp._OutputFile.get()  + ": Primers / Primers (align method: " 
	                             +IPrgPar_uArr->_cp.TAMeth.ToString() +  " ). Multiplex PCR.") ;

	auto res= microArrayProg ( IPrgPar_uArr, 
                               pr	, 
                               *IPrgPar_uArr->_cp._pSeqTargets.get(), 
                               t_0  , 300 
                              )  ;

    IPrgPar_uArr->_rtbl->TitTable(IPrgPar_uArr->_rtbl->TitTable()+ ". Multiplex PCR.");
	return res ; 
}


	int MultiplexPCRProg ( CProgParam_MultiplexPCR *IPrgPar_uArr)  
{

    IPrgPar_uArr->Check_NNp_Targets_probes (IPrgPar_uArr->_probesMS.get());

	return MultiplexPCRProg ( IPrgPar_uArr, *IPrgPar_uArr->_probesMS.get())  ;

	
}

