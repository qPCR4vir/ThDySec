/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2019
*
* @file ThDySec\src\ThDy_programs\ThDy_MultiplexPCR.cpp
*
* @brief
*
*/

#include "ThDy_programs/prog_comm_functions.h"

int microArrayProg( CProgParam_microArray *IPrgPar_uArr,
                    CMultSec &pr, 
                    CMultSec &tg, 
                    time_t t_0,  
                    int MAxGrDegTg, 
                    const std::string& of_x =""	);

void CreateComplProbes(	CMultSec		&pr	)   /// revise !!!!!!!!!!!!
{
    auto &ms = pr.MSecL();
    for ( auto CurMSec_i = ms.begin(); CurMSec_i != ms.end(); CurMSec_i++)
    ///\todo recursive ?! use only selected   ?!!!
	{
        CMultSec *CurMSec = CurMSec_i->get();
        if (CurMSec->_name == "compl")        // if there was already one compl group
		{
            CurMSec->_parentMS->DeleteMSec(CurMSec_i);
		}
		else if (CurMSec->Selected())         // recursively CreateComplProbes for all the other groups
			CreateComplProbes(*CurMSec);
	}

	auto cms = *pr.AddMultiSec("compl");

	for (auto &CurSec : pr.SecL()) 			// recorre todos las sondas
	{	CSec &s = *CurSec ;
		if( s.Selected())
            cms->AddFreeSec(std::shared_ptr<CSec>(s.Clone(DNAstrand::rev_compl)));
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

   IPrgPar_uArr->_rtbl->TitTable( IPrgPar_uArr->_cp._OutputFile.get()  + ": Primers / Primers (align method: "
	                             +IPrgPar_uArr->_cp.TAMeth.ToString() +  " ). Multiplex PCR.") ;
    IPrgPar_uArr->_rtbl_self = std::move(IPrgPar_uArr->_rtbl);

	auto res= microArrayProg ( IPrgPar_uArr, 
                               pr	, 
                               *IPrgPar_uArr->_cp._pSeqTargets, 
                               t_0  , 300 
                              )  ;

    IPrgPar_uArr->_rtbl->TitTable(IPrgPar_uArr->_rtbl->TitTable()+ ". Multiplex PCR.");
	return res ; 
}


	int MultiplexPCRProg ( CProgParam_MultiplexPCR *IPrgPar_uArr)  
{

    IPrgPar_uArr->Check_NNp_Targets_probes (IPrgPar_uArr->_probesMS-get());

	return MultiplexPCRProg ( IPrgPar_uArr, *IPrgPar_uArr->_probesMS)  ;

	
}

