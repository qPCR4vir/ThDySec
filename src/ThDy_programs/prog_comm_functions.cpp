/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\src\ThDy_programs\prog_comm_functions.cpp
*
* @brief A few global functions for hybridization
*
*/

//#include "StdAfx.h"
#pragma unmanaged
#include "ThDy_programs/prog_comm_functions.h"
using namespace std;




unique_ptr<ThDyAlign> Create_ThDyAlign(const ThDyCommProgParam& _cp, LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar)
{
	unique_ptr<ThDyAlign>	apAl;
	switch (	_cp._TAMeth )
	{	case TAMeth_Tm: default:	apAl.reset(new ThDyAlign_Tm( MaxLenSond ,  MaxLenTarg, NNpar) );   break;
		case TAMeth_Fract:			apAl.reset(new FracTDAlign ( MaxLenSond ,  MaxLenTarg, NNpar) );   break;
		case TAMeth_G:				apAl.reset(new ThDyAlign_G ( MaxLenSond ,  MaxLenTarg, NNpar) );   break;
	}

	apAl->SetTa		 (			CtoK(	_cp._Ta	    ));		// OK
	return apAl;
}

/// Set CSec and trigget the aligment in an existing ThDyAlign with output a line of Tm, G and Position to a table and files.
inline void Hybrid(CSec &s, CSec &t, 	ThDyAlign &Al,	ofstream &osTm,
														ofstream &osG,
														ofstream &osPos,
														ofstream &osPl_Tm,
														ofstream &osPl_G,
														ofstream &osAl,
														Table    *rtbl 	/*,
														CTable<Temperature> *tlTm,
														CTable<Energy>	*tlG,
														CTable<SecPos> *tlPos*/)
{	Al.Align( &s, &t);					//  virtual !!!
	Al.SelectOptParam(Al.Ta());				//	FrAl.GetOptHit();

	if (osTm	) osTm		<<sep	<<	KtoC(Al.Tm())			;	if (rtbl) *rtbl << TmGPos ( KtoC( Al.Tm() ), Al.G()	/ 1000	,Al._maxgloj);
	if (osG		) osG		<<sep	<<		 Al.G()	/ 1000		;	/*if (tlG)  *tlG  <<  Al.G()	/ 1000		;*/
	if (osPos	) osPos		<<sep	<<		 Al._maxgloj		;	//if (tlPos) *tlPos<<  Al._maxgloj			;	// pos del 5'
	if (osPl_Tm	) osPl_Tm	<<"\t"	<<  KtoC(Al.Tm())			;
	if (osPl_G	) osPl_G	<<"\t"	<<  	 Al.G()	/ 1000		;

	print_ThDyAlign (osAl, Al);
	//Al.Export_DPMz_Pre(osAl);


}

/// Align all CSec vs all in the list os sequences probe and target in an existing ThDyAlign with output a line of Tm, G and Position to a table and files.
void HybridPr(CMultSec &pr, CSec &t, 	ThDyAlign &Al,	ofstream &osTm,
														ofstream &osG,
														ofstream &osPos,
														ofstream &osPl_Tm,
														ofstream &osPl_G,
														ofstream &osAl,
														Table *rtbl 	/*,
														CTable<Temperature> *tlTm,
														CTable<Energy>	*tlG,
														CTable<SecPos> *tlPos*/)
{	// recorre todos los primers de nuevo
	for (auto& CurSec : pr.SecL())
	{
		CSec &s = *CurSec;
		if(!s.Selected()) continue;
		if ( ! s.NonDegSet()  ) 
			{	
				Hybrid (s, t, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
			} else 
				{	
                    CMultSec *nds=s.NonDegSet().get() ;
					for ( auto& nds_CurSec :  nds->SecL() ) // recorre todos las var no deg
					{	
                        CSec &s = *nds_CurSec ;
						Hybrid (s, t, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
					}// recorre todos las var no deg
				} 
	}
	for (auto& CurMSec : pr.MSecL())  
        HybridPr (*CurMSec, t, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
}			
