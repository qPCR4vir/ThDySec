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
inline TmGPos Hybrid(CSec &primer, CSec &target, 	ThDyAlign &Al,	ofstream &osTm,
														ofstream &osG,
														ofstream &osPos,
														ofstream &osPl_Tm,
														ofstream &osPl_G,
														ofstream &osAl,
														Table    *rtbl 	/*,
														CTable<Temperature> *tlTm,
														CTable<Energy>	*tlG,
														CTable<SecPos> *tlPos*/)
{	Al.Align( &primer, &target);					//  virtual !!!
	Al.SelectOptParam(Al.Ta());				//	FrAl.GetOptHit();

	TmGPos tmgp ( KtoC( Al.Tm() ), Al.G()	/ 1000	,Al._maxgloj);
	if (osTm	) osTm		<<sep	<<	KtoC(Al.Tm())			;	if (rtbl) *rtbl << tmgp;
	if (osG		) osG		<<sep	<<		 Al.G()	/ 1000		;	/*if (tlG)  *tlG  <<  Al.G()	/ 1000		;*/
	if (osPos	) osPos		<<sep	<<		 Al._maxgloj		;	//if (tlPos) *tlPos<<  Al._maxgloj			;	// pos del 5'
	if (osPl_Tm	) osPl_Tm	<<"\t"	<<  KtoC(Al.Tm())			;
	if (osPl_G	) osPl_G	<<"\t"	<<  	 Al.G()	/ 1000		;

	print_ThDyAlign (osAl, Al);
	//Al.Export_DPMz_Pre(osAl);

    return tmgp;
}

/// Align all CSec vs all in the list os sequences probe and target in an existing ThDyAlign with output a line of Tm, G and Position to a table and files.
void HybridPr(CMultSec &primers, CSec &target, 	ThDyAlign &Al,	ofstream &osTm,
														ofstream &osG,
														ofstream &osPos,
														ofstream &osPl_Tm,
														ofstream &osPl_G,
														ofstream &osAl,
														Table *rtbl 	/*,
														CTable<Temperature> *tlTm,
														CTable<Energy>	*tlG,
														CTable<SecPos> *tlPos*/)
{
    // todo take 2 tables - the second for one value per primer.NonDegSet, not per NonDeg variant.

    // recorre todos los primers de nuevo
	for (auto& uptr_primer : primers.SecL())
	{
		CSec &primer = *uptr_primer;
		if(!primer.Selected()) continue;
		if ( ! primer.NonDegSet()  )
			{	
				Hybrid (primer, target, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
			} else 
				{	
                    CMultSec *nds=primer.NonDegSet().get() ;

                    // todo record were it will begin to write to the table
                    // auto &Num = rtbl->Next();
                    // auto [rrow, rcol] = rtbl->get_pos();
                    // index row = rrow, col = rcol;
                    TmGPos tmgp_final; bool first = true;
					for ( auto& nds_CurSec :  nds->SecL() ) // recorre todos las var no deg
					{	
                        CSec &ndeg_primer_variant = *nds_CurSec ;
                        TmGPos tmgp = Hybrid (ndeg_primer_variant, target, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
                        if (first)
                        {
                            tmgp_final = tmgp;
                            first = false;
                        }
                        else
                        {
                            tmgp_final._G  = std::min(tmgp_final._G, tmgp._G);
                            tmgp_final._Tm = std::max(tmgp_final._Tm, tmgp._Tm);
                        }
					}// recorre todos las var no deg

					// todo squash all the previous into one row.
					*rtbl << tmgp_final;
				} 
	}
	for (auto& CurMSec : primers.MSecL())
        HybridPr (*CurMSec, target, 	Al, osTm, osG,osPos,osPl_Tm,osPl_G,osAl, rtbl);
}			
