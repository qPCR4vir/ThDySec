/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\src\ThDySec\sec_rang.cpp
*
* @brief 
*/

#ifdef WINDOWS_FORM_GUI
#include "stdafx.h"
#pragma unmanaged
#endif

//#include <stdio.h>
//#include <string.h>
//#include <assert.h>
//#include <iostream>
//#include <iomanip>
//#include <memory>
//#include <math.h>
//#include <list>
//#include <stack>

using namespace std ; 

#include "ThDySec/sec_mult.h"
#include "ThDySec/sec_rang.h"
#include "ThDySec/common.h" 
using namespace DegCod;

CSecAl::CSecAl(CSec &sec, long LenAlign)
		: _Sec(sec), 
		  _inAlp_B(sec.Len()), 
		  _inBp_Al(LenAlign)
	   {}
CSecAl::CSecAl(CSec &sec) 
		: _Sec(sec), 
		  _inAlp_B(sec.Len())
	  {}

//!  Inicializa array de posibles cand en esta sec. ;---------------------  CSec------------>  CSecCand  ---------------------
		CSecCand::CSecCand(CSec &sec, 	SondeLimits sL)		  	
					:	_Sec(sec), 	
						_rg(sec.Len()+2), 	/// \todo REVISE!! vector of empty unique ptr to rang	 
						_NumPosCand(0), 
						_NumCand(0),
						_NumCandExact(0)							

{	
	//	assert (_rg);	trow exeption !!
	//for (fi=0; fi< sL._L.Min() ; fi++)   ;		/// \todo REVISE!!         // se salta las primeras pos
														// al comienzo fi = L_min	
	for (long fi=sL._L.Min()	; fi<=sec.Len(); fi++) // fi - final base of candidate, recorre toda la sec
	{	
		long        pi = fi - sL._L.Max() +1	; if (pi < 1 ) pi=1;
		CRangBase R(pi	,fi - sL._L.Min() +1 )  ;				assert( fi>R.Max() );		assert( pi<=R.Max() );

		for (long i0=R.Min() ; i0<=R.Max() ; ++i0 )
		{	if ( sL._Tm.inRang( sec.Tm(i0,fi) )    &&    sL._G.inRang( sec.G(i0,fi) ) ) 	// usa calculos de Tm basados en NNpar
			{	 R.adjustCur(i0);
				 _NumCandExact++;	//	Seria lo correcto, pero da problemas cuando existe una "burbuja" de Tm, 
									//  cosa que un "rango" no puede considerar, tendria que hacerce con un conjunto. 
									//  Osea estoy dejando dentro del rango, interiormente la posibilidad de aceptar sondas 
									//	con Tm fuera de rango
			}
		}
		if (R.hasMatch())
		{		_NumPosCand++ ; 
				_NumCand+= R.NumMatch();				
				_rg[fi].reset(new CRang (R.MatchRange()));    /// \todo REVISE!! 				
                                                    assert( fi>R.Max() );	assert( fi>_rg[fi]->Max() );	
		} // else	_rg[fi]=0 ;
	}
	_NumPosCandIn=_NumPosCand ; _NumCandIn=_NumCand ;
}

//NumRang<long> p(fi-sL._L.Max() +1 , fi-sL._L.Min() +1 ) , pcur;
		//long pi		=fi-sL._L.Max() +1		, pf		=fi-sL._L.Max() +1 ;   // al comienzo pf = L_max -L_min +1
		//long picur	=pf+1					, pfcur		=pi	;
//	assert(((clog<< "\nCreating: "<<sec._name<<" \t,#pos: \t"<< _NumPosCand<<", #cand: \t"<< _NumCand),1));
							//if(picur>i0)					picur=i0;
							//if(pfcur<i0)					pfcur=i0;
//(R._pfcur - R._picur + 1);




/// Despues de comparar dos seq "candidatos" analiza los rangos de una y los colapsa 
long	CSecCand::ColapseRangs(bool colapse) // hacer otra variante para "no colapse" rang by self hybri - no sec str!
{								// _NumPosCand = 0 ;
								// anadir 2 parametros: pos de com y fin de zona "efectiva" o cubierta por sec pareada
								// fuera de esa zona no colapsar
	long NumPosCand(0), NumCand(0); // 2011-05-16. Para resolver el problema del conteo de pos y cand.

											//    pi      pf                fi
											//----|++++++++|-----------------|--------			El rango inicial, y como va quedando
											// pfcur                        fi					El rango para calculo ("cur"), antes del comienzo	
											//---|----------|----------------|--------			Asi se queda si no hibridan entre si las sec en esta zona,
											//             picur            fi					y entonces "colapsa" el rango
											//            pfcur             fi					En este caso encontro 5 "cand" comunes	
											//--------|+++|------------------|--------			
											//       picur                  fi					
	
	for (long fi=0; fi<=_Sec.Len(); fi++)
	{	if (!_rg[fi]) continue ;		 
											
		CRang &R (*_rg[fi]);				assert(fi>R.Max());

		if ( R.hasMatch()  )	
		{	
			NumPosCand++ ; NumCand += R.NumMatch() ;						
			if (colapse) 
			{	
				R.SchrinkToMatch();
			}
			R.IncrMatchs();
			R.open();
		}else															
		{
			if (colapse) {	_rg[fi].reset() ;} 	
		}
	}

	// if (colapse)											//   REVISAR    !!!!!!!!!!!
	{	/*_NumPosCandIn =*/ _NumPosCand = NumPosCand ; 			
		/*_NumCandIn    =*/ _NumCand    = NumCand ;   
	}

	return _NumCand ;
}


// _NumPosCand-- ;	//	assert(((clog<<"\n("<<fi<<"->"<<_rg[fi]->_pi<<"-"<<_rg[fi]->_pf<<")
//	---- colapsing:"<<_rg[fi]->_picur<<"-"<<_rg[fi]->_pfcur),1));
// amplitud : pfcur - picur
//	En realidad solo puede ser =0 y no <0, o >0 claro
			//for (int pi_pos=R._pi    ; pi_pos <= R._picur ; pi_pos++ ) 
			//R.match[pi_pos - R._pi]++;
			//for (int pf_pos=R._pfcur ; pf_pos <= R._pf    ; pf_pos++ ) 
			//R.match[pf_pos - R._pi]++;
// en la ultima comparacion no se confir cand en esta fi	//_NumPosCand--;
			//	assert(((clog<< "\nColapsing "<<_Sec._name<<" \t,#pos\t"<< _NumPosCand<<"\t, #cand: "<< _NumCand),1)); ;


char *	CSecAl::CopyAlignedSecChar(long Al_pBeg, long Al_pEnd, char *CharSec)	// CUIDADO !! asume suficiente espacio !!
{	//assert (Al_pBeg<=Al_pEnd);																			// "EXPERIMENTAL" ---------------------
	assert (Al_pBeg>0); //assert (Al_pEnd>0); 
	assert (CharSec);
	long j, i;
	for ( i=Al_pBeg, j=0; i<=Al_pEnd; i++,j++)
		if (_inAlp_B[i] == -1)	CharSec[j] = nu2dgba[0] ;
		else					CharSec[j] =_Sec[_inAlp_B[i]] ; // Para las incerciones uso -1 en el Al
							
	CharSec[j]=0;
	return CharSec;
}
char *	CSecAl::GetAlignedSecChar(long Al_pBeg, long Al_pEnd)  // "regala" esta memoria, no olvide delete !
{	long l=Al_pEnd-Al_pBeg+1;																				// "EXPERIMENTAL" ---------------------
	assert (l>0);
	char *Al_sec= new char[l+1];
	assert (Al_sec);
	return CopyAlignedSecChar(Al_pBeg,  Al_pEnd, Al_sec);
}

