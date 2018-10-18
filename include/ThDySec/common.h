/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  C:\Prog\ThDySec\include\ThDySec\common.h
*
* @brief Basic definitions and declarations. 
* Commonly used by all: user interface and programs. Almost primary, only depends on ProgParam\include\common_basics.h for NumRang<>.
*/

#ifndef _COMMON_H
#define _COMMON_H

#pragma unmanaged	
#pragma warning( disable : 4996 )

#include <cctype> 
#include <string>
#include <sstream> 
#include <stdexcept>
//#include "common_basics.h"
#include "..\..\ProgParam\include\common_basics.h"

template <class Num>
inline bool IsEq (Num f1, Num f2, float prec= 0.01f )	// ; bool	IsEq (float f1, float f2, float prec)   ----> dif/med < prec ????
{	float fm = float( abs(f1) + abs(f2) )/ 2.0f ;		if (fm == 0) return true ;
	return abs(float(f1-f2)) / fm  < prec 	;
}

using Energy        = float	  ;	// ------------   	Energy	--	Entropy	--	Temperature	--	SecPos --	LonSecPos and Rang s
using Entropy       = Energy  ;
using Temperature   = float	  ;
using SecPos        = int	  ;
using LonSecPos     = long	  ;


using  EnergyRang       = NumRang<Energy>       ;
using  EntropyRang      = NumRang<Entropy>		;
using  TemperatureRang  = NumRang<Temperature>  ;
using  SecPosRang       = NumRang<SecPos>		;
using  LonSecPosRang    = NumRang<LonSecPos>	;



extern char sep[];
enum SaltCorrection {NoSelect  =-1 , StLucia  =0 , Owczarzy  =1 }; ///< ya se puede usar StLucia inicializando todo en el constructor. Parcialmente implementado cambio de Conc
enum SaltCorrecP	{NoSMSelect=-1 , SMStLucia=0 , SMOwczarzy=1 }; ///< \todo Unificar !!

enum AlignMeth		{TAMeth_Tm =0  , TAMeth_G    , TAMeth_Fract };

enum class	 DNAstrand		{plus	, minus, direct	, rev	, complem, rev_compl	} ;
constexpr const char *DNAstrandName[]=	{""		, "(c)", ""		, "(r)"	, "(i)", "(c)"		} ;

enum  fileFormat   {fasta =1 , csv =1u<<1, f2 =1u<<2, f3 =1u<<3} ; ///< for OUTPUT !!!! they can be combined with "OR"


struct TmGPos 
{
    Temperature _Tm; 
    Energy      _G; 
    SecPos      _Pos;

	TmGPos(Temperature Tm, Energy G, SecPos Pos)
        :_Tm(Tm), _G(G), _Pos(Pos) 
    {}

	TmGPos() = default ;
};

extern const 	EnergyRang		G_def  ;	//  G_def (-5,-1), 
extern const 	TemperatureRang	Tm_def ;	//  Tm_def(57,63) ;
extern const 	SecPosRang		L_def  ;	//  L_def (20,35) ;

class SondeLimits		// ------------------------------------------------------------------	SondeLimits	-------------------------------
{public: 
	EnergyRang		_G;
	TemperatureRang	_Tm ;
	SecPosRang		_L ;
	SondeLimits():_G(G_def),_Tm(Tm_def),_L(L_def){}
	SondeLimits(EnergyRang G, TemperatureRang Tm, SecPosRang L):_G(G),_Tm(Tm), _L(L){}
};
inline Temperature	  CtoK	(Temperature Tm)	{return Tm + Temperature(273.15f) ;}
inline Temperature	  KtoC	(Temperature Tm)	{return Tm - Temperature(273.15f) ;}

inline TemperatureRang CtoK	( TemperatureRang T  ) {	     T.Set( CtoK (     T.Min()) , CtoK (     T.Max()) ) ; return T ;}
inline TemperatureRang KtoC	( TemperatureRang T  ) {	     T.Set( KtoC (     T.Min()) , KtoC (     T.Max()) ) ; return T ;}

inline SondeLimits    CtoK	( SondeLimits    sL ) { sL._Tm.Set( CtoK (sL._Tm.Min()) , CtoK (sL._Tm.Max()) ) ; return sL;}
inline SondeLimits    KtoC	( SondeLimits    sL ) { sL._Tm.Set( KtoC (sL._Tm.Min()) , KtoC (sL._Tm.Max()) ) ; return sL;}
inline SondeLimits    convCtoK_ctok	( SondeLimits    sL ) { sL._Tm.Set( CtoK (sL._Tm.Min()) , CtoK (sL._Tm.Max()) )		;	
															sL._G.Set ( sL._G.Min()*1000    ,  sL._G.Max()*1000   )		; return sL;}
inline SondeLimits    convKtoC_ktoc	( SondeLimits    sL ) { sL._Tm.Set( KtoC (sL._Tm.Min()) , KtoC (sL._Tm.Max()) )		;	
															sL._G.Set ( sL._G.Min()/1000    ,  sL._G.Max()/1000   )		; return sL;}


#endif 