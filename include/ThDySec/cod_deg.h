/**
* Copyright (C) 2009-2016, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDySec\cod_deg.h
*
* @brief  to avoid billions of repited conversions nucleotide letter - code and make ease find complement, etc.
*
* se usa asi: 
* \code{.cpp}
* for (Base b=0; b <  n_dgba				; b++)			// recorre las bas deg, aqui solo para
*				                                        	// mostrar todas las variante
*
*     for (Base c=0; c <= grad_deg[ nu2dgba[b] ]  ; c++)	// recorre las bas no deg de la base deg b
*     {		dg2ban[b][c]  ;    					            // cod corto de cada caracter a generar
*   		dg2ba [b][c]  ;						            // caracter a generar
*           db2nu[dg2ba [b][c]] ;						    // cod deg (largo) de cada caracter a generar				
*           b & db2nu[dg2ba [b][c]] ;	
*   		dg2bkn[b][c]  ;						            // cod K de cada caracter a generar
*     }
* \codeend
*/

#ifndef _COD_DEG_H
#define _COD_DEG_H
#include <limits.h>
#include <ctype.h>
#pragma warning( disable : 4996 ) // ?

/// Use schort nt codes to avoid billions of repited conversions nucleotide letter <-> code and make ease to construct complement, etc.
namespace DegCod
{
using Base = unsigned char ;///< A simple char can be neagtive which can be a problem when using as index for an array.
using Code = Base;
inline Base base(char b){return static_cast<Base>(b);}

///\todo ? Introducir Letter->UCHAR_MAX , Base -> "ACGT<>-" , DegBase -> "-GCSTKYBARMVWDHN"
///\todo ?  las 4 bases en el orden complementario de Kaderali. Cambiar a algo asi "-TGCA." o "GCTA-." o "GCTAx-." Calcular si no es mejor "-GCSTKYBARMVWDHN"
///\todo ?  usar basek[]=".ACGT$"   or basek[]="-ACGT$"  --- efecto del '-' ????

constexpr 
    Base	nu2ba		[]="GCTA"				,  ///< las 4 bases, cod corto. de "numero a base".
			nu2c_ba		[]="CGAT"				,  // 
			n_ba		  =sizeof(nu2ba)-1      ,			

			basek		[]="-ACGT$"				,  ///< las 4 bases en el orden de Kaderali. ? Coservar este orden  => b+cb=5 ?
			basek_c		[]="-TGCA$"				,  ///< las 4 bases en el orden complementario de Kaderali.

			n_basek		=sizeof(basek)-1		,

			nu2dgba		[]="-GCSTKYBARMVWDHN"	,  ///< "convierte" el "codigo numerico" en letra del deg cod; lo contrario de db2nu[]. de "numero a 
			nu2c_dgba   []="-CGSAMRVTYKBWHDN"	,  ///< devuelve el complementario. Comentar aqui codigo binario.
			n_dgba		  =sizeof(nu2dgba)-1	;  ///< '\020 es 16 ?????? 			//n_dgba		=16 ; 


extern Base is_base		[UCHAR_MAX],		///< <> 0  si base. =base, pero para U, =T 
			is_degbase	[UCHAR_MAX],		///< <> 0  si letra valida (cualquiera del cod deg, mayuscula o minuscula
											///< + insercion '-').=base, pero para U, =T 
			c_degbase	[UCHAR_MAX],		///< devuelve base complementaria, tambien para codigo deg. 
											///< El resto no lo modifica.
			grad_deg	[UCHAR_MAX],		///< grado de degeneracion: 0- no permitido, 1-base, 
											///<						  2-Y,R,K....,     4-N.
			is_GC		[UCHAR_MAX];		///< 1-G or C, 0-lo demas.                    Y que con los deg ??????????????????

extern Code bk2nu		[UCHAR_MAX],		///< codigo corto de Kadelari
			bk2c_nu		[UCHAR_MAX],		///< codigo corto de Kadelari complementario a la base
			bkn2c_nu	[n_basek],			///< codigo corto de Kadelari complementario
			ba2nu		[UCHAR_MAX],		///< codigo corto
			ba2c_nu		[UCHAR_MAX],		///< codigo corto complementario a la base
			ban2c_nu	[n_ba],				///< codigo corto complementario
			dbn2c_nu	[n_dgba],			///< cod numerico complementario a la base deg
			db2c_nu		[UCHAR_MAX],		///< cod numerico complementario a la base deg
			db2nu		[UCHAR_MAX];		///< dada la letra devuelve el cod numerico; lo contrario 
											///< de degcod[]; 0- no validos, pero tambien el '-'.
											///< si las bases a y b pueden match se comprueba: 
											///< if (base2bin[a] & base2bin[b]) ...   : ejemplo - si match 
											///< R y G, R y K, ...
											///< (base2bin[a] & base2bin[b]) da el cod de la min letra comun:
											///< R y G = G, R y K = G
											///< para calcular consenso: degcod[base2bin[a] | base2bin[b])]
extern Base dg2ba		[n_dgba][n_ba],		///< caracter a generar,para generar todas las variantes de una base deg
			dg2ban		[n_dgba][n_ba],		///< array de arrays (para cada base deg), 
			dg2bkn		[n_dgba][n_ba];		///< con el que generar no deg variantes

extern const Base	*c_base,*c_basek;
											

class CInit_Cod_Deg
{    public:
		CInit_Cod_Deg()	;
} ;

extern	CInit_Cod_Deg Init_Cod_Deg;
long	CountDegBases  ( const char *sec);
Base    *Generate_DegSec( const char *sec, bool rev, bool complem, long l=0) ;// , long l=0) ;
inline char *Generate_DegSec_char( const char *sec, bool rev, bool complem, long l=0){return (char *)Generate_DegSec(sec, rev,  complem,  l);};

}
//class CBase	// inventar algo para lograr que las converciones sean automaticas: Base - char
//{public:
//	Base b;
//	CBase (char c): b( /*(Base)*/c ){}
//	Base& operator Base(){return b;}
//	char& operator char(){return b;}//(char &)
//};


//namespace Cod_Deg
//{	enum BaseCod {A='A', C='C', T='T', G='G', a='A', C='C', T='T', G='G'};
//};
//extern Base *Generate_DegSec( const char *sec, bool rev, bool compl, long l=0) ; //, long l=0) ;
#endif


