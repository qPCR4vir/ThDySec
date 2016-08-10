/**
* Copyright (C) 2009-2016, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDySec\th_dy_param.cpp
*
* @brief An implementation of the Nearest Neighbor Model Parameters
*
*        Attempt to be a simplification of what
*        developed and reported Santa Lucia: http://www.annualreviews.org/doi/abs/10.1146/annurev.biophys.32.110601.141800
*        This representation is based on the ideas and code of Kaderali (http://bioinformatics.oxfordjournals.org/content/21/10/2375.abstract)
*        but with many modifications, so that the original authors have no responsability on the many erros,
*        simplifications or inconsistencies I have introduced (most files and class names were changed to avoid confusion with originals).
*
*        The original source file had the following header:
*
* //=============================================================================
* // Module:        nnparams.cpp
* // Project:       Diploma Thesis - Probe Selection for DNA Microarrays
* // Type:          implementation - Nearest Neighbor Model / Parameters
* // Language:      c++
* // Compiler:      microsoft visual c++ 6.0, unix/linux gcc
* // System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
* // Database:      none
* // Description:   class CNNParams - Nearest Neighbor Model / Parameters
* // Author:        kaderali
* // Date:          12/2000
* // Copyright:     (c) L. Kaderali, 9/2000 - 12/2000
* //
* // Revision History
* // $              00sep07 LK : created
* //                00dec17 LK : modified to use new nn parameters
* //                00dec29 LK : modified to include dangling end parameters
* //                01jan09 LK : included CalcSelfTM
* //                01feb07 LK : optimized
* // #$
* //=============================================================================
*
* Which is accesible under GNU GPL at: http://dna.engr.uconn.edu/?page_id=85
*
*/

#ifdef WINDOWS_FORM_GUI
#include "stdafx.h"
#pragma unmanaged
#endif

#include <memory.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include "ThDySec/th_dy_param.h"

using namespace std;
const     EnergyRang        G_def (-5,-1) ;    //  G_def , 
const     TemperatureRang    Tm_def (57,63);    //  Tm_def ;
const     SecPosRang        L_def  (20,35);    //  L_def  ;




void    CSaltCorrNN::InitSaltNNMatriz()
{    
    //Copy_oridS(_dS) ;    //Copy_oridH(_dH) ;

    switch (_SaltCorr) 
    {    case StLucia:
            InitStLuciaSaltNNMatriz();
            break ;
        case Owczarzy:
            //Copy_oridS(_dS) ;//InitOwczarzySaltNNMatriz();
            break ;
    }
}

void    CSaltCorrNN::InitStLuciaSaltNNMatriz()
{    Copy_oridS(_dS) ;    

    Base a_1, a, b_1, b ;
    //for (a_1 =0; a_1<=5; a_1++)            // alternativa a memcopy()
    //    for (a=0; a<=5; a++)            
    //        for (b_1=0; b_1<=5; b_1++)    
    //            for (b=0; b<=5; b++)
    //                _dS    [a_1][a] [b_1][b] = GetOriEntr(a_1, a, b_1, b) ;
                



    float correction = 0.5f* 0.368f * float(log (_ConcSalt)) ;
//    Base a_1, a, b_1, b ;

    for (a_1 =0; a_1<=4; a_1++)            // no $-extremo en a_1
    for (a     =1; a  <=4; a++)            // solo bases en a
    for (b_1 =0; b_1<=5; b_1++)    // todos en b_1
    for (b     =0; b  <=5; b++)    // todos en b
        _dS [a_1][a][b_1][b] += correction;         
    
    ;

    for (b   =0; b  <=4;   b++)                // no $-extremo en b
    for (b_1 =1; b_1<=4; b_1++)        // solo bases en b_1
    for (a_1 =0; a_1<=5; a_1++)    // todos en a_1
    for (a   =0; a  <=5;   a++)    // todos en a
        _dS [a_1][a] [b_1][b] += correction;         
}
    //  Entropy dS for given NN pair from parameter table
    //  (    a_1 , a,         b_1 , b ) :   Pairs to look up in form  
    //   5'- a_1 - a -3'/ 3'- b_1 - b -5'

void    CSaltCorrNN::InitOwczarzySaltNNMatriz() // recalcular solo cuando se "pone" nueva GCp
{    float LogSC = log(_ConcSalt);
    Base a_1, a, b_1, b ;
    Copy_oridS(_dS) ;    
    for (a_1 =0; a_1<=5; a_1++)            
    for (a   =0; a  <=5; a  ++)            
    for (b_1 =0; b_1<=5; b_1++)    
    for (b   =0; b  <=5; b  ++)    
        _dS [a_1][a] [b_1][b] += _oridH [a_1][a] [b_1][b]
                                            * ((4.29f * _GCp -3.95f )*(1e-5f)*LogSC+ (9.4e-6f)*LogSC*LogSC);
}    

void    CSaltCorrNN::SetStLuciaSaltCorr( float C1, float C2, float Csalt)  ///  \todo
{    //kplus = Csalt ;    
//    kfac  = 0.368 * log (kplus) ;
    ;
}

float    CSaltCorrNN::UpdateGC(float GC1, long l1,float GC2, long l2)
{    if(l1<l2) // cambiar para que acepte CSec ???, o solo que acepte el GCp ya calculado 
        _GCp = GC1;     // ver si % o solo fraccion -> /100
    else 
        if (l2<l1) 
            _GCp = GC2; 
        else 
            _GCp = ( GC1*l1 + GC2*l2 )/ (l1+l2);
    return  _GCp ;
}

bool    CSaltCorrNN::SetConc(float C1,  float C2, float CationConc)
{    // float LogSC_alt= log(_ConcSalt) ;

    _ConcSalt=CationConc ;
    //float LogSC    = log(_ConcSalt) ;

    if ( ChangeConc(C1, C2) )
    {    InitSaltNNMatriz() ;
        return true ;
    } 
    return false ;
//    if ( IsEq( LogSC,LogSC_alt ) ) // ver si aporta mucho recalcular solo la parte que cambio con COriNN::UpdateMatriz_f_e_el()

}

bool    COriNN::ChangeConc(float C1,  float C2  ) // hubo que update ??
{    _RlogC    = R * log( (C1>C2)?C1-C2/2:C2-C1/2  ) ;

    if ( IsEq (_RlogC, forbidden_entropy)  )
        return false ;

    forbidden_entropy    =_RlogC        ;
    // recalcular los elementos de la matriz original que dependen de 'forbidden_entropy'
    UpdatedSMatriz_forb_entr_elem();
    return true ;
}


bool COriNN::LoadNNParam(istream &isTDP)  
{    char sep[]=";";   // usar el sep global !!!!
    string separ ;
    // usar mejor: istream&  ignore ( streamsize n = 1, int delim = EOF );
    getline  ( isTDP , separ ) ;    

    if ( string::npos == separ.find("5\'XX3\'/3\'YY5")   )     return false; // Error

    Base a_1, a, b_1, b ;
            
    while ( isTDP.good() )
    {    isTDP    >> a_1 ;    isTDP.ignore(1000,sep[0])  ;
        isTDP    >> a ;        isTDP.ignore(1000,sep[0])   ;
        isTDP    >> b_1;        isTDP.ignore(1000,sep[0])    ;
        isTDP    >> b;        isTDP.ignore(1000,sep[0])  ;  isTDP.ignore(1000,sep[0])  ;

        isTDP    >> _oridH [ bk2nu[a_1]  ][ bk2nu[a]  ] [ bk2nu[b_1]  ][bk2nu[b] ]   ;  isTDP. ignore(1000,sep[0]) ;
        isTDP    >> _oridS [ bk2nu[a_1]  ][ bk2nu[a]  ] [ bk2nu[b_1]  ][bk2nu[b] ]  ;
        getline (isTDP, separ) ;    
    }
    return true;

}

bool CSaltCorrNN::LoadNNParam(istream &isTDP)  
{    
    if ( ! COriNN::LoadNNParam(isTDP) ) return false ;
    InitSaltNNMatriz() ;
    return true ;

}

ostream &operator<<(ostream &stream, const CSaltCorrNN &sp) 
{    char sep[]=";" ;
    stream    
            <<"X1" << sep  <<"X2" << sep  <<"Y2" << sep  <<"Y1" << sep  
            <<"5\'XX3\'/3\'YY5" << sep 
            << "OridH"        << sep
            << "OriddS"        << sep 
            << "StLuc.dS"    << sep
            << "StLuc.dH"    << '\n'  ;

    Base a_1, a, b_1, b ;
            
    for (a_1 =0; a_1<=5; a_1++)            
    for (a     =0; a  <=5;   a++)            
    for (b_1 =0; b_1<=5; b_1++)    
    for (b     =0; b  <=5;   b++)
                {    stream    << basek[a_1] << sep  <<basek[a] << sep  <<basek[b_1] << sep  <<basek[b] << sep  
                            << basek[a_1]<< basek[a] <<'/'        << basek[b_1]<< basek[b]<< sep ;
                    stream    << sp._oridH [a_1][a] [b_1][b]        << sep 
                            << sp.GetOriEntr(a_1, a, b_1, b)    << sep ;
                    stream    << sp._dS    [a_1][a] [b_1][b]        << sep 
                            << sp._oridH [a_1][a] [b_1][b]        << '\n' ;
                }
    return stream;

}


void    COriNN::UpdatedSMatriz_forb_entr_elem() // despues de esto hay que update cualquier SaltCorr o al menos una parte
{        
    Base x, y;// , a, b;                            //        basek[]=".ACGT$"    .-0, A-1, C-2, G-3, T-4, $-5 , g=0 /* gap */ , e=5 /* extremo, end */
        /// \todo  instead of $...$  try to use <....>
        // Set all parameters to zero!?     why 0 and not forbidden ?!                    //  erased by Ariel         
        //memset(_oridH,0,sizeof(_oridH));    //memset(_oridS,0,sizeof(_oridS));

      const Base A = bk2nu['A'],            /// \todo adjust!?  Revisar esto si cambia el cod_deg  !!!!!!!!!!!!!!!!!!!!!!
                 C = bk2nu['C'], 
                 T = bk2nu['T'], 
                 G = bk2nu['G'], 
                 g = bk2nu['-'], /* gap */
                 e = bk2nu['$'], /* extremo, end */
               all = n_basek,
             first = 0, last       = all-1,
        first_base = 1,    last_base  = 4 ;
    
    for     ( x=first; x<=last ;x++)        //  bases plus   - and $                                                    // Ariel 
    {    for ( y=first; y<=last ;y++)        //  forbid $./XY etc.
        {    ndS(e,g,x,y)=forbidden_entropy;   //   $-/XY                        // Ariel 
        //    ndS(g,e,x,y)=forbidden_entropy;   //   -$/XY                        // Ariel - caso especial: comienzo de sec
        //    ndS(x,y,g,e)=forbidden_entropy;   //   XY/-$                        // Ariel - caso especial: comienzo de sec
            ndS(x,y,e,g)=forbidden_entropy;   //   XY/$-                        // Ariel
        }
    }
    
    for (      x=first_base; x<=last_base; x++)
    {    for (  y=first_base; y<=last_base; y++)            // Set all X-/Y-, -X/Y- and X-/-Y so, that TM will be VERY small! 
        {    ndS(g,x,y,g)=forbidden_entropy;
            ndS(x,g,g,y)=forbidden_entropy;
            ndS(x,g,y,g)=forbidden_entropy;
            ndS(g,x,g,y)=forbidden_entropy;   //   -X/-Y                        // Ariel
                                        // forbid X-/Y$ and X$/Y- etc., i.e. terminal must not be paired with gap!
            ndS(x,e,y,g)=forbidden_entropy;
            ndS(x,g,y,e)=forbidden_entropy;
            ndS(e,x,g,y)=forbidden_entropy;
            ndS(g,x,e,y)=forbidden_entropy;
                                        // forbid X$/-Y etc.
            ndS(x,e,g,y)=forbidden_entropy;
            ndS(x,g,e,y)=forbidden_entropy;
            ndS(e,x,y,g)=forbidden_entropy;
            ndS(g,x,y,e)=forbidden_entropy;
        }                                // also, forbid x-/-- and --/x-, i.e. no two inner gaps paired
        ndS(x,g,g,g)=forbidden_entropy;
        ndS(g,g,x,g)=forbidden_entropy;
        ndS(x,g,g,e)=forbidden_entropy;  // x-/-$
        ndS(e,g,g,x)=forbidden_entropy;
        ndS(x,g,g,e)=forbidden_entropy;
        ndS(g,x,e,g)=forbidden_entropy;
        ndS(g,e,g,x)=forbidden_entropy;   //   -$/-X                        // Ariel
        ndS(x,g,e,g)=forbidden_entropy;   //   X-/$-     ----- se puede quitar                        // Ariel
    }
    ndS(g,g,g,g)=forbidden_entropy;    // forbid --/--
    ndS(e,g,g,g)=forbidden_entropy;
    ndS(g,g,e,g)=forbidden_entropy;
    ndS(g,e,e,g)=forbidden_entropy;
}// ver cuales de estos se tienen que corregir con las sales.



void    COriNN::InitOriNNMatriz()
{    Base x,y,a,b;                                    //        basek[]=".ACGT$"    .-0, A-1, C-2, G-3, T-4, $-5 , g=0 /* gap */ , e=5 /* extremo, end */
        //  En lugar de $...$  tratar de poner <....>
    // Set all parameters to zero!     y porque 0 y no forbidden ?!                    //  borrado por Ariel         
    //memset(_oridH,0,sizeof(_oridH));    //memset(_oridS,0,sizeof(_oridS));

      const Base A = bk2nu['A'],            //  Revisar esto si cambia el cod_deg  !!!!!!!!!!!!!!!!!!!!!!
                 C = bk2nu['C'], 
                 T = bk2nu['T'], 
                 G = bk2nu['G'], 
                 g = bk2nu['-'], /* gap */
                 e = bk2nu['$'], /* extremo, end */
               all = n_basek,
             first = 0, last       = all-1,
        first_base = 1,    last_base  = 4 ;


    for (x = first; x < all; x++)            // Prohibir todo lo no permitido                                                // Ariel 
    for (y = first; y < all; y++)                                                // Ariel 
    for (a = first; a < all; a++)                                                // Ariel 
    for (b = first; b < all; b++)                                                // Ariel 
    {    ndH(x,y,a,b)= 0                    ;    ndS(x,y,a,b) = 0                  ; }    //   AB/XY                    // Ariel
                //{    ndH(x,y,a,b)=forbidden_enthalpy;    ndS(x,y,a,b)=forbidden_entropy; }    //   AB/XY                    // Ariel

    // y de paso verificar todos estos datos. Dar posibilidad de ajustar solo algunos parametros (correcciones)
    //        S        forbidden_entropy    (_RlogC    ),        // OJO !! dependencia de parte de la matriz dS de las conc ADN
    //        H        forbidden_enthalpy    ( 1e18f    ),        // initialize parameter table! MUY GRANDE


    for     (x = first;x < all; x++)                //  bases plus   - and $                                                    // Ariel 
    {    for (y = first;y < all; y++)                //  forbid $./XY etc.   ? <./XY  or  >./XY  ?
        {    ndH(e,g,x,y)=forbidden_enthalpy;    ndS(e,g,x,y)=forbidden_entropy;   //   $-/XY                        // Ariel 
        //    ndH(g,e,x,y)=forbidden_enthalpy;    ndS(g,e,x,y)=forbidden_entropy;   //   -$/XY                        // Ariel - caso especial: comienzo de sec
        //    ndH(x,y,g,e)=forbidden_enthalpy;    ndS(x,y,g,e)=forbidden_entropy;   //   XY/-$                        // Ariel - caso especial: comienzo de sec
            ndH(x,y,e,g)=forbidden_enthalpy;    ndS(x,y,e,g)=forbidden_entropy;   //   XY/$-                        // Ariel
        }
    }

    for (x=first_base;x<=last_base;x++)                    // solo bases, no - or $
    {    for (y=first_base;y<=last_base;y++)                // Set all X-/Y-, -X/Y- and X-/-Y so, that TM will be VERY small! ??? H/S ??
        {    ndH(g,x,y,g)=forbidden_enthalpy;    ndS(g,x,y,g)=forbidden_entropy;   //   -X/Y- 
            ndH(x,g,g,y)=forbidden_enthalpy;    ndS(x,g,g,y)=forbidden_entropy;   //   X-/-Y 
            ndH(x,g,y,g)=forbidden_enthalpy;    ndS(x,g,y,g)=forbidden_entropy;   //   X-/Y-
            ndH(g,x,g,y)=forbidden_enthalpy;    ndS(g,x,g,y)=forbidden_entropy;   //   -X/-Y                        // Ariel
                                        // forbid X-/Y$ and X$/Y- etc., i.e. terminal must not be paired with gap!
            ndH(x,e,y,g)=forbidden_enthalpy;    ndS(x,e,y,g)=forbidden_entropy;   //   X$/Y$ 
            ndH(x,g,y,e)=forbidden_enthalpy;    ndS(x,g,y,e)=forbidden_entropy;   //   X-/Y$ 
            ndH(e,x,g,y)=forbidden_enthalpy;    ndS(e,x,g,y)=forbidden_entropy;   //   $X/-Y
            ndH(g,x,e,y)=forbidden_enthalpy;    ndS(g,x,e,y)=forbidden_entropy;   //   -X/$Y
                                        // forbid X$/-Y etc.
            ndH(x,e,g,y)=forbidden_enthalpy;    ndS(x,e,g,y)=forbidden_entropy;   //   X$/-Y
            ndH(x,g,e,y)=forbidden_enthalpy;    ndS(x,g,e,y)=forbidden_entropy;   //   X-/$Y
            ndH(e,x,y,g)=forbidden_enthalpy;    ndS(e,x,y,g)=forbidden_entropy;   //   X$/Y-
            ndH(g,x,y,e)=forbidden_enthalpy;    ndS(g,x,y,e)=forbidden_entropy;   //   -X/Y$
        }                                // also, forbid x-/-- and --/x-, i.e. no two inner gaps paired
        ndH(x,g,g,g)=forbidden_enthalpy;        ndS(x,g,g,g)=forbidden_entropy;   //   -X/--
        ndH(g,g,x,g)=forbidden_enthalpy;        ndS(g,g,x,g)=forbidden_entropy;   //   --/X-
        ndH(x,g,g,e)=forbidden_enthalpy;        ndS(x,g,g,e)=forbidden_entropy;   //   X-/-$
        ndH(e,g,g,x)=forbidden_enthalpy;        ndS(e,g,g,x)=forbidden_entropy;   //   $-/-X     ----- se puede quitar
        ndH(g,e,x,g)=forbidden_enthalpy;        ndS(x,g,g,e)=forbidden_entropy;   //   -$/X-    
        ndH(g,x,e,g)=forbidden_enthalpy;        ndS(g,x,e,g)=forbidden_entropy;   //   -X/$-     ----- se puede quitar    
        ndH(g,e,g,x)=forbidden_enthalpy;        ndS(g,e,g,x)=forbidden_entropy;   //   -$/-X                        // Ariel
        ndH(x,g,e,g)=forbidden_enthalpy;        ndS(x,g,e,g)=forbidden_entropy;   //   X-/$-                     // Ariel    ----- se puede quitar    
    }
    ndH(g,g,g,g)=forbidden_enthalpy;        ndS(g,g,g,g)=forbidden_entropy;    // forbid   --/--
    ndH(e,g,g,g)=forbidden_enthalpy;        ndS(e,g,g,g)=forbidden_entropy;    // forbid   $-/--     ----- se puede quitar
    ndH(g,g,e,g)=forbidden_enthalpy;        ndS(g,g,e,g)=forbidden_entropy;    // forbid   --/$-     ----- se puede quitar
    ndH(g,e,e,g)=forbidden_enthalpy;        ndS(g,e,e,g)=forbidden_entropy;    // forbid   -$/$-     ----- se puede quitar

    for             (x=first_base; x<=last_base; x++)            // Interior loops (double Mismatches)    iloop_entropy(-0.97f) ?*1000?,  iloop_enthalpy        ( 0.00f    ),
        for         (y=first_base; y<=last_base; y++)
            for     (a=first_base; a<=last_base; a++)
                for (b=first_base; b<=last_base; b++)
                    // AT and CG pair, and as A=1, C=2, G=3, T=4 this means
                    // we have Watson-Crick pairs if (x+a==5) and (y+b)==5.
                    if (!((x+a==5)||(y+b==5)))                        // innecesario !!!?????????????????
                    {    ndH(x,y,a,b) = iloop_enthalpy;// No watson-crick-pair, i.e. double mismatch!
                        ndS(x,y,a,b) = iloop_entropy;// set enthalpy/entropy to loop expansion!
                    }
    for     (x=first_base; x<=last_base; x++)            // xy/-- and --/xy (Bulge Loops of size > 1)
        for (y=first_base; y<=last_base; y++)
        {    ndH(x,y,g,g) = bloop_enthalpy;            ndS(x,y,g,g) = bloop_entropy;    //        bloop_entropy    (-1.30f    ),    // xy/-- and --/xy (Bulge Loops of size > 1)
            ndH(g,g,x,y) = bloop_enthalpy;            ndS(g,g,x,y) = bloop_entropy;    //        bloop_enthalpy    ( 0.00f    ),
        }
    // x-/ya and xa/y- as well as -x/ay and ax/-y
    // bulge opening and closing parameters with
    // adjacent matches / mismatches
    // obulge_mism and cbulge_mism chosen so high to avoid 
    //     AAAAAAAAA
    //     T--G----T
    // being better than
    //     AAAAAAAAA
    //     TG------T
    for         (x=first_base; x<=last_base; x++)
        for     (y=first_base; y<=last_base; y++)
            for (a=first_base; a<=last_base; a++)
            {    if (x+y==5)                            // other base pair matches!         :                     H de -2660 a -2660,      S de -14.22 a -14.22       
                {    ndH(x,g,y,a)=obulge_match_H;  ndS(x,g,y,a)=obulge_match_S;// bulge opening    obulge_match_H(-2.66f * 1000),    obulge_match_S(-14.22f),
                    ndH(x,a,y,g)=obulge_match_H;  ndS(x,a,y,g)=obulge_match_S;
                    ndH(g,x,a,y)=cbulge_match_H;  ndS(g,x,a,y)=cbulge_match_S;// bulge closing    cbulge_match_H(-2.66f * 1000),    cbulge_match_S(-14.22f),
                    ndH(a,x,g,y)=cbulge_match_H;  ndS(a,x,g,y)=cbulge_match_S;
                }    else                            // mismatch in other base pair!
                {    ndH(x,g,y,a)=obulge_mism_H;   ndS(x,g,y,a)=obulge_mism_S;// bulge opening
                    ndH(x,a,y,g)=obulge_mism_H;   ndS(x,a,y,g)=obulge_mism_S;
                    ndH(g,x,a,y)=cbulge_mism_H;   ndS(g,x,a,y)=cbulge_mism_S;// bulge closing
                    ndH(a,x,g,y)=cbulge_mism_H;   ndS(a,x,g,y)=cbulge_mism_S;
                }
            }
    // Watson-Crick pairs (note that only ten are unique, as obviously  :        H de -10600 a -7200,           S de -24.4 a -19.9   (PM)
    // 5'-AG-3'/3'-TC-5'  =  5'-CT-3'/3'-GA-5' etc.
    ndH(A,A,T,T)=-7.6f*1000;  ndS(A,A,T,T)=-21.3f;   // AA/TT 04
    ndH(A,C,T,G)=-8.4f*1000;  ndS(A,C,T,G)=-22.4f;   // AC/TG adapted GT/CA
    ndH(A,G,T,C)=-7.8f*1000;  ndS(A,G,T,C)=-21.0f;   // AG/TC adapted CT/GA
    ndH(A,T,T,A)=-7.2f*1000;  ndS(A,T,T,A)=-20.4f;   // AT/TA 04 //    basek[]=".ACGT$" :.-0, A-1, C-2, G-3, T-4, $-5 , g=0 /* gap */ , e=5 /* extremo, end */
    ndH(C,A,G,T)=-8.5f*1000;  ndS(C,A,G,T)=-22.7f;   // CA/GT 04
    ndH(C,C,G,G)=-8.0f*1000;  ndS(C,C,G,G)=-19.9f;   // CC/GG adapted GG/CC
    ndH(C,G,G,C)=-10.6f*1000; ndS(C,G,G,C)=-27.2f;   // CG/GC 04
    ndH(C,T,G,A)=-7.8f*1000;  ndS(C,T,G,A)=-21.0f;   // CT/GA 04
    ndH(G,A,C,T)=-8.2f*1000;  ndS(G,A,C,T)=-22.2f;   // GA/CT 04
    ndH(G,C,C,G)=-9.8f*1000;  ndS(G,C,C,G)=-24.4f;   // GC/CG 04
    ndH(G,G,C,C)=-8.0f*1000;  ndS(G,G,C,C)=-19.9f;   // GG/CC 04
    ndH(G,T,C,A)=-8.4f*1000;  ndS(G,T,C,A)=-22.4f;   // GT/CA 04
    ndH(T,A,A,T)=-7.2f*1000;  ndS(T,A,A,T)=-21.3f;   // TA/AT 04
    ndH(T,C,A,G)=-8.2f*1000;  ndS(T,C,A,G)=-22.2f;   // TC/AG adapted GA/CT
    ndH(T,G,A,C)=-8.5f*1000;  ndS(T,G,A,C)=-22.7f;   // TG/AC adapted CA/GT
    ndH(T,T,A,A)=-7.6f*1000;  ndS(T,T,A,A)=-21.3f;   // TT/AA adapted AA/TT

    // A-C Mismatches (Values for pH 7.0)                                :        H de -700 a +7600,           S de -3.8 a 20.2     (MM)
    ndH(A,A,C,T)=7.6f*1000;   ndS(A,A,C,T)=20.2f;    // AA/CT
    ndH(A,A,T,C)=2.3f*1000;   ndS(A,A,T,C)=4.6f;     // AA/TC
    ndH(A,C,C,G)=-0.7f*1000;  ndS(A,C,C,G)=-3.8f;    // AC/CG
    ndH(A,C,T,A)=5.3f*1000;   ndS(A,C,T,A)=14.6f;    // AC/TA
    ndH(A,G,C,C)=0.6f*1000;   ndS(A,G,C,C)=-0.6f;    // AG/CC
    ndH(A,T,C,A)=5.3f*1000;   ndS(A,T,C,A)=14.6f;    // AT/CA
    ndH(C,A,A,T)=3.4f*1000;   ndS(C,A,A,T)=8.0f;     // CA/AT
    ndH(C,A,G,C)=1.9f*1000;   ndS(C,A,G,C)=3.7f;     // CA/GC
    ndH(C,C,A,G)=5.2f*1000;   ndS(C,C,A,G)=14.2f;    // CC/AG
    ndH(C,C,G,A)=0.6f*1000;   ndS(C,C,G,A)=-0.6f;    // CC/GA
    ndH(C,G,A,C)=1.9f*1000;   ndS(C,G,A,C)=3.7f;     // CG/AC
    ndH(C,T,A,A)=2.3f*1000;   ndS(C,T,A,A)=4.6f;     // CT/AA
    ndH(G,A,C,C)=5.2f*1000;   ndS(G,A,C,C)=14.2f;    // GA/CC
    ndH(G,C,C,A)=-0.7f*1000;  ndS(G,C,C,A)=-3.8f;    // GC/CA
    ndH(T,A,A,C)=3.4f*1000;   ndS(T,A,A,C)=8.0f;     // TA/AC
    ndH(T,C,A,A)=7.6f*1000;   ndS(T,C,A,A)=20.2f;    // TC/AA

    // C-T Mismatches                                                    :        H de -1500 a +5200,          S de -6.2 a 13.5     (MM)
    ndH(A,C,T,T)=0.7f*1000;   ndS(A,C,T,T)=0.2f;     // AC/TT
    ndH(A,T,T,C)=-1.2f*1000;  ndS(A,T,T,C)=-6.2f;    // AT/TC
    ndH(C,A,T,T)=1.0f*1000;   ndS(C,A,T,T)=0.7f;     // CA/TT
    ndH(C,C,G,T)=-0.8f*1000;  ndS(C,C,G,T)=-4.5f;    // CC/GT
    ndH(C,C,T,G)=5.2f*1000;   ndS(C,C,T,G)=13.5f;    // CC/TG
    ndH(C,G,T,C)=-1.5f*1000;  ndS(C,G,T,C)=-6.1f;    // CG/TC
    ndH(C,T,G,C)=-1.5f*1000;  ndS(C,T,G,C)=-6.1f;    // CT/GC
    ndH(C,T,T,A)=-1.2f*1000;  ndS(C,T,T,A)=-6.2f;    // CT/TA
    ndH(G,C,C,T)=2.3f*1000;   ndS(G,C,C,T)=5.4f;     // GC/CT
    ndH(G,T,C,C)=5.2f*1000;   ndS(G,T,C,C)=13.5f;    // GT/CC
    ndH(T,A,C,T)=1.2f*1000;   ndS(T,A,C,T)=0.7f;     // TA/CT
    ndH(T,C,C,G)=2.3f*1000;   ndS(T,C,C,G)=5.4f;     // TC/CG
    ndH(T,C,A,T)=1.2f*1000;   ndS(T,C,A,T)=0.7f;     // TC/AT
    ndH(T,G,C,C)=-0.8f*1000;  ndS(T,G,C,C)=-4.5f;    // TG/CC
    ndH(T,T,C,A)=0.7f*1000;   ndS(T,T,C,A)=0.2f;     // TT/CA
    ndH(T,T,A,C)=1.0f*1000;   ndS(T,T,A,C)=0.7f;     // TT/AC

    // G-A Mismatches
    ndH(A,A,G,T)=3.0f*1000;   ndS(A,A,G,T)=7.4f;     // AA/GT
    ndH(A,A,T,G)=-0.6f*1000;  ndS(A,A,T,G)=-2.3f;    // AA/TG
    ndH(A,C,G,G)=0.5f*1000;   ndS(A,C,G,G)=3.2f;     // AC/GG
    ndH(A,G,G,C)=-4.0f*1000;  ndS(A,G,G,C)=-13.2f;   // AG/GC
    ndH(A,G,T,A)=-0.7f*1000;  ndS(A,G,T,A)=-2.3f;    // AG/TA
    ndH(A,T,G,A)=-0.7f*1000;  ndS(A,T,G,A)=-2.3f;    // AT/GA
    ndH(C,A,G,G)=-0.7f*1000;  ndS(C,A,G,G)=-2.3f;    // CA/GG
    ndH(C,G,G,A)=-4.0f*1000;  ndS(C,G,G,A)=-13.2f;   // CG/GA
    ndH(G,A,A,T)=0.7f*1000;   ndS(G,A,A,T)=0.7f;     // GA/AT
    ndH(G,A,C,G)=-0.6f*1000;  ndS(G,A,C,G)=-1.0f;    // GA/CG
    ndH(G,C,A,G)=-0.6f*1000;  ndS(G,C,A,G)=-1.0f;    // GC/AG
    ndH(G,G,A,C)=-0.7f*1000;  ndS(G,G,A,C)=-2.3f;    // GG/AC
    ndH(G,G,C,A)=0.5f*1000;   ndS(G,G,C,A)=3.2f;     // GG/CA
    ndH(G,T,A,A)=-0.6f*1000;  ndS(G,T,A,A)=-2.3f;    // GT/AA
    ndH(T,A,A,G)=0.7f*1000;   ndS(T,A,A,G)=0.7f;     // TA/AG
    ndH(T,G,A,A)=3.0f*1000;   ndS(T,G,A,A)=7.4f;     // TG/AA

    // G-T Mismatches
    ndH(A,G,T,T)=1.0f*1000;   ndS(A,G,T,T)=0.9f;     // AG/TT
    ndH(A,T,T,G)=-2.5f*1000;  ndS(A,T,T,G)=-8.3f;    // AT/TG
    ndH(C,G,G,T)=-4.1f*1000;  ndS(C,G,G,T)=-11.7f;   // CG/GT
    ndH(C,T,G,G)=-2.8f*1000;  ndS(C,T,G,G)=-8.0f;    // CT/GG
    ndH(G,A,T,T)=-1.3f*1000;  ndS(G,A,T,T)=-5.3f;    // GA/TT
    ndH(G,C,T,G)=-4.4f*1000;  ndS(G,C,T,G)=-12.3f;   // GC/TG
    ndH(G,G,C,T)=3.3f*1000;   ndS(G,G,C,T)=10.4f;    // GG/CT
    ndH(G,G,T,C)=-2.8f*1000;  ndS(G,G,T,C)=-8.0f;    // GG/TC
//    ndH(G,G,T,T)=5.8f*1000;   ndS(G,G,T,T)=16.3f;    // GG/TT   ???
    ndH(G,T,C,G)=-4.4f*1000;  ndS(G,T,C,G)=-12.3f;   // GT/CG
    ndH(G,T,T,A)=-2.5f*1000;  ndS(G,T,T,A)=-8.3f;    // GT/TA
//    ndH(G,T,T,G)=4.1f*1000;   ndS(G,T,T,G)=9.5f;     // GT/TG   ???
    ndH(T,A,G,T)=-0.1f*1000;  ndS(T,A,G,T)=-1.7f;    // TA/GT
    ndH(T,C,G,G)=3.3f*1000;   ndS(T,C,G,G)=10.4f;    // TC/GG
    ndH(T,G,A,T)=-0.1f*1000;  ndS(T,G,A,T)=-1.7f;    // TG/AT
    ndH(T,G,G,C)=-4.1f*1000;  ndS(T,G,G,C)=-11.7f;   // TG/GC
//    ndH(T,G,G,T)=-1.4f*1000;  ndS(T,G,G,T)=-6.2f;    // TG/GT
    ndH(T,T,A,G)=-1.3f*1000;  ndS(T,T,A,G)=-5.3f;    // TT/AG
    ndH(T,T,G,A)=1.0f*1000;   ndS(T,T,G,A)=0.9f;     // TT/GA
//    ndH(T,T,G,G)=5.8f*1000;   ndS(T,T,G,G)=16.3f;    // TT/GG

    // A-A Mismatches
    ndH(A,A,A,T)=4.7f*1000;   ndS(A,A,A,T)=12.9f;    // AA/AT
    ndH(A,A,T,A)=1.2f*1000;   ndS(A,A,T,A)=1.7f;     // AA/TA
    ndH(A,C,A,G)=-2.9f*1000;  ndS(A,C,A,G)=-9.8f;    // AC/AG
    ndH(A,G,A,C)=-0.9f*1000;  ndS(A,G,A,C)=-4.2f;    // AG/AC
    ndH(A,T,A,A)=1.2f*1000;   ndS(A,T,A,A)=1.7f;     // AT/AA
    ndH(C,A,G,A)=-0.9f*1000;  ndS(C,A,G,A)=-4.2f;    // CA/GA
    ndH(G,A,C,A)=-2.9f*1000;  ndS(G,A,C,A)=-9.8f;    // GA/CA
    ndH(T,A,A,A)=4.7f*1000;   ndS(T,A,A,A)=12.9f;    // TA/AA

    // C-C Mismatches
    ndH(A,C,T,C)=0.0f*1000;   ndS(A,C,T,C)=-4.4f;    // AC/TC
    ndH(C,A,C,T)=6.1f*1000;   ndS(C,A,C,T)=16.4f;    // CA/CT
    ndH(C,C,C,G)=3.6f*1000;   ndS(C,C,C,G)=8.9f;     // CC/CG
    ndH(C,C,G,C)=-1.5f*1000;  ndS(C,C,G,C)=-7.2f;    // CC/GC
    ndH(C,G,C,C)=-1.5f*1000;  ndS(C,G,C,C)=-7.2f;    // CG/CC
    ndH(C,T,C,A)=0.0f*1000;   ndS(C,T,C,A)=-4.4f;    // CT/CA
    ndH(G,C,C,C)=3.6f*1000;   ndS(G,C,C,C)=8.9f;     // GC/CC
    ndH(T,C,A,C)=6.1f*1000;   ndS(T,C,A,C)=16.4f;    // TC/AC

    // G-G Mismatches
    ndH(A,G,T,G)=-3.1f*1000;  ndS(A,G,T,G)=-9.5f;    // AG/TG
    ndH(C,G,G,G)=-4.9f*1000;  ndS(C,G,G,G)=-15.3f;   // CG/GG
    ndH(G,A,G,T)=1.6f*1000;   ndS(G,A,G,T)=3.6f;     // GA/GT
    ndH(G,C,G,G)=-6.0f*1000;  ndS(G,C,G,G)=-15.8f;   // GC/GG
    ndH(G,G,C,G)=-6.0f*1000;  ndS(G,G,C,G)=-15.8f;   // GG/CG
    ndH(G,G,G,C)=-4.9f*1000;  ndS(G,G,G,C)=-15.3f;   // GG/GC
    ndH(G,T,G,A)=-3.1f*1000;  ndS(G,T,G,A)=-9.5f;    // GT/GA
    ndH(T,G,A,G)=1.6f*1000;   ndS(T,G,A,G)=3.6f;     // TG/AG

    // T-T Mismatches
    ndH(A,T,T,T)=-2.7f*1000;  ndS(A,T,T,T)=-10.8f;   // AT/TT
    ndH(C,T,G,T)=-5.0f*1000;  ndS(C,T,G,T)=-15.8f;   // CT/GT
    ndH(G,T,C,T)=-2.2f*1000;  ndS(G,T,C,T)=-8.4f;    // GT/CT
    ndH(T,A,T,T)=0.2f*1000;   ndS(T,A,T,T)=-1.5f;    // TA/TT
    ndH(T,C,T,G)=-2.2f*1000;  ndS(T,C,T,G)=-8.4f;    // TC/TG
    ndH(T,G,T,C)=-5.0f*1000;  ndS(T,G,T,C)=-15.8f;   // TG/TC
    ndH(T,T,A,T)=0.2f*1000;   ndS(T,T,A,T)=-1.5f;    // TT/AT
    ndH(T,T,T,A)=-2.7f*1000;  ndS(T,T,T,A)=-10.8f;   // TT/TA

    // Dangling Ends            //        basek[]=".ACGT$"    .-0, A-1, C-2, G-3, T-4, $-5 , g=0 /* gap */ , e=5 /* extremo, end */
    ndH(e,A,A,T)=-0.7f*1000;  ndS(e,A,A,T)=-0.8f;    // $A/AT
    ndH(e,A,C,T)=4.4f*1000;   ndS(e,A,C,T)=14.9f;    // $A/CT
    ndH(e,A,G,T)=-1.6f*1000;  ndS(e,A,G,T)=-3.6f;    // $A/GT
    ndH(e,A,T,T)=2.9f*1000;   ndS(e,A,T,T)=10.4f;    // $A/TT
    ndH(e,C,A,G)=-2.1f*1000;  ndS(e,C,A,G)=-3.9f;    // $C/AG
    ndH(e,C,C,G)=-0.2f*1000;  ndS(e,C,C,G)=-0.1f;    // $C/CG
    ndH(e,C,G,G)=-3.9f*1000;  ndS(e,C,G,G)=-11.2f;   // $C/GG
    ndH(e,C,T,G)=-4.4f*1000;  ndS(e,C,T,G)=-13.1f;   // $C/TG
    ndH(e,G,A,C)=-5.9f*1000;  ndS(e,G,A,C)=-16.5f;   // $G/AC
    ndH(e,G,C,C)=-2.6f*1000;  ndS(e,G,C,C)=-7.4f;    // $G/CC
    ndH(e,G,G,C)=-3.2f*1000;  ndS(e,G,G,C)=-10.4f;   // $G/GC
    ndH(e,G,T,C)=-5.2f*1000;  ndS(e,G,T,C)=-15.0f;   // $G/TC
    ndH(e,T,A,A)=-0.5f*1000;  ndS(e,T,A,A)=-1.1f;    // $T/AA
    ndH(e,T,C,A)=4.7f*1000;   ndS(e,T,C,A)=14.2f;    // $T/CA
    ndH(e,T,G,A)=-4.1f*1000;  ndS(e,T,G,A)=-13.1f;   // $T/GA
    ndH(e,T,T,A)=-3.8f*1000;  ndS(e,T,T,A)=-12.6f;   // $T/TA
    ndH(A,e,T,A)=-2.9f*1000;  ndS(A,e,T,A)=-7.6f;    // A$/TA
    ndH(A,e,T,C)=-4.1f*1000;  ndS(A,e,T,C)=-13.0f;   // A$/TC
    ndH(A,e,T,G)=-4.2f*1000;  ndS(A,e,T,G)=-15.0f;   // A$/TG
    ndH(A,e,T,T)=-0.2f*1000;  ndS(A,e,T,T)=-0.5f;    // A$/TT
    ndH(A,A,e,T)=0.2f*1000;   ndS(A,A,e,T)=2.3f;     // AA/$T
    ndH(A,A,T,e)=-0.5f*1000;  ndS(A,A,T,e)=-1.1f;    // AA/T$
    ndH(A,C,e,G)=-6.3f*1000;  ndS(A,C,e,G)=-17.1f;   // AC/$G
    ndH(A,C,T,e)=4.7f*1000;   ndS(A,C,T,e)=14.2f;    // AC/T$
    ndH(A,G,e,C)=-3.7f*1000;  ndS(A,G,e,C)=-10.0f;   // AG/$C
    ndH(A,G,T,e)=-4.1f*1000;  ndS(A,G,T,e)=-13.1f;   // AG/T$
    ndH(A,T,e,A)=-2.9f*1000;  ndS(A,T,e,A)=-7.6f;    // AT/$A
    ndH(A,T,T,e)=-3.8f*1000;  ndS(A,T,T,e)=-12.6f;   // AT/T$
    ndH(C,e,G,A)=-3.7f*1000;  ndS(C,e,G,A)=-10.0f;   // C$/GA
    ndH(C,e,G,C)=-4.0f*1000;  ndS(C,e,G,C)=-11.9f;   // C$/GC
    ndH(C,e,G,G)=-3.9f*1000;  ndS(C,e,G,G)=-10.9f;   // C$/GG
    ndH(C,e,G,T)=-4.9f*1000;  ndS(C,e,G,T)=-13.8f;   // C$/GT
    ndH(C,A,e,T)=0.6f*1000;   ndS(C,A,e,T)=3.3f;     // CA/$T
    ndH(C,A,G,e)=-5.9f*1000;  ndS(C,A,G,e)=-16.5f;   // CA/G$
    ndH(C,C,e,G)=-4.4f*1000;  ndS(C,C,e,G)=-12.6f;   // CC/$G
    ndH(C,C,G,e)=-2.6f*1000;  ndS(C,C,G,e)=-7.4f;    // CC/G$
    ndH(C,G,e,C)=-4.0f*1000;  ndS(C,G,e,C)=-11.9f;   // CG/$C
    ndH(C,G,G,e)=-3.2f*1000;  ndS(C,G,G,e)=-10.4f;   // CG/G$
    ndH(C,T,e,A)=-4.1f*1000;  ndS(C,T,e,A)=-13.0f;   // CT/$A
    ndH(C,T,G,e)=-5.2f*1000;  ndS(C,T,G,e)=-15.0f;   // CT/G$
    ndH(G,e,C,A)=-6.3f*1000;  ndS(G,e,C,A)=-17.1f;   // G$/CA
    ndH(G,e,C,C)=-4.4f*1000;  ndS(G,e,C,C)=-12.6f;   // G$/CC
    ndH(G,e,C,G)=-5.1f*1000;  ndS(G,e,C,G)=-14.0f;   // G$/CG
    ndH(G,e,C,T)=-4.0f*1000;  ndS(G,e,C,T)=-10.9f;   // G$/CT
    ndH(G,A,e,T)=-1.1f*1000;  ndS(G,A,e,T)=-1.6f;    // GA/$T
    ndH(G,A,C,e)=-2.1f*1000;  ndS(G,A,C,e)=-3.9f;    // GA/C$
    ndH(G,C,e,G)=-5.1f*1000;  ndS(G,C,e,G)=-14.0f;   // GC/$G
    ndH(G,C,C,e)=-0.2f*1000;  ndS(G,C,C,e)=-0.1f;    // GC/C$
    ndH(G,G,e,C)=-3.9f*1000;  ndS(G,G,e,C)=-10.9f;   // GG/$C
    ndH(G,G,C,e)=-3.9f*1000;  ndS(G,G,C,e)=-11.2f;   // GG/C$
    ndH(G,T,e,A)=-4.2f*1000;  ndS(G,T,e,A)=-15.0f;   // GT/$A
    ndH(G,T,C,e)=-4.4f*1000;  ndS(G,T,C,e)=-13.1f;   // GT/C$
    ndH(T,e,A,A)=0.2f*1000;   ndS(T,e,A,A)=2.3f;     // T$/AA
    ndH(T,e,A,C)=0.6f*1000;   ndS(T,e,A,C)=3.3f;     // T$/AC
    ndH(T,e,A,G)=-1.1f*1000;  ndS(T,e,A,G)=-1.6f;    // T$/AG
    ndH(T,e,A,T)=-6.9f*1000;  ndS(T,e,A,T)=-20.0f;   // T$/AT
    ndH(T,A,e,T)=-6.9f*1000;  ndS(T,A,e,T)=-20.0f;   // TA/$T
    ndH(T,A,A,e)=-0.7f*1000;  ndS(T,A,A,e)=-0.7f;    // TA/A$
    ndH(T,C,e,G)=-4.0f*1000;  ndS(T,C,e,G)=-10.9f;   // TC/$G
    ndH(T,C,A,e)=4.4f*1000;   ndS(T,C,A,e)=14.9f;    // TC/A$
    ndH(T,G,e,C)=-4.9f*1000;  ndS(T,G,e,C)=-13.8f;   // TG/$C
    ndH(T,G,A,e)=-1.6f*1000;  ndS(T,G,A,e)=-3.6f;    // TG/A$
    ndH(T,T,e,A)=-0.2f*1000;  ndS(T,T,e,A)=-0.5f;    // TT/$A
    ndH(T,T,A,e)=2.9f*1000;   ndS(T,T,A,e)=10.4f;    // TT/A$
    return;
}




