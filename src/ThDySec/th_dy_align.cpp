/**
* Copyright (C) 2009-2016, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\src\ThDySec\th_dy_align.cpp
*
* @brief  Thermodynamic Alignment Algorithm
*
* This representation is based on the ideas and code of Kaderali et al. (http://bioinformatics.oxfordjournals.org/content/21/10/2375.abstract)
* but with many modifications, so that the original authors have no responsability on the many erros,
* simplifications or inconsistencies I have introduced (most files and class names were changed to avoid confusion with originals).
* 
* The original source file had the following header:
*
*     //=============================================================================
*     // Module:        thermalign.cpp
*     // Project:       Diploma Thesis - Probe Selection for DNA Microarrays
*     // Type:          implementation - Thermodynamic Alignment.
*     // Language:      c++
*     // Compiler:      microsoft visual c++ 6.0, unix/linux gcc
*     // System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
*     // Database:      none
*     // Description:   class CThermAlign - Thermodynamic Alignment Algorithm
*     // Author:        kaderali
*     // Date:          9/2000 - 01/2001
*     // Copyright:     (c) L. Kaderali, 9/2000 - 01/2001
*     //
*     // Revision History
*     // $ 00sep04 LK : created
*     // $ 00dec30 LK : changed to do local alignment of probe against
*     //                one entire sequence
*     // $ 01feb07 LK : optimized!
*     // $ 01aug06 LK : corrected, included salt and concentration input;
*     // $              made true local alignment (problem with initial mismatches!)
*     // #$
*     //=============================================================================
*
* Which is accesible under GNU GPL at: http://dna.engr.uconn.edu/?page_id=85
*
*/

#ifdef WINDOWS_FORM_GUI
#include "stdafx.h"
#pragma unmanaged
#endif

#include "ThDySec/th_dy_align.h"
#include <set>
using namespace std;

char sep[]=";" ;
//            Step: st_ :      0,1,2,3,4,5,6,7,8
const int ThDyAlign::sti []={0,2,2,2,2,1,1,0};
const int ThDyAlign::stj []={0,0,1,0,2,2,2,2};
const int ThDyAlign::sti1[]={0,1,1,1,1,1,0,0};
const int ThDyAlign::stj1[]={0,0,0,1,1,1,1,1};

    //enum Step {    st_a=1, st_i=1, st_b=2, st_j=2, st_d=3, NoStep=0, 
    //            st_0=0, st_1, st_2, st_3, st_4, st_5, st_6, st_7,
    //            c_01=01, c_02=2,c_04=4, c_05=5, c_06=6, c_08=8, c_09=9, c_10=10, 
    //            
    //          } *_pre/*,*_pre0, *_pre1, *_pre2*/ ;   // back pointer matrix a:i+1 dirc sonda, b:j+1 dir tg, d:diag

ThDyAlign::ThDyAlign(LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar, float InitMax) // ------------------------ ThDyAlign    --------
            :_NNpar        (NNpar) ,                                //_sd(nullptr),
             _dH0(nullptr),_dH1(nullptr),_dH2(nullptr),_dS0(nullptr),_dS1(nullptr),_dS2(nullptr),_pre(nullptr),
             _InitMax    (InitMax),
             _LenSond(MaxLenSond+2),  _LenTarg(MaxLenTarg+2), 
             _Ta(CtoK(60)),                                        // _Ta(NNpar->kein_Tm), --asi estaba. Se justifica?
             _Ta_lastCalc(-1),                                    // todavia no se ha calculado. Aqui estara la Ta del ultimo calculo
             _Tm_sig(CtoK(15))                                    //_Tm_min(CtoK(57)), _Tm_max(CtoK(63)), _TableSize ( 0 )
{    
    force_ResizeTable();        
}

void ThDyAlign::ResizeTable(LonSecPos LenSond, LonSecPos LenTarg)        //    --------------------------------------------------------------    ResizeTable    ---------
{    
    bool initBorder = (LenSond+2 > _LenSond) ;    // Si la nueva sonda es mayor que las anteriores (tabla mas "ancha") ajustar los bordes
    _LenSond     = LenSond +2 ;                    // Para que incluya los ´$´    ?
    _LenTarg     = LenTarg +2 ;                    //     inline float &dS0(long i, long j)const{return _dS0[i +j*_LenSondPlus1]; }

    if (_TableSize >= (_LenSond+1)*(_LenTarg+1) )     // la tabla actual alcanza (debe ser la regla). No cambio tamano. Inicializo solo si es mas ancha que antes
    {    if (  initBorder)             InitBorder();        
        return ;
    }
    force_ResizeTable() ;
}

void ThDyAlign::force_ResizeTable()                        //    --------------------------------------------------------------    force_ResizeTable    ---------
{    delete[] _dH0;    delete[] _dH1;    delete[] _dH2;
    delete[] _dS0;    delete[] _dS1;    delete[] _dS2;    
    delete[] _pre;                                    // seguro ????????
    _dH0 = _dH1 = _dH2 = _dS0 = _dS1 = _dS2 = nullptr;    _pre = nullptr;

    _TableSize = (_LenSond+1)*(_LenTarg+1);
    try {
            _dH0 = new Energy  [_TableSize];    _dH1 = new Energy  [_TableSize];    _dH2 = new Energy  [_TableSize];
            _dS0 = new Entropy [_TableSize];    _dS1 = new Entropy [_TableSize];    _dS2 = new Entropy [_TableSize];
            _pre = new Step    [_TableSize];
        }
    catch (...){    delete[] _dH0;        delete[] _dH1;        delete[] _dH2;
                    delete[] _dS0;        delete[] _dS1;        delete[] _dS2;    
                    delete[] _pre; 
                    _dH0 = _dH1 = _dH2 = _dS0 = _dS1 = _dS2 = nullptr;    _pre = nullptr;
                    
                    throw ;        // Queda preparado para hacer algo mejor
                } 
    InitBorder();    // ella se ocupa de : _LenSondPlus1= _LenSond + 1;    
}

void    ThDyAlign::Use (CSec  *sonde, CSec *target)    // antes hacia mas cosas.. porque?        //    ---------------------    Use            ---------
{    
    _sd = sonde  ;         _tg = target ;
    ResizeTable();
    _NNpar->UpdateGC( _sd->GCpercent(), _sd->Len() ,    _tg->GCpercent(), _tg->Len() );    // verificar 
}

void    ThDyAlign::InitBorder()                        //    --------------------------------------------------------------    InitBorder    ---------
{    float IniEnt = _NNpar->GetInitialEntropy();    //     -5.9f+_RlogC;    
    _LenSondPlus1= _LenSond + 1;        // NO QUITAR !! Las dim de la matriz son (ls+3)X(lt+3) : la pos 0  queda "reservada, y se considferan los 2 '$' (inic y final)                                        
    for (long i=0; i<=_LenSond; i++)                                // recorre la sonda "comparandola" con el '$' del target
    {    dH0(i,0)=dH1(i,0)=dH2(i,0)=        0.0f;                        // seguro ??????
        dS0(i,0)=dS1(i,0)=dS2(i,0)=        IniEnt;                     // seguro ??????
        pre(i,0)=                        st_1;                          // seguro ??????        NoStep; 
     }     
     for  (long j=1; j<=_LenTarg; j++)                                // recorre el target  "comparandolo" con el '$' de la sonda
     {    dH0(0,j)=dH1(0,j)=dH2(0,j)=        0.0f;                          // seguro ??????
        dS0(0,j)=dS1(0,j)=dS2(0,j)=        IniEnt;                      // seguro ??????
        pre(0,j)=                        st_7;                         /*        pre0(0,j) =         pre1(0,j) =         pre2(0,j) = */        // seguro ??????     NoStep;     
     }
}

void    ThDyAlign::Run(LonSecPos tg_j_StartPos /*= 1*/) // tg_j_StartPos    :  pos de com del target !? -------------------------    Run    ----------------
{   _maxglo = _InitMax;        // max del parametro rector :  como G o Tm
    _maxgloi  = 0;            
    _maxgloj  = 0;            //    _maxglot  = 0;            // ?????????????????// HACE FALTA ??????????????????

    _THits    = _HitsOK =0;        // Aqui ???

    
    register LonSecPos  i_1,  i,   j_1,  j ;        //  i y j en la DPMz usan las pos i-2 y j-2 en la sec !!!!!! 
    register Base        a_2,  a_1,   b_2,  b_1 ;    // las bases en esas posiciones -- exactamente como -- de la pos anterior
             Base        gap= bk2nu[ basek[0] ];        // basek[]="-ACGT$"            ,  // las 4 bases en el orden K. ?Coservar este orden  => b+cb=5 ??
    register float S0, S1, S2, H0, H1, H2 ;
//    long ls = _sd->_len + 2 , lt = _tg->_len + 2 ;    // lo mismo que _LenSond= LenSond +2 ;// Para que incluya los ´$´ ? y _LenTarg= LenTarg +2 ;
    Energy  forbidden_enthalpy = _NNpar->forbidden_enthalpy/10;        // MUY grande:  1e+18
    Entropy    forbidden_entropy  = _NNpar->forbidden_entropy/10 ;        // forbidden_entropy    = (_RlogC    ),    comparar con  IniEnt =     -5.9f+_RlogC;
    Entropy IniEnt = _NNpar->GetInitialEntropy();                    //     -5.9f+_RlogC;    
    float max1, max2, max0 ;
    Step st0,st1,st2; 
    b_1= gap ;                        // Las dim de la matriz son (lonS+3)X(lonT+3)=(ls+1)X(lt+1) : la pos 0  queda "reservada, y se consideran los 2 '$' (inic y final)
                                    // en la posicion i se concerva el resultado de "optimizar" sumandole a la i-1 el par i_2,i_1, donde 1 es la primera letra y ls pasa $
                                    // osea para el ultimo i=ls, el ultimo par es ls_2,ls_1 = lonS,lonS+1 = a$ donde a es la ultima letra de la sonda
    for (j = tg_j_StartPos ; j<= _LenTarg ; j++ )        // j=0 --> NoStep en la initTable
    {    b_2 = b_1 ;                  j_1 = j-1 ;    // b_2 corresponde a j-2.  Al inicio j=1: pos -1 !!, antes del '$' !! y se asume =0 osea '.'
        b_1   = bkn2c_nu    [_tg->_b [j_1]] ;    // b_1 corresponde a j-1.  Al inicio j=1: pos  0 !!, o sea '$' !! .   Comenzamos con '.$' .   // '.' y '$' y tambien '-' quedan igual
                                                // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
        a_1= gap ;                                    // y los ultimos: cuando i=ls: i_1=ls-1: $
        for (i = 1 ; i<= _LenSond ; i++ )
        {    //RestHS (i);
            a_2 = a_1 ;     i_1 = i-1 ;    // se justifica esto ????
            a_1 = _sd->_b [i_1] ;
            // -------------------------------------------------------- Diagonal cell: Alignment of sd[i] with tg[j]. Viniendo por la diagonal
            S0 = dS0 (i_1, j_1);        S1 = dS1 (i_1, j_1);        S2 = dS2 (i_1, j_1);    // al comienzo : 0,0 : i,0 : 0,j -> IniEnt
            H0 = dH0 (i_1, j_1);        H1 = dH1 (i_1, j_1);        H2 = dH2 (i_1, j_1);    // al comienzo : 0,0 : i,0 : 0,j -> 0
            if (H0<forbidden_enthalpy)        // OJO !!no tiene else ????? parar? // Diagonal  --- st_d,    st_4
            {    S0 += _NNpar->GetEntr (a_2, a_1, b_2, b_1);      // al comienzo : 1,1 -> ('.$','.$') -> 0: i,1 -> ('aa','.$') -> : 1,j -> ('.$','bb') -> 0 
                 H0 += _NNpar->GetEnth (a_2, a_1, b_2, b_1);     //   
            } 

            if (H1<forbidden_enthalpy)// OJO !caracter basek '-',' ' o gap en b, j no avanza y, a o i avanzaron -- st_a, st_3
            {    S1 += _NNpar->GetEntr (a_2, a_1, gap, b_1);        H1 += _NNpar->GetEnth (a_2, a_1, gap, b_1);       }     
            if (H2<forbidden_enthalpy) // hueco en a_1, i no avanza y, b o j avanzaron --- st_b, st_5
            {    S2 += _NNpar->GetEntr (gap, a_1, b_2, b_1);        H2 += _NNpar->GetEnth (gap, a_1, b_2, b_1);         }     

            float p0=(this->*CalcParam)(S0,H0);                float p1=(this->*CalcParam)(S1,H1);                float p2=(this->*CalcParam)(S2,H2);
            // por ejemplo : CalcTM    (float S,float H)const{return (S>=0 || H>=forbidden_enthalpy || (H/S)<0 ) ? 0         :  (H/S)    ; }
            // CalcParamG(S,H){return -_NNpar->CalcG (S, H);} CalcG (S,H,Ta){return (S>=0 || H>=forbidden_enthalpy/10000) ? -forbidden_freeEnerg : +(H - Ta*S);}
            if ((p0 >=p1)&&(p0 >=p2))            // Tm max  -> 0// Diagonal ? --- st_d pre(i_1,j_1) = st_3 ;    
            {                dS0(i,j) = S0;        dH0(i,j) = H0;        max0= p0;        st0 = st_4 ;    
            } else    if ((p1 >=p2))                // Tm max  -> 1
                    {        dS0(i,j) = S1;        dH0(i,j) = H1;        max0= p1;        st0  = st_3 ;    // Diagonal ? --- st_a st_3 
                    }    else                 // Tm max  -> 2
                        {    dS0(i,j) = S2;        dH0(i,j) = H2;        max0= p1;        st0  = st_5 ;    // Diagonal ? --- st_b= st_5
                        }
            
   //         if( dH0(i,j)>=forbidden_enthalpy/1000000 )    // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            //{    dS0(i,j) = IniEnt; dH0(i,j)=0 ; st0 =NoStep ;}

            if( dH0(i,j)==0 )    // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            {    dS0(i,j) = IniEnt;        /*st0 =NoStep ;*/}


            // -------------------------Cell superior o 4--------- Middle? cell: Alignment of b- seq1[i] with '-' en a // st_b
            // los gap son solo interiores
            // we do not allow $- in the beginning, therefore set to forbidden if  j=1    ?????????
            // also, disallow -$ in the end, therefore forbid if j=seq2len-1   . Seguro???????????
            if ( (j==1) ||  (j==_LenTarg-1) )            // pasar a ....?
            {    S0 =    S1 = _NNpar->forbidden_entropy;                H0 =     H1 = _NNpar->forbidden_enthalpy;
            } else 
            {    S0 = dS0 (i_1, j);                S1 = dS1 (i_1, j);
                H0 = dH0 (i_1, j);                H1 = dH1 (i_1, j);
                if (H0<forbidden_enthalpy)                    // OJO !! para que sirve si no tiene else ????? parar? pasar al otro ciclo ??
                {   S0 += _NNpar->GetEntr (a_2, a_1, b_1, gap);       H0 += _NNpar->GetEnth (a_2, a_1, b_1, gap);            }     
                if (H1<forbidden_enthalpy)
                {   S1 += _NNpar->GetEntr (a_2, a_1, gap, gap);        H1 += _NNpar->GetEnth (a_2, a_1, gap, gap);            }     
            }
            p0 = (this->*CalcParam)  ( S0, H0 );        p1 = (this->*CalcParam)  ( S1, H1 );
            if (p0 >=p1)            // Tm max  -> 0
            {    dS1(i,j) = S0;                dH1(i,j) = H0;                max1 = p0;                st1 = st_2 ;    // ? = st_2
            } else                    // Tm max  -> 1
            {    dS1(i,j) = S1;                dH1(i,j) = H1;                max1 = p1;                st1 = st_1 ;    //? = st_1
            }
   //         if( dH1(i,j)>=forbidden_enthalpy/1000000 )    // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            //{    dS1(i,j) = IniEnt; dH1(i,j)=0 ; st1 =NoStep ;}    

            if( dH1(i,j)==0 ) // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            {    dS1(i,j) = IniEnt;        /*st1 =NoStep ;*/}    
            

            // -------------------------------------- Lower cell: Alignment of '-' with seq2[j]
            // we do not allow $- in the beginning, therefore set to forbidden if  i=1    ?????????
            // also, disallow -$ in the end, therefore forbid if i=seq1len-1   . Seguro???????????
            if ( (i==1) ||  (i==_LenSond-1) )                    //  AQUI estaba el error MIO de las asimetrias ?!!!!!!!
            {    S0 =    S2 = _NNpar->forbidden_entropy;                H0 =     H2 = _NNpar->forbidden_enthalpy;
            } else 
            {    S0 = dS0 (i, j_1);                S2 = dS2 (i, j_1);
                H0 = dH0 (i, j_1);                H2 = dH2 (i, j_1);
                if (H0<forbidden_enthalpy)                    // OJO !! para que sirve si no tiene else ????? parar? pasar al otro ciclo ??
                {   S0 += _NNpar->GetEntr (a_1, gap, b_2, b_1);               H0 += _NNpar->GetEnth (a_1, gap, b_2, b_1);                }     
                if (H2<forbidden_enthalpy)
                {   S2 += _NNpar->GetEntr (gap, gap, b_2, b_1);                H2 += _NNpar->GetEnth (gap, gap, b_2, b_1);                }     
            }
            p0 = (this->*CalcParam)  ( S0, H0 );            p2 = (this->*CalcParam)  ( S2, H2 );            
            if (p0 >=p2)            // Tm max  -> 0
            {    dS2(i,j) = S0;                dH2(i,j) = H0;                max2 = p0;                st2 = st_6 ;    // ?
            } else                    // Tm max  -> 2
            {    dS2(i,j) = S2;                dH2(i,j) = H2;                max2 = p2;                st2 = st_7 ;    // ?
            }
   //         if( dH2(i,j)>=forbidden_enthalpy/1000000 )    // set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            //{    dS2(i,j) = IniEnt; dH2(i,j)=0 ; st2 =NoStep ;}    

            if( dH2(i,j)==0 )// set to zero if alignment so far is mismatches only (tm will be zero -> no need to reset maxglo)!
            {    dS2(i,j) = IniEnt;        /*st2 =NoStep ;*/}    

            {    float H,S;
            if(max0>=max1 && max0>=max2)
            { _max=max0; pre(i,j)=st0 ; H=dH0(i,j); S=dS0(i,j);}
            else    if (max1>max2)                // aqui se decide en que sec se insertan los gap falsos iniciales !!!!!!!! sd:'>' or tg:'>='
                    { _max=max1; pre(i,j)=st1;H=dH1(i,j); S=dS1(i,j);}
                    else
                        { _max=max2; pre(i,j)=st2;H=dH2(i,j); S=dS2(i,j);}

            if ( _NNpar->CalcTM (S, H) >= _Tm_sig && _NNpar->CalcG  (S, H) <= _G_sig)        //check if Hit found!      anadir G :  
            {     _THits++ ; if(AddIfHit(i,j)) _HitsOK++;    }    //_Hits.Add( new CHit(i,j,_max) );//check if local maximum found!
            }

            // check if global maximum found!
            // sin ((i==ls) || (j==lt))  seria un max mas global: sin el condicionamiento de llegar hasta el final de la sec.
            // A veces el final de la sec disminulle en lugar de aumentar la Tm, o el "Param", si tienen MM al final

            if (     (i==_LenSond)   &&  ( _max >=_maxglo)   )//  ?????????//     y no para ???  y else ???????????????????
                     if ( (_max >_maxglo) || (j    <_maxgloj))             {    _maxglo=_max;            _maxgloi=i;                _maxgloj=j;        }
            if (     (j==_LenTarg)   &&  ( _max >=_maxglo)   )//  ?????????//     y no para ???  y else ???????????????????
                     if ( (_max >_maxglo) || (i    <_maxgloi))             {    _maxglo=_max;            _maxgloi=i;                _maxgloj=j;        }

            //if (     ((i==ls) || (j==lt))   &&  ( _max >_maxglo)   )//  ?????????//     y no para ???  y else ???????????????????
            //{    _maxglo=_max;            _maxgloi=i;                _maxgloj=j;        }        
        }
    }
    _Ta_lastCalc= _NNpar->_Ta;
}

//bool    ThDyAlign::AddIfHit(long fi, long fj)   
bool    ThDyAlign_TmHits::AddIfHit(LonSecPos fi, LonSecPos fj)   // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
{    LonSecPos    i  = fi, j=fj, l=0, i0=-1, j0=-1;
    Energy        H= Get_H(i,j,step(i,j)) ,H0;
    Entropy        S= Get_S(i,j,step(i,j)) ,S0;
    Temperature    Th=0;                                    //Th= Get_Tm(i,j);
    Step st; 
    while (    Th<_Tm_sig)
    {    st= step(i,j);
        if (! st ) 
            break;
        i-= sti1[ st ];    // i,j ????,     i-= sti[ st ];  2 pasos de una vez   ---comparar !!!!!!!
        j-= stj1[ st ];

        if (i < 2 || j < 2 ) break ;                                    // i,j>1 ???? de todas formas no se pueden calc las Tm de las sondas
        S0= Get_S(i,j,step(i,j)) ;        H0=Get_H(i,j,step(i,j)) ;        // no se toma el maximo, sino segun el step
        Th= CalcParamTm( S-S0+_NNpar->GetInitialEntropy(),  H-H0);        // pos inic NOT INCLUSIV !!!!!!!!
        l+=1;                //l+=2;  2 pasos de una vez   ---comparar !!!!!!!
    }        

    Temperature Ti, Tj, Th_max=Th ;
    while (        i  >=2                        // esto puede no cumplirse nunca. De todas formas sale por el break
            &&    j  >=2
            /*&&  Th>=_Tm_sig*/)
    {    // ---fi, pf i y j en la DPMz usan las pos i-2, j-2, fi-2 y fj-2 en la sec !!!!!! ademas i - NOT INCLUSIV aqui
        Ti=_sd->Tm(i+1 -2, fi -2),                Tj=_tg->Tm(j+1 -2, fj -2) ;     
        if    ( Ti >= _Tm_max && Tj >= _Tm_max )        // ambas Tm de "sondas" demaciado altas (hibridacion "artificialmente" larga) 
            break ;
        if    (_Tm_min <= Ti && Ti <= _Tm_max)    // Tm de sonda en sd: en rango
            i0= i;                                // i0, j0 comienzo de la MAYOR de las sondas en rango terminada en fi,fj
        if    (_Tm_min <= Tj && Tj <= _Tm_max)    // Tm de sonda en tg: en rango
            j0= j;
        l+=1;//l+=2;  2 pasos de una vez   ---comparar !!!!!!!
        
        if( Th_max < Th ) Th_max=Th ;
        S0= Get_S(i,j,step(i,j)) ; H0=Get_H(i,j,step(i,j)) ;
        Th= CalcParamTm( S-S0+_NNpar->GetInitialEntropy(),  H-H0);
    // i,j ????
        st= step(i,j);
        if (! st ) 
            break;
        i-= sti1[ st ];
        j-= stj1[ st ];

    }

    if (i0==-1 && j0==-1) return false;
    _Hits.emplace_back( /*new CHit(*/ fi, fj, i0, j0,l, H-H0,S-S0,Th_max,step(fi,fj) /*)*/ );        //    i,j ???? Th_max??
    return true;
}
ThDyAlign::~ThDyAlign()
{    delete[] _dH0;
    delete[] _dH1;
    delete[] _dH2;
    delete[] _dS0;
    delete[] _dS1;
    delete[] _dS2;    
    delete[] _pre; // seguro ????????
    ClearHits() ;  // seguro ????????

}

Energy        ThDyAlign::Get_H         (LonSecPos i, LonSecPos j, Step st) const    // ------alternativa a  Get_X    --    da param en (i,j) segun step()i,j) !!
{    switch (st)
    {    case st_0:     return 0; 
        case st_1:
        case st_2:    return dH1(i,j);
        case st_3:    
        case st_4:    
        case st_5:    return dH0(i,j); 
        case st_6:    
        case st_7:    return dH2(i,j); 
        default: 
         return 0; 
    }
}
Entropy        ThDyAlign::Get_S         (LonSecPos i, LonSecPos j, Step st) const
{    switch (st)
    {    case st_0:     return 0; 
        case st_1:
        case st_2:    return dS1(i,j);
        case st_3:    
        case st_4:    
        case st_5:    return dS0(i,j); 
        case st_6:    
        case st_7:    return dS2(i,j); 
                default: return 0; 
    }
}
ThDyAlign::Step    ThDyAlign::Get_pre     (LonSecPos i, LonSecPos j, Step st) const 
{    
    return step(i,j);
    /*switch (st)
    {    case st_0:     return NoStep; 
        case st_1:
        case st_2:    return pre1(i,j);
        case st_3:    
        case st_4:    
        case st_5:    return pre0(i,j); 
        case st_6:    
        case st_7:    return pre2(i,j); 
                default: return NoStep; 
    }*/
}
void        ThDyAlign::SelectOptParam(LonSecPos i, LonSecPos j, Temperature Ta )
{    switch (step(i,j))
    {    case st_0:     _optS =0;    _optH = 0;    _optTm= 0 ; return ;        // poner valor "absurdo", nulo, o error o accert o expection
        case st_1:
        case st_2:    _optS = dS1(i,j);    _optH = dH1(i,j);    break ;
        case st_3:    
        case st_4:    
        case st_5:    _optS = dS0(i,j);    _optH = dH0(i,j);    break ;
        case st_6:    
        case st_7:    _optS = dS2(i,j);    _optH = dH2(i,j);    break ;
                default: return ;                                         // poner valor "absurdo", nulo, o error o accert o expection
    }
    _optTm = _NNpar->CalcTM( _optS  , _optH); 
    _optG  = _NNpar->CalcG ( _optS,   _optH ,Ta) ; // +(_optH-Ta*_optS);  
}

CHit        *ThDyAlign::GetOptHit()
{return new CHit(*this);}

Energy        ThDyAlign::Get_H_max_para    (LonSecPos i, LonSecPos j) const    // --------- Get_H    -----    da param en (i,j) para Tm max !!!!    ?????????
{    Temperature tm0 = _NNpar->CalcTM( dS0(i,j)  ,dH0(i,j));
    Temperature tm1 = _NNpar->CalcTM( dS1(i,j)  ,dH1(i,j));
    Temperature tm2 = _NNpar->CalcTM( dS2(i,j)  ,dH2(i,j));

    if ((tm0 >=tm1)&&(tm0 >=tm2))        return dH0(i,j);    // Tm max  -> 0
    if ( tm1 >=tm2)                        return dH1(i,j);    // Tm max  -> 1
                                        return dH2(i,j);    // Tm max  -> 2
}
Entropy        ThDyAlign::Get_S_max_para    (LonSecPos i, LonSecPos j) const    // --------- Get_S    -----    da param en (i,j) para Tm max !!!!    ?????????
{    Temperature tm0 = _NNpar->CalcTM( dS0(i,j)  ,dH0(i,j));
    Temperature tm1 = _NNpar->CalcTM( dS1(i,j)  ,dH1(i,j));
    Temperature tm2 = _NNpar->CalcTM( dS2(i,j)  ,dH2(i,j));

    if ((tm0 >=tm1)&&(tm0 >=tm2))          return dS0(i,j);    // Tm max  -> 0
    if ( tm1 >=tm2)                        return dS1(i,j);    // Tm max  -> 1
                                        return dS2(i,j);    // Tm max  -> 2
}
Temperature    ThDyAlign::Get_Tm_max_para    (LonSecPos i, LonSecPos j) const// --------- Get_Tm    -----    da param en (i,j) para Tm max !!!!    ?????????
{    Temperature tm0 = _NNpar->CalcTM( dS0(i,j)  ,dH0(i,j));
    Temperature tm1 = _NNpar->CalcTM( dS1(i,j)  ,dH1(i,j));
    Temperature tm2 = _NNpar->CalcTM( dS2(i,j)  ,dH2(i,j));

    if ((tm0 >=tm1)&&(tm0 >=tm2))          return tm0<0 ? 0: tm0;    // Tm max  -> 0
    if (tm1 >=tm2)                        return tm1<0 ? 0: tm1;    // Tm max  -> 1
                                        return tm2<0 ? 0: tm2;    // Tm max  -> 2
}
Energy        ThDyAlign::Get_G_max_para    (LonSecPos i, LonSecPos j, Temperature Ta )const // ------- + Get_G    ----    da param en (i,j) para Tm max !!!!    ?????????
{    Temperature tm0 = _NNpar->CalcTM( dS0(i,j)  ,dH0(i,j));
    Temperature tm1 = _NNpar->CalcTM( dS1(i,j)  ,dH1(i,j));
    Temperature tm2 = _NNpar->CalcTM( dS2(i,j)  ,dH2(i,j));

    if ((tm0 >=tm1)&&(tm0 >=tm2))         return +(dH0(i,j)-Ta*dS0(i,j));        // Tm max  -> 0
    if (tm1 >=tm2)                        return +(dH1(i,j)-Ta*dS1(i,j));        // Tm max  -> 1
                                        return +(dH2(i,j)-Ta*dS2(i,j));        // Tm max  -> 2
}
void        ThDyAlign::SelectOptParam_max_para(LonSecPos i, LonSecPos j, Temperature Ta )// ---- SelectOptParam ---    set  _optParam en (i,j) para Tm max !!!!    ?????????
{    Temperature tm0 = _NNpar->CalcTM(dS0(i,j),dH0(i,j));                        // determine optimum values for dG and dH
    Temperature tm1 = _NNpar->CalcTM(dS1(i,j),dH1(i,j));
    Temperature tm2 = _NNpar->CalcTM(dS2(i,j),dH2(i,j));

            if ((tm0 >=tm1)&&(tm0 >=tm2))    { _optS = dS0(i,j);    _optH = dH0(i,j);    _optTm= tm0 ;    }         // Tm max  -> 0
    else    if (tm1 >=tm2)                    { _optS = dS1(i,j);    _optH = dH1(i,j);    _optTm= tm1 ;    }        // Tm max  -> 1 
    else                                     { _optS = dS2(i,j);    _optH = dH2(i,j);    _optTm= tm2 ;    }        // Tm max  -> 2
    _optG  = +(_optH-Ta*_optS);    //     _optTm = _optTm<0 ? 0: _optTm;
}



void    FracTDAlign::BeginAlign (CSec  *sonde, CSec *target)    
{    CalcParam = &ThDyAlign::CalcParamTm;// "Tm" alignment !!! Verificar entonces si se necec lo del funt pointer  !!!!!!!!!!!!!!!!!
    _InitMax=_NNpar->kein_Tm ;   
    if (_Ta==_NNpar->kein_Tm) _Ta=_NNpar->_Ta;// temperatura del ultimo "exp" (o calculo!!!)     def:_Ta(CtoK(60)) 
    else        _NNpar->SetTa(_Ta);            //  !!!!!!
                                        
    ThDyAlign::Align (sonde, target);                    // "Tm" alignment !!! Verificar entonces si se necec lo del funt pointer  !!!!!!!!!!!!!!!!!

    CalcParam = &ThDyAlign::CalcParamG;// "Tm" alignment !!! Verificar entonces si se necec lo del funt pointer  !!!!!!!!!!!!!!!!!
    _InitMax=_NNpar->forbidden_freeEnerg ;    
    _finisch = false ;
    _iterations = 1 ;
}

void    FracTDAlign::iterate        (Temperature ta)
{    assert ( _iterations ); 
    assert ( _fixedNumIter>=0 ); 

    _NNpar->SetTa (ta);                //    InitBorder    ();
    Run();
    SelectOptParam();        // Aqui se selecciona una nueva Tm
    
    _iterations++;
}

bool    FracTDAlign::NotFinisch    () 
{    
    if (    abs( G() ) < _maxG_der /*|| abs( G(Tm()) ) < _maxG_der */    )    _finisch=true ; // ya estamos demaciado cerca de Tm
    if (    IsEq(last_Ta(),Tm(),0.00001f)    )    _finisch=true ; // ya estamos demaciado cerca de Tm - definir mejor !!!!!!!
    if (    _iterations>=_maxNumIt    )    _finisch=true ; // forced _maxNumIt 
    if (  _fixedNumIter ) 
    if (  _iterations>=_fixedNumIter)    _finisch=true ;
        else                            _finisch=false; // forced _fixedNumIter     
    
    return !_finisch; 
}

CHit::CHit (ThDyAlign &Al) : _i(Al._maxgloi), _j(Al._maxgloj), _max(Al._maxglo) // Hit optimo
{    Al.SelectOptParam();
    _H = Al.H();
    _S = Al.S();
    _G = Al.G(); // a la Ta en que se calculo
    _Tm= Al.Tm(); ;
}

void CHitAligned::ExtractAligment(ThDyAlign &Al)
            // Las dim de la matriz son (lonS+3)X(lonT+3)=(ls+1)X(lt+1) : la pos 0  queda "reservada, y se consideran los 2 '$' (inic y final)
            // en la posicion i se concerva el resultado de "optimizar" sumandole a la i-1 el par i_2,i_1, donde 1 es la primera letra y ls pasa $
            // osea para el ultimo i=ls, el ultimo par es ls_2,ls_1 = lonS,lonS+1 = a$i donde a es la ultima letra de la sonda
            // j=0 --> NoStep en la initTable
            // b_2 corresponde a j-2.  Al inicio j=1: pos -1 !!, antes del '$' !! y se asume =0 osea '.'
            // b_1 corresponde a j-1.  Al inicio j=1: pos  0 !!, o sea '$' !! .   Comenzamos con '.$' .   // '.' y '$' y tambien '-' quedan igual
            // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
            // y los ultimos: cuando i=ls: i_1=ls-1: $

{    
     //_sd;  _tg ;  _st ;
    _mt= _mm= _sgap= _tgap= 0 ;    // count sonde and target - matchs , mistmatch, and gaps
    ThDyAlign::Step st; 
    Base gap=     nu2dgba    [0] ; //    nu2dgba    []="-GCSTKYBARMVWDHN"    ,  // "convierte" el "codigo numerico" en letra del deg cod; lo contrario de db2nu[]. de "numero a 

    for (_l=0, _i0=_i, _j0=_j; _i0 >=1 && _j0>=1 ; _l++)        // >=0   or   >0   ????  da el $  ??      Calcula longitud del alignment
        if (! ( st= Al.step(_i0, _j0) ) ) { break;}        //_l++; no funciono
        else     {    _i0 -= ThDyAlign::sti1[ st ];  _j0 -= ThDyAlign::stj1[ st ];}  // i-= sti[ st ];  2 pasos de una vez ---comparar !!            

    _sd.resize( _l+1) ;     
    _tg.resize( _l+1) ; 
    _st.resize( _l  ) ; 

    _sd[_l] =                 _tg[_l]  = 0 ;

    for (long a=_l-1, i=_i, j=_j ; a >=0  ; a--)
    {    _st[a]= st =Al.step(i, j);

        if ( ThDyAlign::sti1[ st ] )        {_sd[a]= (*Al._sd)[--i]    ;            }    // {_sd[a]=Al._sd->_c[i]    ; i-- ;        }
        else                                {_sd[a]=gap                ; _sgap++;    }

        if ( ThDyAlign::stj1[ st ] )        {_tg[a]=(*Al._tg)[--j]    ;            }    // {_tg[a]=Al._tg->_c[j]    ; j-- ;        }
        else                                {_tg[a]=gap                ; _tgap++;    }
    
        if (_sd[a] == _tg[a])    {if (_sd[a]!=basek[5])            _mt++ ;    }    //    basek[]=".ACGT$"    ,  // las 4 bases en el orden de Kadelari. 
        else                                                    _mm++ ;
    }

    ;
}


void CHitAligned::ReCalcule( std::shared_ptr<CSaltCorrNN>  NNpar )
{     
    Base a_1 , 
         a    =bk2nu[_sd[0]], 
         b_1 , 
         b    =bk2c_nu[_tg[0]] ;

    _Hr= 0; 
    _Sr= NNpar->GetInitialEntropy();

    for (_i=_j= _i0= _j0=0 ;  _i+LonSecPos( 1 ) <_sd.length() && _j+1<_tg.length()    ;  ++_i, ++_j )
    {    
        a_1 = a; 
        a   = bk2nu  [ _sd[_i+1] ];    
        b_1 = b; 
        b   = bk2c_nu[ _tg[_j+1] ];

        _Sr += NNpar->GetEntr (a_1, a, b_1, b); 
        _Hr += NNpar->GetEnth (a_1, a, b_1, b); 
    }
    _l    = _i ;
    _Tmr  = NNpar->CalcTM( _Sr, _Hr) ;
    _Gr      = NNpar->CalcG    ( _Sr, _Hr) ;
}


void CMSecCand::Use(CMultSec* MSec)
{    _MSec=MSec; 
    if (! _TDATmC ) 
        _TDATmC.reset(new ThDyAlign_TmCand(_MSec->_Global._Len.Max(), _MSec->_NNPar));  // Por  que????
    else    _TDATmC->ResizeTable(  _MSec->_Global._Len.Max()  ,   _MSec->_Global._Len.Max()   ) ;
    
    _TNumPosCand= _TNumCand = 0;        // AQUI ???????????
    _NumPosCand = _NumCand  = 0;        // AQUI ???????????

    _TDATmC->SetSig(_Tm_sig, _G_sig);    //_TDATmC->SetTmLimits(_Tm_sig, _Tm_min, _Tm_max);
}


std::unique_ptr<CSecCand> CMSecCand::AddBeging(CSec &sec)
{    
    auto newtg = std::make_unique<CSecCand>(sec, _sL);
    _TNumCand    +=newtg->_NumCand;    // sobreestima la cantidad total de candidatos (porque se pueden repetir en varias sec.)
    _TNumPosCand+=newtg->_NumPosCand;
    _NumCand    +=newtg->_NumCand;    // sobreestima la cantidad total de candidatos (porque se pueden repetir en varias sec.)
    _NumPosCand +=newtg->_NumPosCand;
    ++_NSecCand ;
    return newtg;
}

void   CMSecCand::FindCommon    (CSecCand  &newtg, CSecCand &curtg, bool design)    
{    long     dec_NumCand        = newtg._NumCand    + curtg._NumCand ; 
    long     dec_NumPosCand    = newtg._NumPosCand + curtg._NumPosCand ; 

    _TDATmC->FindCommon(&newtg, &curtg, design );

    if (design)                                    //  REVISAR   !!!!!!!!!!!!!!!!!!!!!!!!!
    {
    dec_NumCand        -= newtg._NumCand     + curtg._NumCand ; 
     dec_NumPosCand    -= newtg._NumPosCand + curtg._NumPosCand ; 

        _TNumCand        -=  dec_NumCand ;   
        _NumCand        -=  dec_NumCand ;   
        _TNumPosCand    -=    dec_NumPosCand ;   
        _NumPosCand        -=    dec_NumPosCand ;   
    }
        // sobreestima la cantidad total de candidatos (porque se pueden repetir en varias sec.)
        //_TNumPosCand+=curtg._NumPosCand;
}

/// use sec to create a new CSecCand, compare it with all the previous and the append to the list
CSecCand *CMSecCand::Add(CSec &sec)
{    
    CSecCand *newtg = new CSecCand( sec ,    _sL        );/*            _G_min,        _G_max,                                 
                                                                    _Tm_min,    _Tm_max,
                                                                    _L_min,        _L_max    */    

    _TNumCand    +=newtg->_NumCand;    // sobreestima la cantidad total de candidatos (porque se pueden repetir en varias sec.)
    _TNumPosCand+=newtg->_NumPosCand;
    _NumCand    +=newtg->_NumCand;    // sobreestima la cantidad total de candidatos (porque se pueden repetir en varias sec.)
    _NumPosCand +=newtg->_NumPosCand;

    for ( auto &Cur:  _LSecCand )     
        _TDATmC->FindCommon(newtg, Cur.get() );

    _LSecCand.emplace_back(newtg);

    return newtg;
}
        //_TNumCand    +=tg->_NumCand;
        //_TNumPosCand+=tg->_NumPosCand;
    //_TNumCand    +=newtg->_NumCand;
    //_TNumPosCand+=newtg->_NumPosCand;


/// ---fi (sonda) y fj (target) en la DPMz usan las pos finales-inclus fi-2 y fj-2 en la sec !!!!!!
bool    ThDyAlign_TmCand::AddIfHit(long fi, long fj)
{    
    if (fi<2 || fj<2) return false;   assert (fi>=2);assert (fj>=2);// como pueden ser < 2 y tener Tm > sig ???? 

    auto    &ri = _cs->_rg[fi-2] ,   // Candidato - Sonda  -termina en fi  (o sea -- i )
            &rj = _ct->_rg[fj-2] ;   // Candidato - Target -termina en fj  (o sea -- j )
    if(! ri && ! rj) return false;                        // En al menos una de las sec todavia esta disp
//  puede haber una sec intermedia bena para todas---    // algun cand en estas pos. Eso incluye que no estan demaciado cerca del comienzo de las sec.
        // solucion TEMPORAL, no muy eficiente : mejor duplicar codigo para cada caso de que una de los rang=0
        
    LonSecPos i  = fi, j=fj    /*,     l=0, i0=-1, j0=-1*/;    
    if (ri && rj)    
    {    CRangBase &Ri(*ri),&Rj(*rj); Ri.schift(2-1); Rj.schift(2-1);
        //CRangBaseSchift sRi(*ri, 2-1),sRj(*rj, 2-1); //use ri, rj directly. Dont worry more about schifting

        Step st;                             // i-sonda, j-target. pos de comienzo de la zona de hibrid analizandose
        while (i > Ri.Max() && j > Rj.Max() )    // saltar rapido el tramo de Al con corresp sondas "cortas" en ambas sec  (i > pfi && j > pfj )    
        {    st= step(i,j);                    // hasta que alguno de los dos entre en rango
            if (! st ) 
                {Ri.schift(-2+1);Rj.schift(-2+1);return false;}
            i-= sti1[ st ];    //     i-= sti[ st ];  2 pasos de una vez   ---comparar !!!!!!!
            j-= stj1[ st ];        //            l+=1;    
        }
        Energy        H= Get_H(fi,fj,step(fi,fj)) ,H0,Gh;   // Pasar S y H como argumentos ???? - acaban de ser calculados
        Entropy        S= Get_S(fi,fj,step(fi,fj)), S0  , S00=_NNpar->GetInitialEntropy();
        Temperature    last_Th_OK=0, Th=0;        
        while (i >= Ri.Min() || j >= Rj.Min() )            // mientras que alguno de los dos este en rango//assert(i < 0 || j < 0 );
        {    
            if (i < 0 || j < 0 ) break ;            // y el otro no se haga "neg." no debe pasar nunca. Dejarlo en 1 ?, ...or 2  ???
            S0= Get_S(i,j,step(i,j)) ;
            H0= Get_H(i,j,step(i,j)) ;                // no se toma el maximo, sino segun el step ??????????????????
            Th= CalcParamTm( S-S0 +S00,  H-H0); Gh= CalcParamG( S-S0 +S00,  H-H0);  // pos inic NOT INCLUSIV !!!!!!!!
            if (Th>_Tm_sig   && Gh<_G_sig)                            // Th -OK
            {    last_Th_OK=Th;
                Ri.addMatch(i);
                Rj.addMatch(j);
            }                        //l+=1;        //l+=2;  2 pasos de una vez   ---comparar !!!!!!!        
            st= step(i,j);
            if (! st ) 
                break;
            i-= sti1[ st ];
            j-= stj1[ st ];
        }
        if (last_Th_OK==0) {Ri.schift(-2+1);Rj.schift(-2+1);return false;}
        Ri.schift(-2+1);
        Rj.schift(-2+1);
        return true;
    }
    if (!rj)    
    {    assert (ri);            //   ?????
        CRangBase &Ri(*ri); Ri.schift(2-1); 
        Step st;     
        while (i > Ri.Max() && j > 0 ) // i && j > pfj )            // saltar rapido el tramo de Al con corresp sonda i "corta" 
        {    st= step(i,j);                                    // hasta que i entre en rango
            if (! st ) 
                {Ri.schift(-2+1);return false;}
            i-= sti1[ st ];    // i,j ????,     i-= sti[ st ];  2 pasos de una vez   ---comparar !!!!!!!
            j-= stj1[ st ];    //            l+=1;    
        }
        float H= Get_H(fi,fj,step(fi,fj)) ,H0, last_Th_OK=0,  // Pasar S y H como argumentos ???? - acaban de ser calculados
                S= Get_S(fi,fj,step(fi,fj)), S0, Th=0,Gh          , S00=_NNpar->GetInitialEntropy();

        while (i >=  Ri.Min() )            // mientras que alguno de los dos este en rango
        {    if (i < 0 || j < 0 ) break ;    // y el otro no se haga "neg."
            S0= Get_S(i,j,step(i,j)) ;
            H0= Get_H(i,j,step(i,j)) ;        // no se toma el maximo, sino segun el step ??????????????????

            Th= CalcParamTm( S-S0 +S00,  H-H0); Gh= CalcParamG( S-S0 +S00,  H-H0);  // pos inic NOT INCLUSIV !!!!!!!!
            if (Th>_Tm_sig   && Gh<_G_sig)                            // Th -OK
            {    last_Th_OK=Th;
                Ri.addMatch(i);
            }            //            l+=1;                //l+=2;  2 pasos de una vez   ---comparar !!!!!!!        
            st= step(i,j);
            if (! st ) 
                break;
            i-= sti1[ st ];
            j-= stj1[ st ];
        }
        if (last_Th_OK==0) {Ri.schift(-2+1);return false;}
        Ri.schift(-2+1);
        return true;

    } else 
    {    // ri=rj;

        CRangBase &Rj(*rj);  Rj.schift(2-1);

        Step st;     
        while (j > Rj.Max() && i > 0) // i && j > pfj )            // saltar rapido el tramo de Al con corresp sonda i "corta" 
        {    st= step(i,j);                                    // hasta que i entre en rango
            if (! st ) 
                {Rj.schift(-2+1);return false;}
            i-= sti1[ st ];    // i,j ????,     i-= sti[ st ];  2 pasos de una vez   ---comparar !!!!!!!
            j-= stj1[ st ];            //            l+=1;    
        }
        float H= Get_H(fi,fj,step(fi,fj)) ,H0, last_Th_OK=0,  // Pasar S y H como argumentos ???? - acaban de ser calculados
                S= Get_S(fi,fj,step(fi,fj)), S0, Th=0,Gh          , S00=_NNpar->GetInitialEntropy();

        while (j >= Rj.Min() )            // mientras que alguno de los dos este en rango
        {    if (i < 0 || j < 0 ) break ;    // y el otro no se haga "neg."
            S0= Get_S(i,j,step(i,j)) ;
            H0= Get_H(i,j,step(i,j)) ;        // no se toma el maximo, sino segun el step ??????????????????

            Th= CalcParamTm( S-S0 +S00,  H-H0); Gh= CalcParamG( S-S0 +S00,  H-H0);  // pos inic NOT INCLUSIV !!!!!!!!
            if (Th>_Tm_sig   && Gh<_G_sig)                            // Th -OK
            {    last_Th_OK=Th;
                Rj.addMatch(j);
            }            //            l+=1;                //l+=2;  2 pasos de una vez   ---comparar !!!!!!!        
            st= step(i,j);
            if (! st ) 
                break;
            i-= sti1[ st ];
            j-= stj1[ st ];
        }
        if (last_Th_OK==0)     {Rj.schift(-2+1);return false;}
        Rj.schift(-2+1);

        return true;
    }
}

        //LonSecPos   pii, pij, pfi, pfj;        // i-sonda, j-target. En este rango la Tm y "calidad" de las sondas esta garantizada desde el principio
        //LonSecPos  cpii,cpij,cpfi,cpfj;        // c-current - solo resta comprobar la Th. Todas estas son pos inic. Por eso el -1. Mejor cambiar ??? que??
        // pii=ri->_pi   +2-1;        // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
        //cpii=ri->_picur+2-1;        // pos inic NOT INCLUSIV !!!!!!!! y en la sec si !!!
        // pij=rj->_pi   +2-1;        //    pii      pfi               fi
        //cpij=rj->_picur+2-1;        //----|++++++++|-----------------|--------
        //                                    //    pii     pfi               fi
        //                                    //----|++++++++|-----------------|--------            El rango inicial, y como va quedando
        //                                    //   cpfi                       fi                    El rango para calculo ("cur"), antes del comienzo    
        //                                    //----|---------|----------------|--------            Asi se queda si no hibridan entre si las sec en esta zona,
        //                                    //             cpii             fi                    y entonces "colapsa" el rango
        //                                    //            cpfi              fi                    En este caso encontro 5 "cand" comunes    
        //                                    //--------|+++|------------------|--------            
        //                                    //       cpii                   fi                    
        // pfi=ri->_pf   +2-1;
        //cpfi=ri->_pfcur+2-1;
        // pfj=rj->_pf   +2-1;
        //cpfj=rj->_pfcur+2-1;
        //assert(pii>=2);assert(cpii>=2);assert(pij>=2);assert(cpij>=2);
        //assert(pfi>=2);assert(cpfi>=1);assert(pfj>=2);assert(cpfj>=1);
                //if    (pii<=i && i<=pfi)                // y sonda 1 tambien  -- modif los curr rg
                //{            if (cpfi<i)
                //                cpfi=i;
                //    if (i<cpii)
                //        cpii=i;
                //    Ri.adjustCur(i);// ++ri->match[i-pii];            // OTRO Match !!!!!!
                //}
                //if  (pij<=j && j<=pfj)                // y sonda 2 tambien  -- modif los curr rg
                //{            if (cpfj<j)
                //                cpfj=j;
                //    if (j<cpij)
                //        cpij=j;
                //    Rj.adjustCur(j);//++rj->match[j-pij];            // OTRO Match !!!!!!
                //}
        //assert(pii>=2);assert(cpii>=2);assert(pij>=2);assert(cpij>=2);
        //assert(pfi>=2);assert(cpfi>=2);assert(pfj>=2);assert(cpfj>=2);
        //ri->_picur    = cpii  -2+1 ;                // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
        //rj->_picur    = cpij    -2+1 ;        // pos inic NOT INCLUSIV !!!!!!!! y en la sec si !!!
        //ri->_pfcur    = cpfi    -2+1 ;        //    pii      pfi               fi
        //rj->_pfcur    = cpfj    -2+1 ;        //----|++++++++|-----------------|--------
        ////_Hits.Add( new CHit( fi, fj, i0, j0,l, H-H0,S-S0,Th_max,step(fi,fj) ) );        //    i,j ???? Th_max??
        //pi=ri->_pi   +2-1;        // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
     //  cpi=ri->_picur+2-1;        // pos inic NOT INCLUSIV !!!!!!!! y en la sec si !!!
        //pf=ri->_pf   +2-1;
     //  cpf=ri->_pfcur+2-1;
        //LonSecPos   pi, pf, cpi,cpf;                    
                //if    (pi<=i && i<=pf)                // y sonda 1 tambien  -- modif los curr rg
                //{            if (cpf<i)
                //                cpf=i;
                //    if (i<cpi)
                //        cpi=i;
                //    //++ri->match[i-pi];            // OTRO Match !!!!!!
                //}
        //ri->_picur    = cpi  -2+1 ;                // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
        //ri->_pfcur    = cpf  -2+1 ;        //    pii      pfi               fi
        //long   pi, pf, cpi,cpf;                    
        //pi=rj->_pi   +2-1;        // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
     //  cpi=rj->_picur+2-1;        // pos inic NOT INCLUSIV !!!!!!!! y en la sec si !!!
        //pf=rj->_pf   +2-1;
     //  cpf=rj->_pfcur+2-1;
                //if    (pi<=j && j<=pf)                // y sonda 1 tambien  -- modif los curr rg
                //{            if (cpf<j)
                //                cpf=j;
                //    if (j<cpi)
                //        cpi=j;
                //    //++rj->match[j-pi];            // OTRO Match !!!!!!
                //}    
        //rj->_picur    = cpi   -2+1 ;                // ---fi y fj en la DPMz usan las pos fi-2 y fj-2 en la sec !!!!!!
        //rj->_pfcur    = cpf    -2+1 ;        //    pii      pfi               fi

void    ThDyAlign::Export_Hits(ofstream &osHits, char *sep)        // mientras estan conectados al Al
{    
    osHits    << endl        << "Al Len"        << sep  <<"Tm hib"    
            << sep        << "de (i0)"    << sep    << "a (i) "        << sep  <<"Tm snd" 
            << sep        << "j0"            << sep    << "j "            << sep  <<"Tm tgt"     ;


    for ( const CHit &h : _Hits ) //(CHit *)_Hits.goBeging() ; _Hits.NotEnd() ; h=(CHit *)_Hits.goNext())    
    {    
        osHits    << endl    << h._l        << sep << KtoC(h._max)    
                << sep     << h._i0 -1    << sep << h._i -2        << sep  ;
                
                if (h._i0<1)                                    osHits    << "No Tm!!" ;
                else                                            osHits    <<   KtoC(_sd->Tm(h._i0+1 -2, h._i -2 )) ;

        osHits    << sep     << h._j0 -1    << sep << h._j -2        << sep ;
                if (h._j0<1)                                    osHits    << "No Tm!!" ;
                else                                            osHits    <<   KtoC(_tg->Tm(h._j0+1 -2, h._j -2 )) ;                 
    }
}


void CMSecCand::write_probes(    CMultSec         *res       ,
                    const std::string &fileName , 
                    fileFormat         format   /* = fileFormat::fasta*/)
{
    bool    f_fas = format & fileFormat::fasta ,
            f_csv = format & fileFormat::csv   ;

    ofstream    osFasta, osCSV;
    char sep[]=";";                            // cambiar a global sep

    if (f_fas)
        osFasta.open(fileName + ".sonden.fasta");

    if (f_csv)
    {     
        osCSV.open(fileName + ".sonden.csv");
        osCSV << "Num" <<sep<< "SecName" <<sep<<    "Inic" <<sep<< "Fin" <<sep<<    "Len"  <<sep<< "Tm" <<sep<< "Sec"
                << sep <<"H"<< sep <<"S"<< sep <<"G(Ta=" << KtoC(_TDATmC->Ta()) << " gr)" << sep << "No.matchs";
    }

    //if (f_csv)
    //    osCSV << endl << Num << sep << s._Sec.Name() << '.' << pi << '.' << fi << sep << pi << sep << fi << sep << fi - pi + 1 << sep
    //    << KtoC(s._Sec.Tm(pi, fi)) << sep << cur_s
    //    << sep << cand._SdH[cand.Len() - 1] << sep << cand._SdS[cand.Len() - 1] << sep << cand.G() << sep << matchs + 1;

    //if (f_fas)
    //    osFasta << endl << '>' << s._Sec.Name() << '.' << pi << '.' << fi << "  ; Tm=" << KtoC(s._Sec.Tm(pi, fi))
    //    << endl << cur_s;


}

void    CMSecCand::ExportCommonSonden(  bool             colpased,
                                        NumRang<float>   ExtrCovPerc, 
                                        CMultSec         *res       /*= {}*/,
                                        const std::string &fileName /*= ""*/, 
                                        fileFormat         format   /* = fileFormat::fasta*/)
{    

    if (!res) return; ///\todo temp !

    auto n = static_cast<double>(_NSecCand - 1);
    NumRang<int> ExtrCov (static_cast<int>((n * ExtrCovPerc.Min()) /100.0)  , static_cast<int>((n * ExtrCovPerc.Max()) /100.0 )) ;
    

    set <string> CandSet;   //  grant no sequence repetition

    //  perform self alignment to filter possible dimers 
    FracTDAlign fAl( _sL._L.Max() + 1 ,  _sL._L.Max() + 1, _TDATmC->_NNpar);
    fAl.SetTa ( _TDATmC->Ta () );
    long Num{0};
    std::unique_ptr<CSec> cand, c_cand;

/// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::cout << "\n EXPORTING...";
/// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (auto &Cur : _LSecCand)    // all targets
    {    
        CSecCand &s=*Cur ;
        long l= s._Sec.Len() ;
        for (long fi=1 ; fi<=l;fi++)           // for each position in the target (were the probe candidate end in each rang)
        {    
            auto &r=s._rg[fi] ;                /// take the rang of survived probe candidates in this position   \todo review !!!
            if (! r) continue;   
            for (long pi=r->Min(); pi <=r->Max();pi++)    // for each candidate initial position in the rang
            {    
                assert(pi<l);
                int matchs=r->matchs[pi- r->Min()];

                if ( (colpased || !ExtrCov.isIntern (matchs) ) && CandSet.insert(s._Sec.Copy_Seq(pi, fi)).second )
                {    
                    cand  .reset(s._Sec.Clone(pi, fi, DNAstrand::direct)); // (cur_s  , 1,"s", _TDATmC->_NNpar);     
                    c_cand.reset(s._Sec.Clone(pi, fi, DNAstrand::compl ));   
                    
                    fAl.Align(cand.get(), c_cand.get() );
                    ++Num;
                    if (fAl.Tm() < _MaxSelfTm && fAl.G(_TDATmC->Ta()) > _MinSelfG)
                    {
                    // anadir otras comprobaciones: estruct secund (self-align-compl)

                    /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //  std::cout<< "\nSeq: "<<  
                    /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
                        res->AddSec(cand.release())->Name();



                    //if (f_csv)    
                    //    osCSV <<endl<<     Num<<sep<<    s._Sec.Name()<<'.'<<pi<<'.'<<fi<<sep<<    pi<<sep<<     fi<<sep<<    fi-pi+1<<sep
                    //          << KtoC(s._Sec.Tm(pi,fi)) <<sep<< cur_s 
                    //          << sep <<cand._SdH[cand.Len ()-1] << sep <<cand._SdS[cand.Len ()-1]<< sep <<cand.G() << sep << matchs+1; 
                    //if (f_fas)
                    //    osFasta    <<endl << '>' << s._Sec.Name()<<'.'<<pi<<'.'<<fi    <<"  ; Tm="<< KtoC(s._Sec.Tm(pi,fi))  
                    //            <<endl<< cur_s  ; 

                    }
                }
            }
        }

    }
}//            

void print_ThDyAlign (ofstream &osTm,ThDyAlign &Al)
{    
    CHitAligned Hit(Al);

    osTm    << endl            
            << "Match:"            << sep<< Hit._mt            << sep
            << "M Match:"        << sep<< Hit._mm            << sep
            << Al._sd->Name()    << sep<< Al._sd->Len()        << sep
            << "Tm sond:"         << sep<< KtoC(Al._sd->_Tm.Ave())<< sep
            << "Iter:"             << sep
            << "gaps sond:"     << sep<< Hit._sgap            << sep
            << Hit._i0            << sep<< Hit._i                << sep
            << (char*)Hit._sd.c_str()

             << endl            
            << Hit._G            << sep<< KtoC(Hit._Tm)        << sep
            << "Long al:"        << sep<< Hit._l                << sep
            << Al._tg->Name()    << sep<< Al._tg->Len()        << sep
            << "Tm target:"     << sep<< KtoC(Al._tg->_Tm.Ave())<< sep
            << Al.IterationNum()<< sep
            << "gaps tg:"         << sep<< Hit._tgap            << sep
            << Hit._j0        << sep<< Hit._j            << sep
            << (char*)Hit._tg.c_str()

;    // << "G aneling: " << sep<< Al.G()         << sep;
        //return osTm;    
}

ofstream &operator<<(ofstream &stream, ThDyAlign_Tm &TmAl) 
{    stream    << endl<< "T melting Align";     //<< (ThDyAlign)TmAl;
    print_ThDyAlign (stream, TmAl);
    return stream;
}
ofstream &operator<<(ofstream &stream, ThDyAlign_G &G_Al) 
{    stream    << endl<< "G TDDP-Align"; 
    print_ThDyAlign (stream, G_Al);
    return stream;
}
ofstream &operator<<(ofstream &stream, FracTDAlign &FrAl) 
{    stream    <<  endl<< "Frac TDDP: " ;     //<< (ThDyAlign)TmAl;
    print_ThDyAlign (stream, FrAl);
    stream    << FrAl.IterationNum();//
    return stream;
}
void print_Tm        (ofstream &osTm, CMultSec    &pr, int MaxGrDeg=-1, char sep[]=";" )
//void print_Tm (ofstream &osTm, CMultSec    &pr, int MaxGrDeg, char sep[])
{    
    for ( auto& CurSec : pr.SecL() ) 
    {    
        CSec *s = CurSec.get() ;
        if (MaxGrDeg!=-1 && MaxGrDeg < s->Degeneracy() ) continue;
        pr.AddMultiSec (  s->CreateNonDegSet () ) ;
    
        //assert (( (osTm  << endl<< s->_c                 << sep//     << "Tm=" <<'\b'    
        //                 << (s->_Tm.Min() - 273.15) << sep // << " °C"             << " ("  
        //                 << (s->_Tm.Ave() - 273.15)    << sep //  << " °C"  << ") "
        //                 << (s->_Tm.Max() - 273.15) << sep //  << " °C"  
        //                 << s->Name()              << endl //    
        //                                                //     << "\n" 
        //            ) , 1 ) ) ;
    }
}

void    ThDyAlign::Export_DPMz_Tm(ofstream &osDP_mz, char *sep)
{    
    register long  i_1,  i,   j_1,  j ;        // las posiciones
    register Base  a_1,  a,   b_1,  b ;        // las bases en esas posiciones ?? exactamente como ?? de la pos anterior?
//    register float S0, S1, S2, H0, H1, H2 ;
    long ls = _sd->Len() + 2 , lt = _tg->Len() + 2 ;    // ponerlo como miembro ??  incluye 2x'$'
    float forbidden_enthalpy = _NNpar->forbidden_enthalpy,
          forbidden_entropy  = _NNpar->forbidden_entropy ;        // ponerlo como miembro ??


        osDP_mz    << endl << "forb_H:" <<sep <<forbidden_enthalpy<<sep <<sep << "max (Tm) matriz:" 
            << endl << "forb_S:"  <<sep <<forbidden_entropy
            << endl << "Init_H:"  <<sep <<_NNpar->GetInitialEntropy()            
            << endl ;

        a=0; a_1 = a ;i=0;
        osDP_mz<<endl << "j" <<sep<< "b_1" <<sep<< "b"<<sep;
        osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        }

        a=0;  j=0 ;i=0; b_1 = b=0 ;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << KtoC(_NNpar->CalcTM (dS0 (i, j) , dH0 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS1 (i, j) , dH1 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS2 (i, j) , dH2 (i, j) ))<<sep;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
        osDP_mz << KtoC(_NNpar->CalcTM (dS0 (i, j) , dH0 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS1 (i, j) , dH1 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS2 (i, j) , dH2 (i, j) ))<<sep;
        }
    
        
    b=0 ;
    for (j = 1 ; j<= lt ; j++ )
    {    b_1 = b ;                  j_1 = j-1 ;    // se justifica esto ????
        b   = bkn2c_nu    [_tg->_b [j_1]] ;        // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
        a=0; i=0;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << KtoC(_NNpar->CalcTM (dS0 (i, j) , dH0 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS1 (i, j) , dH1 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS2 (i, j) , dH2 (i, j) ))<<sep;

        for (i = 1 ; i<= ls ; i++ )
        {    //RestHS (i);
            a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
        osDP_mz << KtoC(_NNpar->CalcTM (dS0 (i, j) , dH0 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS1 (i, j) , dH1 (i, j) ))<<sep
                << KtoC(_NNpar->CalcTM (dS2 (i, j) , dH2 (i, j) ))<<sep;
        }
    }
}
void    ThDyAlign::Export_DPMz_S(ofstream &osDP_mz, char *sep)
{    
    register long  i_1,  i,   j_1,  j ;        // las posiciones
    register Base  a_1,  a,   b_1,  b ;        // las bases en esas posiciones ?? exactamente como ?? de la pos anterior?
//    register float S0, S1, S2, H0, H1, H2 ;
    long ls = _sd->Len() + 2 , lt = _tg->Len() + 2 ;    // ponerlo como miembro ??  incluye 2x'$'
    float forbidden_enthalpy = _NNpar->forbidden_enthalpy,
          forbidden_entropy  = _NNpar->forbidden_entropy ;        // ponerlo como miembro ??
    osDP_mz    << endl << "forb_H:" <<sep <<forbidden_enthalpy<<sep <<sep << "S matriz:" 
            << endl << "forb_S:"  <<sep <<forbidden_entropy
            << endl << "Init_H:"  <<sep <<_NNpar->GetInitialEntropy()            
            << endl ;

        a=0; a_1 = a ;i=0;
        osDP_mz<<endl << "j" <<sep<< "b_1" <<sep<< "b"<<sep;
        osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        }

        a=0;  j=0 ;i=0; b_1 = b=0 ;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << dS0 (i, j) <<sep<< dS1 (i, j) <<sep<< dS2 (i, j)<<sep;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << dS0 (i, j) <<sep<< dS1 (i, j) <<sep<< dS2 (i, j)<<sep;
        }
    
        
    b=0 ;
    for (j = 1 ; j<= lt ; j++ )
    {    b_1 = b ;                  j_1 = j-1 ;    // se justifica esto ????
        b   = bkn2c_nu    [_tg->_b [j_1]] ;        // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
        a=0; i=0;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << dS0 (i, j) <<sep<< dS1 (i, j) <<sep<< dS2 (i, j)<<sep;

        for (i = 1 ; i<= ls ; i++ )
        {    //RestHS (i);
            a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << dS0 (i, j) <<sep<< dS1 (i, j) <<sep<< dS2 (i, j)<<sep;
        }
    }
}
void    ThDyAlign::Export_DPMz_H(ofstream &osDP_mz, char *sep)
{    register long  i_1,  i,   j_1,  j ;        // las posiciones
    register Base  a_1,  a,   b_1,  b ;        // las bases en esas posiciones ?? exactamente como ?? de la pos anterior?
//    register float S0, S1, S2, H0, H1, H2 ;
    long ls = _sd->Len() + 2 , lt = _tg->Len() + 2 ;    // ponerlo como miembro ??  incluye 2x'$'
    float forbidden_enthalpy = _NNpar->forbidden_enthalpy,
          forbidden_entropy  = _NNpar->forbidden_entropy ;        // ponerlo como miembro ??

        osDP_mz    << endl << "forb_H:" <<sep <<forbidden_enthalpy<<sep <<sep << "H matriz:" 
            << endl << "forb_S:"  <<sep <<forbidden_entropy
            << endl << "Init_H:"  <<sep <<_NNpar->GetInitialEntropy()            
            << endl ;

        a=0; a_1 = a ;i=0;
        osDP_mz<<endl << "j" <<sep<< "b_1" <<sep<< "b"<<sep;
        osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << i <<sep<< basek[a_1] <<sep<< basek[a]<<sep ;
        }

        a=0;  j=0 ;i=0; b_1 = b=0 ;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << dH0 (i, j) <<sep<< dH1 (i, j) <<sep<< dH2 (i, j)<<sep;
        for (i = 1 ; i<= ls ; i++ )
        {    a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << dH0 (i, j) <<sep<< dH1 (i, j) <<sep<< dH2 (i, j)<<sep;
        }
    
        
    b=0 ;
    for (j = 1 ; j<= lt ; j++ )
    {    b_1 = b ;                  j_1 = j-1 ;    // se justifica esto ????
        b   = bkn2c_nu    [_tg->_b [j_1]] ;        // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
        a=0; i=0;
        osDP_mz<<endl << j <<sep<< basek_c[b_1] <<sep<< basek_c[b]<<sep;
        osDP_mz << dH0 (i, j) <<sep<< dH1 (i, j) <<sep<< dH2 (i, j)<<sep;

        for (i = 1 ; i<= ls ; i++ )
        {    //RestHS (i);
            a_1 = a ;     i_1 = i-1 ;    // se justifica esto ????
            a = _sd->_b [i_1] ;
            osDP_mz << dH0 (i, j) <<sep<< dH1 (i, j) <<sep<< dH2 (i, j)<<sep;
        }
    }
}
void    ThDyAlign::Export_DPMz_Pre(ofstream &osDP_mz)
{    register long  i,  j ;        //  i y j en la DPMz usan las pos i-2 y j-2 en la sec !!!!!! 
    register Base  a,   b ;        // las bases en esas posiciones -- exactamente como -- de la pos anterior
    long ls = _sd->Len() + 2 , lt = _tg->Len() + 2 ;    // lo mismo que _LenSond= LenSond +2 ;// Para que incluya los ´$´ ? y _LenTarg= LenTarg +2 ;
        
    osDP_mz<<endl<<"-"<<sep ;
    for ( i  = 0 ; i<= ls ; i++ )
        {    a = _sd->_b [i] ;
            osDP_mz << basek[a] <<sep ;
        }
    b   = bkn2c_nu    [_tg->_b [0]] ;
    osDP_mz<<endl<<basek_c[b]<<sep ;
    for ( i  = 0 ; i<= ls ; i++ )
        osDP_mz <<  i <<sep ;

    for (j = 1 ; j<= lt ; j++ )
    {    b   = bkn2c_nu    [_tg->_b [j]] ;        // se busca la sonda en la sec complementaria a la target (la target se considerara de doble cadena) 
        osDP_mz<<endl <<  basek_c[b]<<sep <<j <<sep;
        for (i = 1 ; i<= ls ; i++ )
            osDP_mz << step (i, j) <<sep;
        
    }
}

//_pre0 = new Step  [TableSize];
    //_pre1 = new Step  [TableSize];
    //_pre2 = new Step  [TableSize];
        // HACE FALTA ??????????????????
    //memset(_dH1,0,_TableSize*sizeof (float));
    //memset(_dH2,0,_TableSize*sizeof (float));
    //memset(_dS0,0,_TableSize*sizeof (float));
    //memset(_dS1,0,_TableSize*sizeof (float));
    //memset(_dS2,0,_TableSize*sizeof (float));
    //memset(_pre,0,_TableSize*sizeof (float));
    //memset(_pre0,0,TableSize*sizeof (Step));
    //memset(_pre1,0,TableSize*sizeof (Step));
    //memset(_pre2,0,TableSize*sizeof (Step));


        //if(f_fas)
        //    sCand->ExportSondenFasta(osFasta,fileName);
        //if(f_fas)
        //    sCand->ExportSondenFasta(osFasta,fileName);

        //assert( (cout    <<endl<<"Comparing: "<<sep<<    newtg->_Sec._name    <<sep<<    newtg->_Sec._len    <<sep<<    newtg->_NumPosCand    <<sep<<    newtg->_NumCand
        //                <<sep  <<       tg->_Sec._name    <<sep<<       tg->_Sec._len    <<sep<<       tg->_NumPosCand    <<sep<<       tg->_NumCand
        //        ,1));

        //_TDATmC->FindCommon(newtg, tg );

        //assert( (cout    <<endl<<"Result: "<<sep<<    newtg->_Sec._name    <<sep<<    newtg->_Sec._len    <<sep<<    newtg->_NumPosCand    <<sep<<    newtg->_NumCand
        //                <<sep  <<       tg->_Sec._name    <<sep<<       tg->_Sec._len    <<sep<<       tg->_NumPosCand    <<sep<<       tg->_NumCand
        //                <<sep    << "Hits T: "            <<sep<<        _TDATmC->_THits
        //                <<sep    << "HitsOK : "            <<sep<<        _TDATmC->_HitsOK
        //        ,1));


//ofstream &operator<<(ofstream &osTm, ThDyAlign &Al) 


//float    ThDyAlign::Get_G    (long i, long j, float Tm )const //float GAlign::GetFreeEnergyK(int i, int j, float t)
//{    return -(Get_H(i,j)- Tm * Get_S(i,j) );   // dG = dH - T * dS 
//}                                              // first, need to determine optimum values for dG and dH


//CHit (long i, long j, long i0, long j0, float H, float S,float max, ThDyAlign::Step st)




    //long            _i,_j, _i0, _j0, _l;
    //DNAstrand        _strnd;
    //ThDyAlign::Step _Step ;
    //float            _H, _S, _G, _Tm, _max ;








