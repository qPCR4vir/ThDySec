/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file ThDySec\src\ThDy_programs\ThDy_SondeDesign.cpp
*
* @brief
*
*/

//#include "StdAfx.h"
#pragma unmanaged

#include "ThDySec\th_dy_align.h"
#include "ThDy_programs/prog_comm_functions.h"


/// add and ThDy compare local target sequences into canditates and recursively add the internal groups of targets 
void FindSonden( CMultSec               *tg,             ///< the target sequences          /*int& tgN,*/ 
                 int&                   compN,           ///< a comparition consecutive number
                 CMSecCand&             msCand,          ///< the current ThDy comparition structure
                 CProgParam_SondeDesign *IPrgPar_SdDes   ///< 
               )
{
    CProgParam_SondeDesign::targets_comp res;

    for ( auto &CurSec : tg->SecL() )  // recorre todos los        ----  targets  ------
    {    CSec &nt = *CurSec ;  // our new target
        
        if ( nt.Degeneracy() > 1 )  continue ;           // No analiza las target deg...por ahora.Facil de ampliar
        if ( ! nt.Selected() )      continue ;           // No analiza las target non selected

        res.target_1_name = nt.Name();

        /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // std::cout<<"\n"<<nt.Name();
        /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        auto ntg = msCand.AddBeging(nt);         // create current ThDy candidate target
        CSecCand &newtg = *ntg;         
        for( auto& curTg : msCand._LSecCand )    //  y lo comp con todos los anteriormente anadidos
        {    
            CSecCand &curtg = *curTg;
            
            res.iteration_num = ++compN;
            res.target_num    = msCand._NSecCand ;
            res.target_2_name = curtg._Sec.Name();
            res.before        = {  msCand._TNumPosCand, msCand._TNumCand
                                 , newtg._NumPosCand,  newtg._NumCand
                                 , curtg._NumPosCand,  curtg._NumCand };

            msCand.FindCommon    ( newtg, curtg, IPrgPar_SdDes->_design )    ;

            res.after         = {  msCand._TNumPosCand, msCand._TNumCand
                                 , newtg._NumPosCand,  newtg._NumCand
                                 , curtg._NumPosCand,  curtg._NumCand };

            res.THits         =   msCand._TDATmC->_THits ;
            res.HitsOK        =      msCand._TDATmC->_HitsOK  ;

            IPrgPar_SdDes->targets_comparitions.push_back( res  )    ;    
        }
        msCand._LSecCand.emplace_back(ntg.release());    // add current target to the end of targets
    }
    for (auto &CurMSec : tg->MSecL() )  // recorre todos los subgroups of targets
    {    
        FindSonden(CurMSec.get(), compN, msCand, IPrgPar_SdDes );
    }
}
std::vector<std::string> CProgParam_SondeDesign::targets_comp::headers
  { 
       "Num T Pos", 
       "Num T Cand", 
         
       "Num Pos", 
       "Num Cand", 
       "Targ Num", 
       "Targ name", 
       "Num Cand", 
       "Num Pos", 
         
       "Num T Hits", 
       "Num Hits OK", 
         
       "Num Pos", 
       "Num Cand", 
       "Iterat#", 
       "Targ name", 
       "Num Cand", 
       "Num Pos", 
         
       "Num T Pos", 
       "Num T Cand" 
};


void write_results(CProgParam_SondeDesign *IPrgPar_SdDes)
{
    IPrgPar_SdDes->_cp.Check_NNp_Targets();   /// \todo review

    std::ofstream osNCand(IPrgPar_SdDes->_cp._OutputFile.get() + ".TgCand.csv");
    osNCand.precision(2);
    osNCand << "\n" << "Num T Pos" << sep << "Num T Cand"         ///\todo use CProgParam_SondeDesign::targets_comp::headers
        << sep << "Targ Num" << sep << "Iterat#"
        << sep << "Targ name" << sep << "Num Pos" << sep << "Num Cand"
        << sep << "Targ name" << sep << "Num Pos" << sep << "Num Cand"
        << sep << "Num T Hits" << sep << "Num Hits OK"
        << sep << "Targ name" << sep << "Num Pos" << sep << "Num Cand"
        << sep << "Targ name" << sep << "Num Pos" << sep << "Num Cand"
        << sep << "Num T Pos" << sep << "Num T Cand"
        << std::fixed;

    //for( ;;)
    //    osNCand << "\n" << msCand._TNumPosCand << sep << msCand._TNumCand
    //    << sep << msCand._NSecCand << sep << compN
    //    << sep << newtg._Sec.Name() << sep << newtg._NumPosCand/*In*/ << sep << newtg._NumCand/*In*/
    //    << sep << curtg._Sec.Name() << sep << curtg._NumPosCand/*In */ << sep << curtg._NumCand/*In*/;

       //<<sep<< msCand._TDATmC->_THits<< sep<< msCand._TDATmC->_HitsOK 


}

int SondeDesignProg ( CProgParam_SondeDesign &IPrgPar_SdDes)
{    
    
    time_t t_0 = time(nullptr);

    IPrgPar_SdDes._cp.Actualice_NNp();  /// \todo review

                                         /// \todo use only the file name not the path (use filesystem)
    IPrgPar_SdDes.probes = *IPrgPar_SdDes._cp.AddSeqGroup(*IPrgPar_SdDes.probes,
                                                           IPrgPar_SdDes._cp._OutputFile.get()   );


    /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::cout << "Num" << sep << "SecName" << sep << "Inic" << sep << "Fin" << sep << "Len" << sep << "Tm" << sep << "Sec"
    //    << sep << "H" << sep << "S" << sep << "G(Ta=" << KtoC(IPrgPar_SdDes->_cp.Ta.get()) << " gr)" << sep << "No.matchs";
    /// \debug   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    time_t t_sec = time(nullptr);

    CMSecCand    msCand    (    
                   convCtoK_ctok( IPrgPar_SdDes->_sL ) ,
                            CtoK(IPrgPar_SdDes->_Tm_sig),         IPrgPar_SdDes->_G_sig      * 1000,
                            CtoK(IPrgPar_SdDes->_MaxSd_nTgTm),     IPrgPar_SdDes->_MinSd_nTgG * 1000,    
                            CtoK(IPrgPar_SdDes->_MaxSelfTm),     IPrgPar_SdDes->_MinSelfG   * 1000     );

    if (!IPrgPar_SdDes->_cp._pSeqTargets->_NNPar)   ///  \todo Make seq specific
        IPrgPar_SdDes->_cp._pSeqTargets->_NNPar=IPrgPar_SdDes->_cp._pSaltCorrNNp;

    msCand.Use                 (IPrgPar_SdDes->_cp._pSeqTargets);  ///  \todo count comparisons for "progress"
    msCand._TDATmC->SetTa (CtoK(IPrgPar_SdDes->_cp._Ta         ));

    time_t t_al_created = time(nullptr);

    int /*tgN(0),*/ compN=0;

    FindSonden(IPrgPar_SdDes->_cp._pSeqTargets, /*tgN,*/ compN, msCand, IPrgPar_SdDes );

    NumRang<float> ExtrCovPerc(IPrgPar_SdDes->Coverage.get());
    if (! IPrgPar_SdDes->common.get())         ExtrCovPerc.SetMax(101.0f);
    if (! IPrgPar_SdDes->unique.get())         ExtrCovPerc.SetMin( -1.0f);

    /// Will return probes with a percent of other-target coverage with is not intern to the range ExtrCovPerc.
    /// That is: probes which hybrid in one target but in not than more than in ExtrCovPerc.Min % of the others, 
    /// and additionally, probes with hybrid in one target and at last in ExtrCovPerc.Max % of the others.

    /// filtre what hybrid into nontargets   ...........    !!!!


    msCand.ExportCommonSonden(  IPrgPar_SdDes->_design, 
                                ExtrCovPerc, 
                                IPrgPar_SdDes->probes,
                                IPrgPar_SdDes->_cp._OutputFile.get().c_str(), 
                                fileFormat ( (int)fasta | (int)csv )
                              );

    time_t t_tm_cal = time(nullptr);
    std::cout<< "\n" << "\n" <<"Time sec= "            << sep<< t_sec            - t_0        
                     << "\n" <<"Time Ob crea="        << sep<< t_al_created    - t_sec        
                     << "\n" <<"Time Tm calc= "        << sep<< t_tm_cal        -t_al_created ;
    return 1;
}




