/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDySec\sec.h
*
* @brief To manipulate DNA sequences for simple thermodynamic modeling of hybridization.
*
* Classes to manipulate DNA sequences. Adapted exclusively to DNA and specifically for thermodynamic calculations.
*
*
*/

#pragma unmanaged    

#ifndef _SEC_H
#define _SEC_H

#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <string>
#include <memory>
#include <vector>

#include <nana/filesystem/filesystem_ext.hpp>


namespace filesystem = std::experimental::filesystem;


#include "sec_basic.h" 
#include "th_dy_param.h"   
#include "common.h" 
 
///\todo crear (adaptar) clase primer derivada de CSec, con pos y Tm en cada sec Target
 
class CMultSec    ;

                   // ---------------------------------------   CSec    ---------------------------------------------------

/// To manipulate DNA sequences (secuencias) for simple thermodynamic modeling of hybridization.

/// Fundamental class to manipulate DNA sequences. Adapted exclusively to DNA and specifically for thermodynamic calculations.
/// Have the sequence in "letter" or nt format, AND in "code" format to avoid millions of repeated conversions. 
/// Remember most characteristics for fast retrieval: length, CG%, degeneracy, etc.
/// Importantly have a sequence of the values of dH and dS from the beginning to each position,
/// which make trivial and fast to find accumulated dH and dS and Tm of any fragment by a simple difference calculation.
/// When degenerated it can have a NonDegSet with all the CSec with degeneracy 1 - each of the variants in the original
/// It track some group of sequences (a "parent" CMultSec) to which it belongs.
/// It can be filtred and can be selected. Normally only selected sequences are used in calculations,
/// and only if the parent group is also selected.
///
/// Some variables have index base [1] while others have [0] in sec. 
///
///      sec_begin                                                                             
///      0                                                                                     sec_end
///      |    1        exp_begin                        exp_end                                |       ---> sq
///      |    |        |                                |                                      |   
///     Z-----ATCGCGTAGCTAGCTAGCTAGCTGACTTGTCTGGTAGCT--GCTATCTAATGCTGATGCTAGTCGATCGTAGCTGC-x----ZX x?
///      |             |                                |
///      1             orig_beging                      orig_end                                       ---> aln
///                    |                                |
///                    1 -fltr_beging                   fltr_end - not counting internal gaps
///     
/// set sq.sq to some original "experimental" CSec if any
///
/// set aln.sq to the parent CMultSec*
class CSec : public CSecBasInfo    
{public:
        int                           x;            ///<  ????
        TemperatureRang               _Tm ;         // float        _Tm, _minTm, _maxTm ;            
        std::shared_ptr<CSaltCorrNN>  _NNpar ;
        float                   _Conc ;            ///< concentration of this molecule. If equal to all others is set to -1 and NNParam is used instead
        std::vector<Code>       _b=   std::vector<Code>   {n_basek-1};            ///< codified sequence, initially? basek
        std::vector<Entropy>    _SdS= std::vector<Entropy>{ _NNpar->GetInitialEntropy()};        ///< dS acumulada. Calcular Delta S sera solo restar la final menos la inicial    
        std::vector<Energy>     _SdH= std::vector<Energy> {  Energy{} };          // 
        CMultSec               *_parentMS{nullptr}    ;   //std::weak_ptr<CMultSec> _parentMS    ;

    /// Some variables have index base [1] while others have [0] in sec. We need to clean the sec, which can contain non bases letter, like tabs, end of line, blancs, etc, but we take into account "-".
    CSec (  const std::string&  sec, 
            int                 id,
            const std::string&  nam,     
            std::shared_ptr<CSaltCorrNN>  NNpar, 
            LonSecPos           lmax=0,    ///< limita la cant de bases originales a leer despues de las primeras secBeg-1 bases 
            LonSecPos           secBeg=1,  ///< base [1] in sec. The first letter in sec to be read. 
            const std::string&  clas="" , 
            float               conc=-1         );

    /// Dummy sequence without actual sequence, but with NNpar and reserved space
    CSec ( long l, std::shared_ptr<CSaltCorrNN>  NNpar) ;


    //long        Len            ()const        {return _SdS.size();} //
    void         CorrectSaltOwczarzy    () ;
    CMultSec    *CreateNonDegSet        () ;              ///< crea todo el set si no existia, solo si existen bases deg: _NDB>0
    CMultSec    *ForceNonDegSet         ();               ///< lo crea siempre, incluso para =1??
    CSec        *GenerateNonDegVariant  (CSec *s, long pos, Base ndb)   ; ///< recursiva
    CSec        *CopyFirstBases         (long pos)   ;    ///< copia parcialmente hasta la pos
    void         CorrectSalt            () { if ( _NNpar->UseOwczarzy () ) CorrectSaltOwczarzy();};
    CSec *Clone(DNAstrand strnd=DNAstrand::direct     ) const override; /// unique_ptr<ISec> crea una copia muy simple. \todo CUIDADO con copias de CSecBLASTHit y otros derivados
    CSec *Clone( long  InicBase,
                 long  EndBase, 
                 DNAstrand   strnd = DNAstrand::direct) const override;

    //virtual CSec*CreateCopy        (DNAstrand strnd=direct) override;//< crea una copia muy simple. CUIDADO con copias de CSecBLASTHit y otros derivados
    //const char    *Get_charSec            ()const{return (const char*)_c.c_str();}  ///   ???????????

    bool        Selected() const;                 //< User-editable    ???????????????????????????????????????????????????????????????????????????
    bool        Selected(bool select)    
                          { 
                              _selected = select; 
                              return Selected(); 
                           }             //< make protected: ?????????????????????????????????????????????????
    Base        operator()        (int i)const{return _b[i];}
    Temperature  Tm    (long pi, long pf   )const;                               ///< Tm de la sonda con sec desde pi hasta pf, inclusive ambas!! 
    Temperature  Tm    (long pi            )const    {return Tm(pi,Len())   ;}   ///< Tm de la sonda con sec desde pi hasta el final, inclusive ambos!!
    Energy        G    (long pi, long pf, float Ta)const;                        ///< G de la sonda con sec desde pi hasta pf, inclusive ambas!! 
    Energy        G    (long pi, float Ta  )const    {return G(pi,Len(), Ta);}   ///< G de la sonda con sec desde pi hasta el final, inclusive ambos!!
    Energy        G    (float Ta           )const    {return G(1 ,Len(), Ta);}   ///< G de la sonda con sec desde inicio hasta el final, inclusive ambos!!
    Energy        G    (long pi, long pf   )const;                               ///< G de la sonda con sec desde pi hasta pf, inclusive ambas!! 
    Energy        G    (long pi            )const    {return G(pi,Len())    ;}   ///< G de la sonda con sec desde pi hasta el final, inclusive ambos!!
    Energy        G    (                   )const    {return G(1,Len())     ;}   ///< G de la sonda con sec desde inicio hasta el final, inclusive ambos!!

     ~CSec() override  ;    
    virtual bool NotIdem(CSec *sec) {return false;}
};

      //<Hsp_num>1</Hsp_num>
      //<Hsp_bit-score>482.786</Hsp_bit-score>
      //<Hsp_score>534</Hsp_score>
      //<Hsp_evalue>3.71782e-133</Hsp_evalue>
      //<Hsp_query-from>1</Hsp_query-from>
      //<Hsp_query-to>267</Hsp_query-to>
      //<Hsp_hit-from>9043</Hsp_hit-from>
      //<Hsp_hit-to>9309</Hsp_hit-to>
      //<Hsp_query-frame>1</>
      //<Hsp_hit-frame>1</>
      //<Hsp_identity>267</>
      //<Hsp_positive>267</>
      //<Hsp_gaps>0</>
      //<Hsp_align-len>267</>
      //<Hsp_qseq>TACAACATGATGGGAAAGAGAGAGAAGAAGCCTGGAGAGTTCGGCAAGGCTAAAGGCAGCAGAGCCATCTGGTTCATGTGGCTGGGGGCTCGCTTCCTGGAGTTTGAAGCTCTCGGATTCCTCAATGAAGACCACTGGCTGGGTAGGAAGAACTCAGGAGGAGGAGTAGAAGGCTTAGGACTGCAGAAGCTTGGGTACATCTTGAAGGAAGTTGGGACAAAGCCTGGAGGAAAGATTTACGCTGATGATACCGCAGGCTGGGACACA</Hsp_qseq>
      //<Hsp_hseq>TACAACATGATGGGAAAGAGAGAGAAGAAGCCTGGAGAGTTCGGCAAGGCTAAAGGCAGCAGAGCCATCTGGTTCATGTGGCTGGGGGCTCGCTTCCTGGAGTTTGAAGCTCTCGGATTCCTCAATGAAGACCACTGGCTGGGTAGGAAGAACTCAGGAGGAGGAGTAGAAGGCTTAGGACTGCAGAAGCTTGGGTACATCTTGAAGGAAGTTGGGACAAAGCCTGGAGGAAAGATTTACGCTGATGATACCGCAGGCTGGGACACA</Hsp_hseq>
      //<>

class CSecBLASTHit : public CSec // ---------------------------------------   CSecBLASTHit    ------------------------------------------------
{public:
    struct Info
    {
        unsigned int    BlastOutput_query_len ;
        unsigned int    Hit_num=0 ;            // para cada hit
        std::string     Hit_id  ;
        std::string     Hit_def ;                // descriptor ??
        std::string     Hit_accession     ;
        long            Hit_len=0 ;
        float           Hsp_bit_score=0 ;
        unsigned int    Hsp_score=0 ;
        double          Hsp_evalue=0 ;
        LonSecPos       Hsp_query_from=0 ;    // dejar signed or unsigned !!!!????
        LonSecPos       Hsp_query_to=0 ;
        LonSecPos       Hsp_hit_from=0 ;
        LonSecPos       Hsp_hit_to=0 ;
        int             Hsp_query_frame=0 ;
        int             Hsp_hit_frame=0 ;
        LonSecPos       Hsp_identity=0 ; // revisar type --- no sera %  : float??
        LonSecPos       Hsp_positive=0 ;
        LonSecPos       Hsp_gaps=0 ;
        LonSecPos       Hsp_align_len=0 ;
        std::string     Hsp_midline ;
        bool            FormatOK=false ;
    };
    Info info;

    CSecBLASTHit(   CSecBLASTHit::Info  &&info,
                    std::string         &&sec,
                    LonSecPos           lmax,                       //long    SecBeg,long    SecEnd,
                    LonSecPos           secBeg,                     //long    SecBeg,long    SecEnd,
                    int                 id,                         //    Hit_num    char    *    nam,    Hit_def
                    std::shared_ptr<CSaltCorrNN>  NNpar,            //    long  l=0,    Hit_len ------> _Hsp_align_len
                    std::string        clas="", 
                    float            conc=-1
                )  :
                        CSec (  std::move(sec),             // sec
                                id,                         // 
                                std::move(info.Hit_accession),   // name
                                NNpar,                      //  . maxlen . secBeg .
                                lmax, 
                                secBeg, 
                                clas,
                                conc ),
                        info(std::move(info)) /*, _SecLim( SecLim ) */
                                               /*,_SecBeg    (SecBeg),     _SecEnd        (SecEnd)*/
                {
                    if ( this->_aln_fragment)   ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
                    {
                        if(this->_aln_fragment->sq.Max())
                            info.Hsp_query_to    = info.Hsp_query_from + this->_aln_fragment->sq.Max()  -1;
                        else
                            info.Hsp_query_to    = info.Hsp_query_from + this->Len() -1;
                        info.Hsp_query_from  = info.Hsp_query_from + this->_aln_fragment->sq.Min()-1;
                    }
                    else
                    {
                        this->_aln_fragment.reset(new Aligned_fragment);      ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
                        info.Hsp_query_to    = info.Hsp_query_from + secBeg + this->Len() -1;
                        info.Hsp_query_from  = info.Hsp_query_from + secBeg-1;
                    }

                    this->_aln_fragment->sq_ref.Set(       info.Hsp_query_from, info.Hsp_query_to);    ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
                    this->_aln_fragment->aln   .set(*this, info.Hsp_query_from, info.Hsp_query_to);
                    this->_aln_fragment->sq    .Set(       info.Hsp_hit_from,   info.Hsp_hit_to  );

                }



    // para cada hit
    //NumRang<long>    _SecLim;                                     //long    _SecBeg;    //long    _SecEnd;
    std::string    Description ()const    override {return _description.length() ? _description : info.Hit_def ; }
};

class CSecGB : public CSec // ---------------------------------------   CSecGB    ------------------------------------------------
{public:
    std::string     _Textseq_id_accession ;    
    std::string     _Org_ref_taxname      ;
    std::string     _Seqdesc_title        ;
    long            _Seq_inst_length      ;        
    
    CSecGB(         std::string     Textseq_id_accession ,    
                    std::string     Org_ref_taxname     ,
                    std::string     Seqdesc_title    ,
                    long            Seq_inst_length     ,    
                    const char    * sec    ,    
                    int             id,                        //    char        *    nam,    Seqdesc_title    ,    
                    std::shared_ptr<CSaltCorrNN>  NNpar,       //    long            l=0,        Seq_inst_length
                    std::string     clas="", 
                    float           conc=-1
                ):
                CSec (sec, id, Textseq_id_accession, NNpar, Seq_inst_length,1, clas, conc ),// actualizar Beg-End
                _Textseq_id_accession   ( std::move(Textseq_id_accession) ) ,
                _Org_ref_taxname        ( std::move(Org_ref_taxname )) ,
                _Seqdesc_title          ( std::move(Seqdesc_title) ) ,                
                _Seq_inst_length        ( Seq_inst_length )                 {}            
    virtual std::string    Description ()const override    {return _description.length() ? _description : _Seqdesc_title ; }

    virtual ~CSecGB(){    }
};

class CSecGBtxt : public CSec // ---------------------------------------   CSecGBtxt    ------------------------------------------------
{public:
    struct Info
    {
        std::string  LOCUS           ;
        LonSecPos    Seq_inst_length ;
        std::string  DEFINITION      ;
        std::string  ACCESSION       ;
        std::string  ORGANISM        ;            // ORIGIN
    };
    Info info;

    CSecGBtxt(      Info&&             inf,
                    std::string&&      sec,
                    int                id,                        //    char        *    nam,
                    std::shared_ptr<CSaltCorrNN>  NNpar,                  //    long            l=0,        Seq_inst_length
                    std::string        clas="", 
                    float              conc=-1
                ):        CSec (std::move(sec), id, inf.LOCUS, NNpar, 0,1, clas, conc ), // actualizar Beg-End
                          info(std::move(inf))       {}

    virtual std::string    Description ()const override    {return info.DEFINITION.empty() ? _description : info.DEFINITION ; }
    
    virtual ~CSecGBtxt() {                }    
};

#endif

//enum fileFormat {fasta=1 , csv=1<<1, f2=1<<2, f3=1<<3} ; // se pueden "OR" los formatos : OUTPUT !!!!!!

///\todo ?? anadir funcion de compactar cod (eliminar los gap y bases deg?). SdH y S se recalculan.
///\todo ?? anadir funcion para regenerar cod no compactado, recordar estado comp/no comp
