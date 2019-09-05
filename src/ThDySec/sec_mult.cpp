/**
* Copyright (C) 2009-2019, Ariel Vina Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2019
*
* @file  ThDySec/src/ThDySec/sec_mult.cpp
*
* @brief 
*/

//#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <math.h>
#include <list>
#include <stack>
#include <algorithm>

#include <ThDySec/sec_mult.h>
#include <ThDySec/common.h>

using namespace DegCod;
namespace fs = std::filesystem;

/// \todo make more efficient and elegant
fs::path unique_filename(fs::path name)
{
    while(fs::exists(name))
    {
        std::string ext = name.extension().string();
        std::string nam = name.stem().string();
        //std::string pat = name.remove_filename() ;
        name = name.remove_filename() / fs::path( nam + "X" + ext) ;
    }
    return name;
}

/// export all seq in this MSec in a file named as this tree-path.
bool    CMultSec::Export_from   ( CMultSec& tree_base, bool only_selected)
{
    fs::path dir, file;
    auto my_tree_path= path();
    std::cout << "my_tree_path: " + my_tree_path << "\n";

    for ( const auto& CurMSec : tree_base.MSecL())           // recorre todos las msec to find what is my tree_base dir
    {    
        auto major_sub_tree_path= CurMSec->path();
        std::cout << "major_sub_tree_path: " + major_sub_tree_path << "\n";       // like: 'All seq\Target seq\'
        std::cout << "major_sub_tree_dir: " + CurMSec->_orig_file.string() << "\n";   // like: '../ThDy/sequences/'

                                                             // found OK only if my_tree_path begin with major_sub_tree_path
        if ( major_sub_tree_path.empty() || my_tree_path.find(major_sub_tree_path))
            continue;

                                                             // todo: find a robust solution for this hack: replace tree base with the original directory base.
        file = dir = my_tree_path.replace(0, major_sub_tree_path.length()-1, CurMSec->_orig_file.string());
        std::cout << "file: " + file.generic_string() << "\n";
        file = file.parent_path().replace_extension("export.fasta");
        std::cout << "file': " + file.generic_string() << "\n";
        dir = file.parent_path();
        std::cout << "dir': " + dir.generic_string() << "\n";

        //if (!fs::is_regular_file(dir))
        fs::create_directories(dir);
        Export_as(unique_filename(file).string(), only_selected);
        return true;
    }
    return false;
}

bool      CMultSec::Export_local_seq   ( CMultSec& base, bool only_selected)
{
    //assert(ms);
    fs::path dir, file;
    auto s= path();
    auto b= base.path();
    if ( s.find(b))  return false;            // finded OK only if s beging with b
    file = dir = s.replace(0, b.length(), base._orig_file.string());
    file.replace_extension("fasta");
    dir.remove_filename();

    fs::create_directories(dir);
    Export_as(unique_filename(file).string(), only_selected);
    return true;
}

CMultSec::CMultSec (const fs::path &path,
                    std::shared_ptr<CSaltCorrNN>  NNpar    , 
                    bool           all_dir  /*= false*/,
                    float          MaxTgId    /*= 100*/,
                    LonSecPosRang  SecLim    /*= LonSecPosRang {1,0}*/,     
                    SecPosRang     SecLenLim/*= SecPosRang    {0,0}*/,
                    bool           loadSec  /*= true*/
                 )  
    :
        _SecLim     (SecLim),
        _MaxTgId    (MaxTgId), 
        _SecLenLim  (SecLenLim),
        _NNPar      (NNpar)/*,
        _orig_file_path       (file)*/
{
    fs::path  itf(path);

    if (all_dir)                    // Load all files and directories recursiverly?
    {
        if (fs::is_regular_file(itf))
            itf = itf.parent_path();

        _name = itf.filename ().string();     // The new MSec take the name of the dir.
        _orig_file = itf;                     // and the _orig_file point to it.

        fs::directory_iterator rdi{ itf }, end;

        for (; rdi != end; ++ rdi)
            AddFreeMultiSec(std::make_shared<CMultSec>(
                    rdi->path(),
                    NNpar,
                    fs::is_directory(rdi->status()),
                    MaxTgId,
                    SecLim,
                    SecLenLim,
                    loadSec));
        return;
    }
    else                           // Load only this file
        if (itf.has_filename())
        {
            _name = itf.filename ().string();     /// The new MSec take the name of the file.
            _orig_file = itf;                 /// and the _Path point directly to the file.
            if (loadSec)
               AddFromFile(itf);  /// will throw if not a file
            return;
        }
    throw std::ios_base::failure(std::string("Non recursive load from a non regular file: no such sequence file: ")
                                + path.string() );
}

int        CMultSec::AddFromFile (const fs::path& file)    // return la cantidad de sec add ------  AddFromFile   ---
{    
    std::ifstream ifile( file );
    if ( ! ifile ) 
        throw std::ios_base::failure("Could not open the sequence file: "+ file.string() );
    return AddFromFile(ifile);      /// \todo: rethrow anadiendo el nombre del file
}

int        CMultSec::AddFromFile (std::ifstream& ifile)        // return la cantidad de sec add -----------  AddFromFile   ---
{        
    //if (  _SecLim.Max() <= _SecLim.Min() )  // _SecLim.SetMax(0) ;  // if ( _SecEnd<=_SecBeg) _SecEnd=0 ;

    int j=0;
    char c1;
    ifile >> std::skipws  >> c1;
    if ( ! ifile.good() )     
    {
        throw std::ios_base::failure(std::string("Could not read the sequence file: ")/*+ file */);
    }

    if( c1 =='>' )                                                                         
        return AddFromFileFASTA (ifile);
                    // esto es FASTA, si no  BLAST o GB o ...
    if (c1 =='<' )
    {
        std::string xml_DOCTYPE ;
                            // <?xml version="1.0"?>
                            // <!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
        if ( ! std::getline (ifile, xml_DOCTYPE,'>')  ||  ! std::getline (ifile, xml_DOCTYPE,'>') )
        {    
            throw std::ios_base::failure(std::string("Could not open the sequence file: ")/*+ file*/ );
        }

        if        (std::string::npos != xml_DOCTYPE.find("DOCTYPE BlastOutput") )
            return AddFromFileBLAST (ifile);
        else if (std::string::npos != xml_DOCTYPE.find("DOCTYPE Bioseq-set" ) )
            return AddFromFileGB (ifile);
        return 0;
    }
    if (c1 =='L' )    // LOCUS ?
    {    
        ifile.putback(c1);                                                
        return AddFromFileGBtxt (ifile);    
    }
    return 0;   // cerr unknow format        
}     

int        CMultSec::AddFromFileFASTA (std::ifstream &ifile)  // -------------------    AddFromFileFASTA   ------------
{    
    LonSecPos secBeg = _SecLim.Min()  ;   // here beginig to read, set to 1 if originaly <1
    if (secBeg < 1) secBeg = 1;

    LonSecPos secEnd = _SecLim.Max()  ;   // here end to read, set to 0 to ignore it, and if < secBeg
    if (secEnd <= secBeg) secEnd = 0;

    LonSecPos lmin = _SecLenLim.Min() ;   // read at least this # of bases, always >=1
    if (lmin < 1 ) lmin = 1;

    LonSecPos lmax = _SecLenLim.Max() ;   // read at most this # of bases, ignore if =0
    if (lmax < lmin ) lmax = 0;

    if (secEnd)
        if (lmax == 0 || lmax > secEnd - secBeg +1)
            lmax = secEnd - secBeg +1;

    if (lmax && lmin > lmax)
        return 0;
    
    int NumSeq = 0;   // number os sequences
    std::string Descriptor  ;
    while (std::getline (ifile, Descriptor) )
    {
        size_t name_end=Descriptor.find_first_not_of(
                                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890!#$()*+-./<=>@[]^_{|}~"    );
        std::string Fasta_NAME = Descriptor.substr(0,name_end );
        Descriptor=Descriptor.substr(Fasta_NAME.length());

        std::string Fasta_SEC ;
        if (!std::getline (ifile, Fasta_SEC,'>'))         break ;
        if (!ifile.eof() && !Fasta_SEC.empty() && Fasta_SEC.back()!= '\n')
            std::cout << "\n Warning: simbol > do not begin the line. Previous seq:" << Fasta_SEC;

        if (     Fasta_SEC.length() < static_cast<std::size_t>( secBeg + lmin -1 )   )    
             continue;

        auto sec = std::make_shared<CSec>(  Fasta_SEC ,
                                            _Local._NSec, 
                                            Fasta_NAME , 
                                            _NNPar,
                                            lmax, 
                                            secBeg  //  \todo: cambiar constr de CSec ???? use make_unique?
                                          ) ;     
                        
        if ( sec->Len() < lmin  )    
            continue; 

        //if ( sec->_aln_fragment)
        //    sec->_aln_fragment->aln.set(*this, sec->_aln_fragment->sq.Min(),
        //                                       sec->_aln_fragment->sq.Max());

        sec->Description(trim_string(Descriptor));

        auto idem=Idem(*sec);

        if (idem != SecL().end()) 
        {
            if ((*idem)->Len() >= sec->Len() ) 
            {
                sec->Selected(false);
                sec->Filtered(true);
                InsertSecAfter (idem, sec) ;    
            }
            else
            {
                (*idem)->Selected(false);
                (*idem)->Filtered(true);
                InsertSec(idem, sec);    // std::move?
            }
        }
        else
            AddSec(sec);   // std::move?

        NumSeq++;        
    }
    return NumSeq;    
}

int        CMultSec::AddFromFileBLAST (std::ifstream &fi) // ----------------  CMultSec::     AddFromFileBLAST  format XML---
{    
    ///\todo adapt to multiquery BLAST

    unsigned int _BlastOutput_query_len ;        // x todos los "hits"
    int             id=0;
    std::string       li ;  //  xml_line

    LonSecPos secBeg = _SecLim.Min()  ;   // here beginig to read, set to 1 if originaly <1
    if (secBeg < 1) secBeg = 1;

    LonSecPos secEnd = _SecLim.Max()  ;   // here end to read, set to 0 to ignore it, and if < secBeg
    if (secEnd <= secBeg) secEnd = 0;

    LonSecPos lmin = _SecLenLim.Min() ;   // read at least this # of bases, always >=1
    if (lmin < 1 ) lmin = 1;

    LonSecPos lmax = _SecLenLim.Max() ;   // read at most this # of bases, ignore if =0
    if (lmax < lmin ) lmax = 0;

    if (secEnd)
        if (lmax == 0 || lmax > secEnd - secBeg +1)
            lmax = secEnd - secBeg +1;

    if (lmax && lmin > lmax)
        return 0;
  

    do  {    if ( !std::getline (fi, li,'>') )  return 0; }   // BLAST format error
    while  (std::string::npos==li.find("BlastOutput_query-len") ); fi>>_BlastOutput_query_len;//  <BlastOutput_query-len>267</BlastOutput_query-len>
    
    do 
    {
        CSecBLASTHit::Info i;

        std::string       sec;                        // para CSec
        std::string       nam;
        long              l=0;
        std::string       clas;

		auto scan = [&](char const*s) { while (getline(fi, li, '>') && li.find(s) == std::string::npos);     };

        scan("Hit_num"      )  ;  fi>>i.Hit_num;                                //  <Hit_num>1</Hit_num>
        scan("Hit_id"       )  ; if(!getline(fi, i.Hit_id, '<') ) return id;    //<Hit_id>gi|84028434|gb|DQ318020.1|</Hit_id>
        scan("Hit_def"      )  ; if(!getline(fi, i.Hit_def,'<') ) return id;    //<Hit_def>Wets NIle virus strain ArB3573/82, complete genome</Hit_def>
        scan("Hit_accession")  ; if(!getline(fi, i.Hit_accession,'<'))return id;//<Hit_def>Wets NIle virus strain ArB3573/82, complete genome</Hit_def>
        scan("Hit_len"      )  ;  fi>>i.Hit_len;                //  <Hit_len>11048</Hit_len>
        scan("Hsp_bit-score")  ;  fi>>i.Hsp_bit_score;          //  <Hsp_bit-score>482.786</Hsp_bit-score>
        scan("Hsp_score"    )  ;  fi>>i.Hsp_score;              //  <Hsp_score>534</Hsp_score>
        scan("Hsp_evalue"   )  ;  fi>>i.Hsp_evalue;             //  <Hsp_evalue>3.71782e-133</Hsp_evalue>
        scan("Hsp_query-from") ;  fi>>i.Hsp_query_from;         //  <Hsp_query-from>1</Hsp_query-from>
        scan("Hsp_query-to" )  ;  fi>>i.Hsp_query_to;           //  <Hsp_query-to>267</Hsp_query-to>
        scan("Hsp_hit-from" )  ;  fi>>i.Hsp_hit_from;           //  <Hsp_hit-from>9043</Hsp_hit-from>
        scan("Hsp_hit-to"   )  ;  fi>>i.Hsp_hit_to;             //  <Hsp_hit-to>9309</Hsp_hit-to>
        scan("Hsp_query-frame");  fi>>i.Hsp_query_frame;        //  <Hsp_query-frame>1</Hsp_query-frame>
        scan("Hsp_identity" )  ;  fi>>i.Hsp_identity;           //  <Hsp_identity>267</Hsp_identity>
        scan("Hsp_positive" )  ;  fi>>i.Hsp_positive;           //  <Hsp_positive>267</Hsp_positive>
        scan("Hsp_gaps"     )  ;  fi>>i.Hsp_gaps;               //  <Hsp_gaps>0</Hsp_gaps>
        scan("Hsp_align-len")  ;  fi>>i.Hsp_align_len;          //  <Hsp_align-len>267</Hsp_align-len>
        scan("Hsp_hseq"     )  ; if(!getline(fi, sec,'<'))         return id;//<Hsp_hseq>TACAACATGATGGGAAAGAGAGAGAAGAAG
        scan("Hsp_midline"  )  ; if(!getline(fi, i.Hsp_midline,'<'))return id;//<Hsp_midline>|||||||||||||||||||||||||||||||||||||||||
                                                            

        LonSecPos  secHitBeg {secBeg - i.Hsp_query_from +1};  // omitir estas primeras bases del Hit (de la parte del hit que tenemos)
        LonSecPos  lmaxHit   {lmax };  // max # de bases a leer del Hit (de la parte del hit que tenemos)

        if (secHitBeg <= 1)
        {
            if(lmaxHit)
                lmaxHit += secHitBeg ;
            if(lmaxHit < lmin)
                continue;
            secHitBeg = 1;
        }
        
        //if ( _Hsp_query_to  < secBeg ) // end even before it need to beging
        //    continue;
        //if ( secEnd && _Hsp_query_from > secEnd  ) // secEnd==0 mean dont cut the seq or beging after it need to end
        //    continue;
        //if (lmax && ((_Hsp_query_to - _Hsp_query_from)>lmax) && (_Hsp_hit_to - _Hsp_hit_from)>lmax)
        //    continue;

        auto secH = std::make_shared<CSecBLASTHit>
			                      (     std::move( i ), // _BlastOutput_query_len ,   ??????????
                                        std::move(sec)        ,
                                        lmaxHit,
                                        secHitBeg,          // _SecBeg, _SecEnd,
                                        id,                 //Hit_num   ???        //    char    *    nam,    Hit_def
                                        _NNPar           /*,  //long l=0,(Hit_len ---> NO ) !!!  -->_Hsp_align_len -OK clas,    conc*/
                                  );

        if ( secH->Len() < lmin  )    
            continue;

        long MaxIdem= long(ceil((secH->info.Hsp_align_len*_MaxTgId)/100.0f));    // max of Id base
        secH->_aln_fragment->aln   .set(*this, secH->info.Hsp_query_from, secH->info.Hsp_query_to);

        if( secH->info.Hsp_identity > MaxIdem )
        {
            secH->Selected(false);
            secH->Filtered(true);
            AddSec(secH);
        }
        else
        {
            auto idem=Idem(*secH);
            if (idem != SecL().end())
            {
                if ((*idem)->Len() >= secH->Len() )
                {
                    secH->Selected(false);
                    secH->Filtered(true);
                    InsertSecAfter(idem, secH);
                }
                else
                {
                    (*idem)->Selected(false);
                    (*idem)->Filtered(true);
                    InsertSec(idem, secH);
                }
            }
            else
               AddSec(secH);
        }
        id++;
    }
    while (fi.good() ); 
    return id; 
}

int        CMultSec::AddFromFileGBtxt (std::ifstream &fi) // ----------------  CMultSec::            AddFromFileGBtxt  -----------------------------
{    
    /// \todo update !!

    std::string li ;
    int         id=0;
    auto        all =std::numeric_limits<std::streamsize>::max();

    LonSecPos secBeg = _SecLim.Min()  ;   // here beginig to read, set to 1 if originaly <1
    if (secBeg < 1) secBeg = 1;

    LonSecPos secEnd = _SecLim.Max()  ;   // here end to read, set to 0 to ignore it, and if < secBeg
    if (secEnd <= secBeg) secEnd = 0;

    LonSecPos lmin = _SecLenLim.Min() ;   // read at least this # of bases, always >=1
    if (lmin < 1 ) lmin = 1;

    LonSecPos lmax = _SecLenLim.Max() ;   // read at most this # of bases, ignore if =0
    if (lmax < lmin ) lmax = 0;

    if (secEnd)
        if (lmax == 0 || lmax > secEnd - secBeg +1)
            lmax = secEnd - secBeg +1;

    if (lmax && lmin > lmax)
        return 0;

    do {
        CSecGBtxt::Info inf;
        std::string ORIGIN ;
        std::string sec ;

        /*auto scan = [&](std::string& v, char const*s)
                        { do {  if (!(fi >> v)) return id;
                                fi.ignore(all, '\n');
                             } while (v!=s);
                        };*/
        // LOCUS       X14383                  6875 bp    RNA     linear   VRL 30-JUL-1991
        do {
            if (!(fi >> inf.LOCUS)) return id;   if (inf.LOCUS=="LOCUS") break;
            fi.ignore(all, '\n');
        } while (true);
        if (!(fi>>inf.LOCUS)) return id;     std::cout<< "\nLOCUS: "<<inf.LOCUS<<std::endl;
        if (!(fi>>inf.Seq_inst_length)) continue;         // log error ?
              fi.ignore(all, '\n');

        // Bunyamwera virus L protein RNA, complete cds.
        if (!(fi>>inf.DEFINITION)) return id;       if ( inf.DEFINITION!="DEFINITION" ) continue; // blank before ? log error ?
        if (!std::getline(fi, inf.DEFINITION)) return id;  std::cout<< "DEFINITION: "<<inf.DEFINITION<<std::endl;

        // ACCESSION   X14383
        do {
            if (!(fi >> inf.ACCESSION)) return id;   if (inf.ACCESSION=="ACCESSION") break;
            inf.DEFINITION += " " + inf.ACCESSION;
            if (!std::getline(fi, inf.ACCESSION)) return id;
            inf.DEFINITION += " " + inf.ACCESSION;
        } while (true);
        if (!(fi>>inf.ACCESSION)) return id;             std::cout<< "ACCESSION: "<<inf.ACCESSION<<std::endl;

        // ORGANISM  Bunyamwera virus
        do {
            if (!(fi >> inf.ORGANISM)) return id;   if (inf.ORGANISM=="ORGANISM") break;
            fi.ignore(all, '\n');
        } while (true);
        if (!std::getline(fi, inf.ORGANISM)) return id;               std::cout<< "ORGANISM: "<<inf.ORGANISM<<std::endl;

        // ORIGIN
        do {
            if (!(fi >> ORIGIN)) return id;
            fi.ignore(all, '\n');
        } while (ORIGIN!="ORIGIN");  std::cout<< "ORIGIN: "<<ORIGIN<<std::endl;
        if (!getline(fi, sec, '/')) return id;
        fi.ignore(all, '\n');     std::cout<< "sec: "<<sec.length()<<std::endl;

        auto secGBtxt = std::make_shared<CSecGBtxt>(std::move(inf),
                                                    std::move(sec),
                                                    id,           //    char        *    nam,    DEFINITION    ,
                                                    _NNPar);

        std::cout<<"Created. Degeneracy: "<<secGBtxt->Degeneracy() << std::endl;

        if (secGBtxt->Len() < lmin)
            continue;

        auto idem = Idem(*secGBtxt);

        std::cout<<"Idem check: "<< std::endl;


        if (idem != SecL().end())
        {
            if ((*idem)->Len() >= secGBtxt->Len())
            {
                secGBtxt->Selected(false);
                secGBtxt->Filtered(true);
                InsertSecAfter(idem, secGBtxt);
            }
            else
            {
                (*idem)->Selected(false);
                (*idem)->Filtered(true);
                InsertSec(idem, secGBtxt);
            }
        }
        else
            AddSec(secGBtxt);

        id++;        
    }
    while (fi.good() );
    return id; 
}

int        CMultSec::AddFromFileGB (std::ifstream &ifile)  // ----------------  CMultSec::            AddFromFileGB  -----------------------------
{
    /// \todo update !!

    int        id=0;
    std::string xml_line ;

    do {    char        *    _Textseq_id_accession    =0 ;    
            char        *    _Org_ref_taxname        =0 ;
            char        *    _Seqdesc_title            =0 ;
            long            _Seq_inst_length        =0 ;    
            // para CSec
            char        *sec=0;            //char        *nam=0;        //long         l=0;        //char        *clas=0;

            do  {    getline (ifile, xml_line,'>') ;    if ( ! ifile.good() ) return id; }     // <Textseq-id_accession>DQ318020</Textseq-id_accession>
            while  (std::string::npos==xml_line.find(    "Textseq-id_accession"    ) );
            getline (ifile, xml_line,'<') ;                _Textseq_id_accession=new char[xml_line.length()+1] ;
            xml_line.copy(_Textseq_id_accession,xml_line.length()) ;    _Textseq_id_accession    [xml_line.length()]=0;

            do  {    getline (ifile, xml_line,'>') ;    if ( ! ifile.good() ) return id; } // <Org-ref_taxname>West Nile virus</Org-ref_taxname>
            while  (std::string::npos==xml_line.find(    "Org-ref_taxname"    ) );
            getline (ifile, xml_line,'<') ;                                _Org_ref_taxname    =new char[xml_line.length()+1] ;
            xml_line.copy(  _Org_ref_taxname  ,xml_line.length()) ;        _Org_ref_taxname    [xml_line.length()]=0;
        // <Dbtag_db>taxon</Dbtag_db>
        // <Object-id_id>11082</Object-id_id>
        // <OrgName_name_virus>West Nile virus</OrgName_name_virus>
        // <OrgMod_subtype value="strain">2</OrgMod_subtype>
        // <OrgMod_subname>ArB3573/82</OrgMod_subname>
        // <OrgMod_subtype value="gb-acronym">32</OrgMod_subtype>
        // <OrgMod_subname>WNV</OrgMod_subname>
        // <OrgName_lineage>Viruses; ssRNA positive-strand viruses, no DNA stage; Flaviviridae; Flavivirus; Japanese encephalitis virus group</OrgName_lineage>
        // <SubSource_subtype value="country">23</SubSource_subtype>
        // <SubSource_name>Central African Republic</SubSource_name>
        // <SubSource_subtype value="other">255</SubSource_subtype>
        // <SubSource_name>lineage 2; SMB pass 4, C6/36 pass 1</SubSource_name>
            do  {    getline (ifile, xml_line,'>') ;    if ( ! ifile.good() ) return id; }        // <Seqdesc_title>Wets NIle virus strain ArB3573/82, complete genome</Seqdesc_title>
            while  (std::string::npos==xml_line.find(    "Seqdesc_title"    ) );
            getline (ifile, xml_line,'<') ;                                _Seqdesc_title    =new char[xml_line.length()+1] ;
            xml_line.copy(  _Seqdesc_title  ,xml_line.length()) ;        _Seqdesc_title    [xml_line.length()]=0;
        // <MolInfo_biomol value="mRNA">3</MolInfo_biomol>
        // <MolInfo_completeness value="complete">1</MolInfo_completeness>
        // <Date-std_year>2006</Date-std_year>
        // <Date-std_month>1</Date-std_month>
        // <Date-std_day>1</Date-std_day>
        // <Seq-inst_mol value="rna"/>

            do  {    std::getline (ifile, xml_line,'>') ;    if ( ! ifile.good() ) return id; }  // GB format error
            while  (std::string::npos==xml_line.find("Seq-inst_length") ) ; ifile>>_Seq_inst_length;        // <Seq-inst_length>11048</Seq-inst_length>
            
        // <Seq-inst_strand value="ss"/>

            do  {    getline (ifile, xml_line,'>') ;    if ( ! ifile.good() ) return id; }    // <IUPACna>AGTAGTTCGCCTGTGTGAGCTGACA.... GGTGCTAGAACACAGGATCT</IUPACna>
            while  (std::string::npos==xml_line.find(  "IUPACna"  ) );
            getline (ifile, xml_line,'<') ;         sec=new char[xml_line.length()+1] ;
            xml_line.copy(sec,xml_line.length()) ;    sec[xml_line.length()]=0;    


            auto secGB = std::make_shared<CSecGB>(      
				                            _Textseq_id_accession,
                                            _Org_ref_taxname    ,
                                            _Seqdesc_title,
                                            _Seq_inst_length     ,
                                            sec    ,    
                                            id,            //    char        *    nam,        Hit_def
                                            _NNPar/*,      //    long            l=0,        (Hit_len ---> NO ) !!!  -->_Hsp_align_len -OK clas, conc*/
                                        );
            delete []sec ;

                if ( secGB->Len() >= _SecLenLim.Min()  )        
                {    
                    auto idem =Idem(*secGB);
                    if (idem != SecL().end())
                    {
                        if ((*idem)->Len() >= secGB->Len())
                        {
                            secGB->Selected(false);
                            secGB->Filtered(true);
                            InsertSecAfter(idem, secGB);
                        }
                        else
                        {
                            (*idem)->Selected(false);
                            (*idem)->Filtered(true);
                            InsertSec(idem, secGB);
                        }
                    }
                    else
                        AddSec(secGB);
                    id++;        
                }
        }
    while (ifile.good() ); 
    return id; 
}
    
int        CMultSec::AddFromFileODT (std::ifstream &ifileODT){return 0;}
int        CMultSec::AddFromFileODS (std::ifstream &ifileODS){return 0;}

CMultSec::SecIt CMultSec::Idem ( CSec &sec ) const  // ------  CMultSec:: NotIdem  --- busqueda trivial de sec identicas -------------
{    
    if ( _MaxTgId >= 100  ) //  no restriction on similarity
        return _LSec.end() ;    

    long LenCandSec=sec.Len() ;     // Lenght of the Candidate Sec (to be in the list, with MaxId)

    long MaxErCS= long(ceil(float(LenCandSec*(100.0f-_MaxTgId) ) / 100.0f)); // min of not Id base to be in the list
    
    for (SecIt CurSec = _LSec.begin(); CurSec != _LSec.end(); ++CurSec) // recorre todos las primeras sec de esta misma ms
    {    
        CSec &s = **CurSec ; 
        if (s.Filtered())
            continue;

        // sec:q           -------------------------------------------------------------------------
        //          ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        // s  :h    -------------------------------------------------------------------
        // aln ------------------------------------------------------------------------------------------------
        //
                                       // b - beging, e - end
        long qb=0, qe=LenCandSec-1;    // q - (sec, CandSec) query sec: is this sec alrready here?
        long hb=0, he=s.Len()-1;       // h - (s) hit sec: this sec from the list hit the query?
        if(sec._aln_fragment && s._aln_fragment)
        {
            NumRang<LonSecPos> *SecRang{}, *SRang{};
            //if (sec._aln_fragment->aln.is_comparable_to( s._aln_fragment->aln ) )

            if (sec._aln_fragment->aln.sq && sec._aln_fragment->aln.sq == s._aln_fragment->aln.sq  )
            {
                SecRang = &sec._aln_fragment->aln ;
                SRang   = &  s._aln_fragment->aln ;
            }else
            if (sec._aln_fragment->consensus.sq && sec._aln_fragment->consensus.sq == s._aln_fragment->consensus.sq  )
            {
                SecRang = &sec._aln_fragment->consensus ;
                SRang   = &  s._aln_fragment->consensus ;
            }else
            if (sec._aln_fragment->bio_ref.sq && sec._aln_fragment->bio_ref.sq == s._aln_fragment->bio_ref.sq  )
            {
                SecRang = &sec._aln_fragment->bio_ref ;
                SRang   = &  s._aln_fragment->bio_ref ;
            }else
            if (sec._aln_fragment->sq_ref.sq && sec._aln_fragment->sq_ref.sq == s._aln_fragment->sq_ref.sq  )
            {
                SecRang = &sec._aln_fragment->sq_ref ;
                SRang   = &  s._aln_fragment->sq_ref ;
            } 

            if (SecRang)
            {
               if (sec._aln_fragment->aln.Min() > s._aln_fragment->aln.Min() )
                   hb = sec._aln_fragment->aln.Min() - s._aln_fragment->aln.Min() ;
               else
                   qb = s._aln_fragment->aln.Min() - sec._aln_fragment->aln.Min() ;

               if (sec._aln_fragment->aln.Max() > s._aln_fragment->aln.Max() )
                   he -= (sec._aln_fragment->aln.Max() - s._aln_fragment->aln.Max()) ;
               else
                   qe -= s._aln_fragment->aln.Max() - sec._aln_fragment->aln.Max() ;
            }
        }
        long l, MaxEr ;  // active lenght, MaxEr - if we detect this number of error the sec are suficient different and are not filtred
        if (s.Len() < LenCandSec)         
        {    
            l=s.Len() ;            
            MaxEr= long(ceil(l*(100-_MaxTgId)  / 100));    
        } 
        else                     
        {    
            l=LenCandSec ;                
            MaxEr= MaxErCS;                        
        }
        l=std::min(he-hb, qe-qb);
        long Er=0 ;
        for (long i=0   ; i <=l; ++i ) 
        {
            Er += ( s[i+hb] != sec[i+qb] ) ;
            if (   Er >  MaxEr   ) 
                break ;
        }
        if (   Er <=  MaxEr   ) 
             return CurSec;
    }
    return _LSec.end();
}

CMultSec::pSec CMultSec::AddSec( CMultSec::pSec sec)
{
    if (!sec) return sec ;
    _LSec.push_back(sec);
    UpdateTotalsMoving ( *sec );
    return sec;
}

CMultSec::pSec CMultSec::InsertSec(CMultSec::SecIt pos, CMultSec::pSec sec)
{
    if (!sec) return sec ;
    _LSec.insert(pos, sec);
    UpdateTotalsMoving ( *sec );
    return sec;
}
CMultSec::pSec CMultSec::InsertSecAfter(CMultSec::SecIt preSec, CMultSec::pSec sec)
{    
    return InsertSec(++preSec, sec);
}

void    CMultSec::UpdateTotalsMoving ( CSec &sec ) 
{    
    if (sec._parentMS == this)                                        // no hay sec o ya estaba aqui
        return;
    CMultSec *parMS    = sec._parentMS;                    // /*._Get()*/
    CMultSec *My_parMS = _parentMS;                    // /*._Get()*/
    bool checkExtr(true) ; 
    CMultSec   *cp;

    if (parMS)
    {
        cp=findComParent( parMS);
        parMS->_Local._NSec--    ;                                // la elimino de la ms orig. 
        //parMS->_Global._NSec -- ;                                // 
        if (parMS->_Local._NSec)
            if (parMS->isLocExtreme(sec))
                RecalExtremes();
            else
                checkExtr=false;
        else
            checkExtr=false;

        for ( /*parMS=parMS->*/_parentMS;  parMS != cp &&  parMS ;  parMS=parMS->_parentMS)        // desde localizacion orig subiendo hasta parent comun
            {
                parMS->_Global._NSec -- ;                                // elimino  de este total.
                if (checkExtr && parMS->_Global._NSec)
                    if (parMS->isGlobExtreme(sec))
                        RecalExtremes();
                    else
                        checkExtr=false;
                else
                    checkExtr=false;
            }
    }else
        cp=nullptr;
    Add2LocalExtreme(sec);
    for (My_parMS ; My_parMS!=cp && My_parMS; My_parMS=My_parMS->_parentMS)  // desde mi hacia arriba hasta el com parent anadiendo
    {
        if (checkExtr)
            My_parMS->Add2GlobalExtreme(sec);
        else
        {
            if (parMS) parMS->_Global._NSec ++;                                // sumo sus s a este total.
        }
    }    
    sec._parentMS = (this) ;                            //* std::weak_ptr<CMultSec> */
}

CMultSec   *CMultSec::findComParent( CMultSec *ms)
{
    if(!ms || ms==this) 
        return ms;
    std::stack<CMultSec*> myTree,        oTree;
    CMultSec            *myPms=this, *oPms=ms;
        myTree.push(myPms);
        oTree.push(oPms);

    do 
    {    myPms=myPms->_parentMS;
        myTree.push(myPms);
    }while (myPms);

    do 
    {    oPms=oPms->_parentMS;
        oTree.push(oPms);
    }while (oPms);

    do 
    {    if(  (oPms = myTree.top())   !=  oTree.top()  )
            return myPms;
        myPms= oPms;
        myTree.pop();
        oTree.pop();
    } while (!myTree.empty() && !oTree.empty() );
    return myPms;
}

CMultSec::MSecIt CMultSec::AddFreeMultiSec(CMultSec::pMSec MultSec)  //--------------------------------------    AddFreeMultiSec    --------------------
{
    MSecIt end = _LMSec.end();
    if (!MultSec) return end;
    MSecIt it = _LMSec.insert(end, MultSec);
    UpdateTotalsMoving ( *MultSec );   // al llamar ya esta la ms movida fisicamente. Falta solo actualizar extremes
    return it;
}
//std::shared_ptr<CMultSec> CMultSec::AddMultiSec (std::shared_ptr<CMultSec> ms )  //--------------------------------------    AddFreeMultiSec    --------------------
//{    if (!ms) return nullptr;    
//    _LMSec.push_back(ms);
//    UpdateTotalsMoving ( ms.get() );   // al llamar ya esta la ms movida fisicamente. Falta solo actualizar extremes
//    return ms.get();
//}
void        CMultSec::UpdateTotalsMoving ( CMultSec &msec ) 
{    
    if (msec._parentMS==this)                    // no hay msec o ya estaba antes en una de mis subtrees inmediatas. 
        return;

    CMultSec *parMS   = msec._parentMS;                    // /*._Get()*/
    CMultSec *My_parMS=      _parentMS;                    // /*._Get()*/
    bool checkExtr(true) ; 
    CMultSec   *cp;
    if (parMS)                                        // no es imprescindible. Anadido solo por claridad de intencion
    {    cp=findComParent( &msec);
        for ( parMS;  parMS != cp  && parMS ;  parMS=parMS->_parentMS)            // desde localizacion orig subiendo hasta parent comun
            {
                parMS->_Global._NSec -= msec._Global._NSec ;            // elimino sus s de este total.
                parMS->_Global._NMSec-= msec._Global._NMSec + 1;        // elimino sus ms de este total.
                parMS->_Local._NMSec--    ;                                // la elimino de esta ms. 
                if (checkExtr && parMS->_Global._NSec)
                    if (parMS->isGlobExtreme(msec))
                        RecalExtremes();
                    else
                        checkExtr=false;
            }
    }else
        cp=nullptr;

    Add2LocalExtreme(msec);
    for (My_parMS ; My_parMS!=cp && My_parMS; My_parMS=My_parMS->_parentMS)  // desde mi hacia arriba hasta el com parent anadiendo
    {
        if (checkExtr)
            My_parMS->Add2GlobalExtreme(msec);
        else
        {
            My_parMS->_Global._NSec += msec._Global._NSec ;            // sumo sus s a este total.
            My_parMS->_Global._NMSec+= msec._Global._NMSec + 1;        // sumo sus ms a este total.
        }
    }
    msec._parentMS = (this) ;                            // std::weak_ptr<CMultSec> 

}

        CMultSec::~CMultSec ()                // funciona bien solo si la lista es "lineal"
{
}    
