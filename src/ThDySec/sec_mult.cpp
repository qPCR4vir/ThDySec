/**
* Copyright (C) 2009-2016, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\src\ThDySec\sec_mult.cpp
*
* @brief 
*/

#ifdef WINDOWS_FORM_GUI
#include "stdafx.h"
#pragma unmanaged
#endif

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <math.h>
#include <list>
#include <stack>
#include <algorithm>

//using namespace std ; 

#include "ThDySec/sec_mult.h"
#include "ThDySec/common.h" 
using namespace DegCod;
using namespace std;   // temp

/// \todo make more efficient and elegant
filesystem::path unique_filename(filesystem::path name)
{
    while(filesystem::exists(name))
    {
        std::string ext = name.extension().string();
        std::string nam = name.stem().string();
        //std::string pat = name.remove_filename() ;
        name = name.remove_filename() / filesystem::path( nam + "X" + ext) ;
    }
    return name;
}

bool    CMultSec::Export_from   ( CMultSec& base, bool only_selected)  
{
    filesystem::path dir, file;
    auto s= path();

    for ( const auto& CurMSec : base.MSecL())		// recorre todos las msec
    {    
        auto b= CurMSec->path();
        if ( b.empty() || s.find(b))  continue;    // finded OK only if s beging with p

        file = dir = s.replace(0, b.length()-1, CurMSec->_Path);
        file.remove_filename().replace_extension("fasta");
        dir.remove_filename().remove_filename();

        filesystem::create_directories(dir);
        Export_as(unique_filename(file).string(), only_selected);
        return true;
    }
    return false;
}

bool      CMultSec::Export_local_seq   ( CMultSec& base, bool only_selected)
{
    //assert(ms);
    filesystem::path dir, file;
    auto s= path();
    auto b= base.path();
    if ( s.find(b))  return false;            // finded OK only if s beging with b
    file = dir = s.replace(0, b.length(), base._Path);
    file.replace_extension("fasta");
    dir.remove_filename();

    filesystem::create_directories(dir);
    Export_as(unique_filename(file).string(), only_selected);
	return true;
}


CMultSec::CMultSec (	const std::string &path	, 
					std::shared_ptr<CSaltCorrNN>  NNpar	, 
					bool           all_dir  /*= false*/,
					float		   MaxTgId	/*= 100*/, 
					LonSecPosRang  SecLim	/*= LonSecPosRang {1,0}*/,	 
                    SecPosRang     SecLenLim/*= SecPosRang    {0,0}*/,
					bool           loadSec  /*= true*/
				 )  
    :	/*_name(trim_string(file)),	*/
	    _SecLim     (SecLim),
	    _MaxTgId    (MaxTgId), 
        _SecLenLim  (SecLenLim),
	    _NNPar      (NNpar)/*,
        _Path       (file)*/
{
	filesystem::path  itf(path);

	if (all_dir)                    // Load all files and directories recursiverly?
    {
        if (filesystem::is_regular_file(itf)) 
            itf.remove_filename();

	    _name = itf.filename ().string();     // The new MSec take the name of the dir.
        _Path = itf.string();                 // and the _Path point to it.

        filesystem::directory_iterator rdi{ itf }, end;

        for (; rdi != end; ++ rdi)
            AddMultiSec(  new CMultSec(  rdi->path().string() , 
                                         NNpar,
                                         filesystem::is_directory(rdi->status()), 
                                         MaxTgId, 
                                         SecLim,   
                                         SecLenLim,
                                         loadSec) );
        return;
	}
	else                           // Load only this file
	    if (itf.has_filename())
	    {
		    _name = itf.filename ().string();     /// The new MSec take the name of the file.
            _Path = itf.string();                 /// and the _Path point directly to the file.
            if (loadSec)
		       AddFromFile(path);  /// will throw if not a file
            return;
	    }
	throw std::ios_base::failure(string("No such sequence file: ")+ path );
    // throw "request a non recursive load from a non regular file";
}

int		CMultSec::AddFromFile (const std::string& file)	// return la cantidad de sec add ------  AddFromFile   ---
{	
    ifstream ifile( file ); 
	if ( ! ifile ) 
	{
	    throw std::ios_base::failure(string("Could not open the sequence file: ")+ file );
	}

        return AddFromFile(ifile); /// \todo: retrow anadiendo el nombre del file
}

int		CMultSec::AddFromFile (ifstream& ifile)		// return la cantidad de sec add -----------  AddFromFile   ---
{		
	//if (  _SecLim.Max() <= _SecLim.Min() ) 
	//	_SecLim.SetMax(0) ;                      // if ( _SecEnd<=_SecBeg) _SecEnd=0 ;

	int j=0;
	char c1;
	ifile>>skipws  >> c1;
	if ( ! ifile.good() ) 	
	{
	    throw std::ios_base::failure(string("Could not read the sequence file: ")/*+ file */);
	}

	if( c1 =='>' ) 																		
		return AddFromFileFASTA (ifile);
					// esto es FASTA, si no  BLAST o GB o ...
	if (c1 =='<' )
	{	
        string xml_DOCTYPE ;
		                    // <?xml version="1.0"?>
		                    // <!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
		if ( ! getline (ifile, xml_DOCTYPE,'>')  ||  ! getline (ifile, xml_DOCTYPE,'>') ) 	
		{	
            throw std::ios_base::failure(string("Could not open the sequence file: ")/*+ file*/ );		
        }

		if		(string::npos != xml_DOCTYPE.find("DOCTYPE BlastOutput") )				
			return AddFromFileBLAST (ifile);
		else if (string::npos != xml_DOCTYPE.find("DOCTYPE Bioseq-set" ) )				
			return AddFromFileGB (ifile);
		return 0;
	}
	if (c1 =='L' )	// LOCUS ?
	{	
        ifile.putback(c1);												
		return AddFromFileGBtxt (ifile);	
	}
	return 0;   // cerr unknow format		
}	 

int		CMultSec::AddFromFileFASTA (ifstream &ifile)  // -------------------    AddFromFileFASTA   ------------
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
	string Descriptor  ;
	while (getline (ifile, Descriptor) )
    {
		size_t name_end=Descriptor.find_first_not_of(
								"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567890!#$()*+-./<=>@[]^_{|}~"    );
		string Fasta_NAME = Descriptor.substr(0,name_end );
		Descriptor=Descriptor.substr(Fasta_NAME.length());

  		string Fasta_SEC ;					
		if (!getline (ifile, Fasta_SEC,'>')) 		break ;		

		if (     Fasta_SEC.length() < secBeg + lmin -1   )	
             continue;

		unique_ptr<CSec> sec (  new CSec(   Fasta_SEC , 
                                            _Local._NSec, 
                                            Fasta_NAME , 
                                            _NNPar,
                                            lmax, 
                                            secBeg  //  \todo: cambiar constr de CSec ???? use make_unique?
                                         )) ;                  assert(sec);
						
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
		        InsertSecAfter (idem, sec.release() ) ;	
            }
            else
            {
				(*idem)->Selected(false);
				(*idem)->Filtered(true);
		        InsertSec(idem, sec.release());
            }
		}
        else
            AddSec(sec.release() );

		NumSeq++;		
	}
	return NumSeq;	
}

int		CMultSec::AddFromFileBLAST (ifstream &fi) // ----------------  CMultSec::            AddFromFileBLAST  -----------------------------
{	
	///\todo adapt to multiquery BLAST

	unsigned int _BlastOutput_query_len ;		// x todos los "hits"
	int			 id=0;
	string       li ;  //  xml_line

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
  

	do  {	if ( !getline (fi, li,'>') )  return 0; }   // BLAST format error
	while  (string::npos==li.find("BlastOutput_query-len") ); fi>>_BlastOutput_query_len;//  <BlastOutput_query-len>267</BlastOutput_query-len>
	
	do 
    {	
        unsigned int	_Hit_num=0 ;			// para cada hit
		std::string	    _Hit_id  ;				
		std::string	    _Hit_def ;				// descriptor ??
		std::string	    _Hit_accession 	;
		long			_Hit_len=0 ;				
		float			_Hsp_bit_score=0 ;
		unsigned int	_Hsp_score=0 ;
		double			_Hsp_evalue=0 ;
		LonSecPos		_Hsp_query_from=0 ;    // dejar signed or unsigned !!!!????
		LonSecPos		_Hsp_query_to=0 ;
		LonSecPos 		_Hsp_hit_from=0 ;
		LonSecPos 		_Hsp_hit_to=0 ;
		int				_Hsp_query_frame=0 ;
		int				_Hsp_hit_frame=0 ;
		LonSecPos 		_Hsp_identity=0 ; // revisar type --- no sera %  : float??
		LonSecPos 		_Hsp_positive=0 ;
		LonSecPos 		_Hsp_gaps=0 ;
		LonSecPos 		_Hsp_align_len=0 ;
		std::string	    _Hsp_midline ;
		bool			_FormatOK=0 ;

		std::string	   sec;						// para CSec
		std::string	   nam;
		long		   l=0;
		std::string	   clas;

	    while(getline(fi,li,'>')&& string::npos==li.find("Hit_num"      ) ) ;  fi>>_Hit_num;                                //  <Hit_num>1</Hit_num>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hit_id"       ) ) ; if(!getline(fi, _Hit_id, '<') ) return id;    //<Hit_id>gi|84028434|gb|DQ318020.1|</Hit_id> 
	    while(getline(fi,li,'>')&& string::npos==li.find("Hit_def"      ) ) ; if(!getline(fi, _Hit_def,'<') ) return id;    //<Hit_def>Wets NIle virus strain ArB3573/82, complete genome</Hit_def>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hit_accession") ) ; if(!getline(fi, _Hit_accession,'<'))return id;//<Hit_def>Wets NIle virus strain ArB3573/82, complete genome</Hit_def>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hit_len"      ) ) ;  fi>>_Hit_len;				//  <Hit_len>11048</Hit_len> 
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_bit-score") ) ;  fi>>_Hsp_bit_score;			//  <Hsp_bit-score>482.786</Hsp_bit-score>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_score"    ) ) ;  fi>>_Hsp_score;              //  <Hsp_score>534</Hsp_score>		
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_evalue"   ) ) ;  fi>>_Hsp_evalue;		    	//  <Hsp_evalue>3.71782e-133</Hsp_evalue>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_query-from")) ;  fi>>_Hsp_query_from;		    //  <Hsp_query-from>1</Hsp_query-from>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_query-to" ) ) ;  fi>>_Hsp_query_to;			//  <Hsp_query-to>267</Hsp_query-to>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_hit-from" ) ) ;  fi>>_Hsp_hit_from;			//  <Hsp_hit-from>9043</Hsp_hit-from>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_hit-to"   ) ) ;  fi>>_Hsp_hit_to;			    //  <Hsp_hit-to>9309</Hsp_hit-to>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_query-frame"));  fi>>_Hsp_query_frame;		//  <Hsp_query-frame>1</Hsp_query-frame>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_identity" ) ) ;  fi>>_Hsp_identity;			//  <Hsp_identity>267</Hsp_identity>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_positive" ) ) ;  fi>>_Hsp_positive;			//  <Hsp_positive>267</Hsp_positive>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_gaps"     ) ) ;  fi>>_Hsp_gaps;	 		    //  <Hsp_gaps>0</Hsp_gaps>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_align-len") ) ;  fi>>_Hsp_align_len;		    //  <Hsp_align-len>267</Hsp_align-len>
	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_hseq"     ) ) ; if(!getline(fi, sec,'<'))         return id;//<Hsp_hseq>TACAACATGATGGGAAAGAGAGAGAAGAAG 
 	    while(getline(fi,li,'>')&& string::npos==li.find("Hsp_midline"  ) ) ; if(!getline(fi, _Hsp_midline,'<'))return id;//<Hsp_midline>|||||||||||||||||||||||||||||||||||||||||
                                                            

        LonSecPos  secHitBeg {secBeg - _Hsp_query_from +1};  // omitir estas primeras bases del Hit (de la parte del hit que tenemos)
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

        std::unique_ptr<CSecBLASTHit> secH
                { new CSecBLASTHit(      _BlastOutput_query_len ,
										// para cada hit
										_Hit_num ,
										std::move(_Hit_id) ,				
										std::move(_Hit_def) ,				
										std::move(_Hit_accession)	,
										_Hit_len ,				
										_Hsp_bit_score ,
										_Hsp_score ,
										_Hsp_evalue ,
										_Hsp_query_from ,
										_Hsp_query_to ,
										_Hsp_hit_from ,
										_Hsp_hit_to ,
										_Hsp_query_frame ,
										_Hsp_hit_frame ,
										_Hsp_identity ,
										_Hsp_positive ,
										_Hsp_gaps ,
										_Hsp_align_len ,
										std::move(_Hsp_midline) ,
										_FormatOK ,
										std::move(sec)		,
                                        lmaxHit,
                                        secHitBeg,          // _SecBeg, _SecEnd,
										id,			     //Hit_num   ???		//	char	*	nam,	Hit_def
										_NNPar           /*,  //long l=0,(Hit_len ---> NO ) !!!  -->_Hsp_align_len -OK clas,	conc*/
										)
                };

	    if ( secH->Len() < lmin  )	
            continue;

        if ( secH->_aln_fragment)   ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
        {
            if(secH->_aln_fragment->sq.Max()) 
                _Hsp_query_to    = _Hsp_query_from + secH->_aln_fragment->sq.Max()  -1;
            else
                _Hsp_query_to    = _Hsp_query_from + secH->Len() -1;
            _Hsp_query_from  = _Hsp_query_from + secH->_aln_fragment->sq.Min()-1;
        }
        else
        {
            secH->_aln_fragment.reset(new Aligned_fragment);      ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
            _Hsp_query_to    = _Hsp_query_from + secHitBeg + secH->Len() -1;
            _Hsp_query_from  = _Hsp_query_from + secHitBeg-1;
        }

        secH->_aln_fragment->sq_ref.Set(       _Hsp_query_from, _Hsp_query_to);    ///\todo review  ACTUALIZE !!!!!!!!!!!!!!!!!!!!
        secH->_aln_fragment->aln   .set(*this, _Hsp_query_from, _Hsp_query_to);
        secH->_aln_fragment->sq    .Set(       _Hsp_hit_from,   _Hsp_hit_to  );

        long MaxIdem= long(ceil((_Hsp_align_len*_MaxTgId)/100.0f));	// max of Id base 

        if( _Hsp_identity > MaxIdem )
        {
			secH->Selected(false);
			secH->Filtered(true);
            AddSec(secH.release() );
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
		            InsertSecAfter(idem, secH.release());
                }
                else
                {
					(*idem)->Selected(false);
					(*idem)->Filtered(true);
		            InsertSec(idem, secH.release());
                }
		    }
            else
               AddSec(secH.release() );
        }
		id++;		
    }
	while (fi.good() ); 
	return id; 
}

int		CMultSec::AddFromFileGBtxt (ifstream &ifile) // ----------------  CMultSec::            AddFromFileGBtxt  -----------------------------
{	const int gb_descr_w=12 ;	char gb_descr[gb_descr_w+1]; gb_descr[gb_descr_w]=0;
	size_t strl;
	int	id	=0     ;
	string txt_line ;

	do {	char	*	LOCUS			=nullptr      ;
			long		Seq_inst_length	=0		;	
			char	*	DEFINITION		=nullptr     ;
			char	*	ACCESSION		=nullptr     ;
			char	*	ORGANISM		=nullptr     ;
			// para CSec
			char	*	sec		=nullptr     ;		//	char		*	nam,	DEFINITION	,	//	long			l=0,		Seq_inst_length

			do  {	ifile>>setw (gb_descr_w)>>gb_descr ;	if ( ! ifile.good() ) return id; } 	// LOCUS       AY702040               10675 bp    RNA     linear   VRL 24-MAR-2005
			while  (strcmp(gb_descr,	"LOCUS"	) );
			ifile>>setw (gb_descr_w)>>gb_descr ;				LOCUS=new char[1+(strl=strlen(gb_descr))] ;
			strncpy(LOCUS,gb_descr,strl) ;						LOCUS	[strl]=0;
			ifile>>Seq_inst_length ;	

			do  {	ifile>>setw (gb_descr_w)>>gb_descr ;	if ( ! ifile.good() ) return id; } 	// DEFINITION  Dengue virus type 2 strain I348600, complete genome.
			while  (strcmp(gb_descr,	"DEFINITION"	) );
			getline (ifile, txt_line) ;								DEFINITION	=new char[txt_line.length()+1] ;
			txt_line.copy(  DEFINITION  ,txt_line.length()) ;		DEFINITION	[txt_line.length()]=0;

			do  {	ifile>>setw (gb_descr_w)>>gb_descr ;	if ( ! ifile.good() ) return id; } 	// ACCESSION   AY702040
			while  (strcmp(gb_descr,	"ACCESSION"	) );
			ifile>>setw (gb_descr_w)>>gb_descr ;			ACCESSION=new char[1+(strl=strlen(gb_descr))] ;
			strncpy(ACCESSION,gb_descr,strl) ;						ACCESSION	[strl]=0;

			do  {	ifile>>setw (gb_descr_w)>>gb_descr ;	if ( ! ifile.good() ) return id; } 	//   ORGANISM  Dengue virus 2
			while  (strcmp(gb_descr,	"ORGANISM"	) );
			getline (ifile, txt_line) ;							ORGANISM	=new char[txt_line.length()+1] ;
			txt_line.copy(  ORGANISM  ,txt_line.length()) ;		ORGANISM	[txt_line.length()]=0;

			do  {	ifile>>setw (gb_descr_w)>>gb_descr ;	if ( ! ifile.good() ) return id; } 	// ORIGIN      
			while  (strcmp(gb_descr,	"ORIGIN"	) );
			getline (ifile, txt_line,'/') ;							sec	=new char[txt_line.length()+1] ;
			txt_line.copy(  sec  ,txt_line.length()) ;				sec	[txt_line.length()]=0;

			CSecGBtxt *secGBtxt= new CSecGBtxt(	LOCUS       ,
												Seq_inst_length,	
												DEFINITION     ,
												ACCESSION      ,
												ORGANISM       ,
												sec	,	
												id,								//	char		*	nam,	DEFINITION	,	
												_NNPar);
				if ( secGBtxt->Len() >= _SecLenLim.Min()   )		
				{	
					CSec *idem=Idem(*secGBtxt);
					InsertSecAfter (secGBtxt  , idem) ;	
					if (idem) 
					{
						secGBtxt->Selected(false);
						secGBtxt->Filtered(true);
					}
					else
						id++;		
				}
				else delete secGBtxt;
			delete []sec ;
		}
	while (ifile.good() ); 
	return id; 
}

int		CMultSec::AddFromFileGB (ifstream &ifile)  // ----------------  CMultSec::            AddFromFileGB  -----------------------------
{	int		id=0;
	string xml_line ;

	do {	char		*	_Textseq_id_accession	=0 ;	
			char		*	_Org_ref_taxname		=0 ;
			char		*	_Seqdesc_title			=0 ;
			long			_Seq_inst_length		=0 ;	
			// para CSec
			char		*sec=0;			//char		*nam=0;		//long		 l=0;		//char		*clas=0;

			do  {	getline (ifile, xml_line,'>') ;	if ( ! ifile.good() ) return id; } 	// <Textseq-id_accession>DQ318020</Textseq-id_accession>
			while  (string::npos==xml_line.find(	"Textseq-id_accession"	) );
			getline (ifile, xml_line,'<') ;				_Textseq_id_accession=new char[xml_line.length()+1] ;
			xml_line.copy(_Textseq_id_accession,xml_line.length()) ;	_Textseq_id_accession	[xml_line.length()]=0;

			do  {	getline (ifile, xml_line,'>') ;	if ( ! ifile.good() ) return id; } // <Org-ref_taxname>West Nile virus</Org-ref_taxname>
			while  (string::npos==xml_line.find(	"Org-ref_taxname"	) );
			getline (ifile, xml_line,'<') ;								_Org_ref_taxname	=new char[xml_line.length()+1] ;
			xml_line.copy(  _Org_ref_taxname  ,xml_line.length()) ;		_Org_ref_taxname	[xml_line.length()]=0;
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
			do  {	getline (ifile, xml_line,'>') ;	if ( ! ifile.good() ) return id; }		// <Seqdesc_title>Wets NIle virus strain ArB3573/82, complete genome</Seqdesc_title>
			while  (string::npos==xml_line.find(	"Seqdesc_title"	) );
			getline (ifile, xml_line,'<') ;								_Seqdesc_title	=new char[xml_line.length()+1] ;
			xml_line.copy(  _Seqdesc_title  ,xml_line.length()) ;		_Seqdesc_title	[xml_line.length()]=0;
		// <MolInfo_biomol value="mRNA">3</MolInfo_biomol>
		// <MolInfo_completeness value="complete">1</MolInfo_completeness>
		// <Date-std_year>2006</Date-std_year>
		// <Date-std_month>1</Date-std_month>
		// <Date-std_day>1</Date-std_day>
		// <Seq-inst_mol value="rna"/>

			do  {	getline (ifile, xml_line,'>') ;	if ( ! ifile.good() ) return id; }  // GB format error
			while  (string::npos==xml_line.find("Seq-inst_length") ) ; ifile>>_Seq_inst_length;		// <Seq-inst_length>11048</Seq-inst_length>
			
		// <Seq-inst_strand value="ss"/>

			do  {	getline (ifile, xml_line,'>') ;	if ( ! ifile.good() ) return id; }	// <IUPACna>AGTAGTTCGCCTGTGTGAGCTGACA.... GGTGCTAGAACACAGGATCT</IUPACna>
			while  (string::npos==xml_line.find(  "IUPACna"  ) );
			getline (ifile, xml_line,'<') ;         sec=new char[xml_line.length()+1] ;
			xml_line.copy(sec,xml_line.length()) ;	sec[xml_line.length()]=0;	


			CSecGB *secGB=  new CSecGB(      _Textseq_id_accession,
											_Org_ref_taxname	,
											_Seqdesc_title,
											_Seq_inst_length	 ,
											sec	,	
											id,			//	char		*	nam,		Hit_def
											_NNPar/*,  	//	long			l=0,		(Hit_len ---> NO ) !!!  -->_Hsp_align_len -OK clas, conc*/
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
				else delete secGB;
		}
	while (ifile.good() ); 
	return id; 
}
	
int		CMultSec::AddFromFileODT (ifstream &ifileODT){return 0;}
int		CMultSec::AddFromFileODS (ifstream &ifileODS){return 0;}

CMultSec::LSec::const_iterator CMultSec::Idem ( CSec &sec )   // ------  CMultSec:: NotIdem  --- busqueda trivial de sec identicas -------------
{	
    if ( _MaxTgId >= 100  ) //  no restriction on similarity
        return _LSec.end() ;    

	long LenCandSec=sec.Len() ;     // Lenght of Candidate Sec (to be in the list, with MaxId)

	long MaxErCS= long(ceil(float(LenCandSec*(100.0f-_MaxTgId) ) / 100.0f)); // min of not Id base to be in the list
	
	for (auto CurSec = _LSec.begin(); CurSec != _LSec.end(); ++CurSec)		// recorre todos las primeras sec de esta misma ms
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

 CSec *	CMultSec::AddSec ( CSec *sec )
{	if (!sec) return nullptr ;
	_LSec.emplace_back(sec);
	UpdateTotalsAdding ( sec );
	return sec;
}
 CSec *	CMultSec::InsertSec(LSec::const_iterator pos, CSec *sec) 
{	if (!sec) return nullptr ;
	_LSec.emplace(pos, sec);
	UpdateTotalsAdding ( sec );
	return sec;
}
 CSec *	CMultSec::InsertSecAfter(LSec::const_iterator preSec, CSec *sec)
{	
	if (!sec) return nullptr ;
	_LSec.emplace(++preSec, sec);
	UpdateTotalsAdding ( sec );
	return sec;
}

void	CMultSec::UpdateTotalsAdding ( CSec *sec ) 
{	
	if (!sec || sec->_parentMS == this)										// no hay sec o ya estaba aqui
		return;
	CMultSec *parMS   =sec->_parentMS;					// /*._Get()*/
	CMultSec *My_parMS=     _parentMS;					// /*._Get()*/
	bool checkExtr(true) ; 
	CMultSec   *cp;

	if (parMS)
	{
	    cp=findComParent( parMS);
		parMS->_Local._NSec--	;								// la elimino de la ms orig. 
		//parMS->_Global._NSec -- ;								// 
		if (parMS->_Local._NSec)
			if (parMS->isLocExtreme(sec))
				RecalExtremes();
			else
				checkExtr=false;
		else
			checkExtr=false;

		for ( /*parMS=parMS->*/_parentMS;  parMS != cp &&  parMS ;  parMS=parMS->_parentMS)		// desde localizacion orig subiendo hasta parent comun
			{
				parMS->_Global._NSec -- ;								// elimino  de este total.
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
	Add2LocalExtreme(*sec);
	for (My_parMS ; My_parMS!=cp && My_parMS; My_parMS=My_parMS->_parentMS)  // desde mi hacia arriba hasta el com parent anadiendo
	{
		if (checkExtr)
			My_parMS->Add2GlobalExtreme(*sec);
		else
		{
			parMS->_Global._NSec ++;								// sumo sus s a este total.
		}
	}	
	sec->_parentMS = (this) ;							//* std::weak_ptr<CMultSec> */
}


CMultSec   *CMultSec::findComParent( CMultSec *ms)
{
	if(!ms || ms==this) 
		return ms;
	std::stack<CMultSec*> myTree,		oTree;
	CMultSec			*myPms=this, *oPms=ms;
		myTree.push(myPms);
		oTree.push(oPms);

	do 
	{	myPms=myPms->_parentMS;
		myTree.push(myPms);
	}while (myPms);

	do 
	{	oPms=oPms->_parentMS;
		oTree.push(oPms);
	}while (oPms);

	do 
	{	if(  (oPms = myTree.top())   !=  oTree.top()  )
			return myPms;
		myPms= oPms;
		myTree.pop();
		oTree.pop();
	} while (!myTree.empty() && !oTree.empty() );
	return myPms;
}


CMultSec *	CMultSec::AddMultiSec ( CMultSec *ms )  //--------------------------------------    AddMultiSec    --------------------
{	if (!ms) return nullptr;	
	_LMSec.emplace_back(ms);
	UpdateTotalsAdding ( ms );   // al llamar ya esta la ms movida fisicamente. Falta solo actualizar extremes
	return ms;
}
void	    CMultSec::UpdateTotalsAdding ( CMultSec *msec ) 
{	
	if (!msec || msec->_parentMS==this)					// no hay msec o ya estaba antes en una de mis subtrees inmediatas. 
		return;

	CMultSec *parMS   =msec->_parentMS;					// /*._Get()*/
	CMultSec *My_parMS=     _parentMS;					// /*._Get()*/
	bool checkExtr(true) ; 
	CMultSec   *cp;
	if (parMS)										// no es imprescindible. Anadido solo por claridad de intencion
	{	cp=findComParent( msec);
		for ( parMS;  parMS != cp   ;  parMS=parMS->_parentMS)			// desde localizacion orig subiendo hasta parent comun
			{
				parMS->_Global._NSec -= msec->_Global._NSec ;			// elimino sus s de este total.
				parMS->_Global._NMSec-= msec->_Global._NMSec + 1;		// elimino sus ms de este total.
				parMS->_Local._NMSec--	;								// la elimino de esta ms. 
				if (checkExtr && parMS->_Global._NSec)
					if (parMS->isGlobExtreme(msec))
						RecalExtremes();
					else
						checkExtr=false;
			}
	}else
		cp=nullptr;

	Add2LocalExtreme(*msec);
	for (My_parMS ; My_parMS!=cp && My_parMS; My_parMS=My_parMS->_parentMS)  // desde mi hacia arriba hasta el com parent anadiendo
	{
		if (checkExtr)
			My_parMS->Add2GlobalExtreme(*msec);
		else
		{
			My_parMS->_Global._NSec += msec->_Global._NSec ;			// sumo sus s a este total.
			My_parMS->_Global._NMSec+= msec->_Global._NMSec + 1;		// sumo sus ms a este total.
		}
	}
	msec->_parentMS = (this) ;							// std::weak_ptr<CMultSec> 

}

		CMultSec::~CMultSec ()				// funciona bien solo si la lista es "lineal"
{	
	//_LSec.Destroy  ()	;
	//_LMSec.Destroy ()	;
	//Remove(); 
	// CSec		*_Consenso ;
}    


//void	CMultSec::RefreshExtremes( CMultSec *ms)
//{
//	if(!ms)
//		return;
//	if(ms->_NSec)
//
//}
	//	
	//while (parMS)									// subo por el tree hasta llegar al root
	//	if (parMS==this)						
	//		return;									// la sec estaba antes en una de mis subtrees. 
	//	else
	//auto TLen= msec->_TLen ;
	//auto TTm = msec->_TTm ;
	////auto TNMS
	//CMultSec *My_parMS=this ;    //_parentMS/*._Get()*/;		// 
	//while (My_parMS)								// subo por el tree hasta llegar al root
	//{	if ( ! My_parMS->_TNSec  )	
	//	{											
	//		My_parMS->_TLen= TLen ;				//   si ademas es la primera del todo inicializar los max totales
	//		My_parMS->_TTm = TTm ;
	//	}else
	//		{	
	//			My_parMS->_TLen.Expand( TLen) ;
	//			My_parMS->_TTm .Expand( TTm ) ;
	//		}	
	//	My_parMS->_TNSec +=msec->_TNSec;			// la elimino de este total. Que hacer con las max??
	//	My_parMS->_TNMSec+=msec->_TNMSec;
	//	TLen=My_parMS->_TLen ;
	//	TTm =My_parMS->_TTm ;
	//}
