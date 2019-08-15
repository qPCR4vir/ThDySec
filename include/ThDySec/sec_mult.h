/**
* Copyright (C) 2009-2019, Ariel Vina Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2009-2019
*
* @file  ThDySec\include\ThDySec\sec_mult.h
*
* @brief 
*/

#pragma unmanaged	

#ifndef _MULTSEC_H
#define _MULTSEC_H

#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <string>
#include <memory>
#include <vector>
#if defined(BOOST_FILESYSTEM_FORCE) || defined(NANA_FILESYSTEM_FORCE)
#  include <nana/filesystem/filesystem_ext.hpp>
#else
#  include <filesystem>
#endif
#include <list>


#include "sec.h" 
#include "sec_rang.h" 

                                        // ---------------------------------------------------------- 	CMultSec    -------------------
  /// Permite hacer grupos de sec o de MultiSec (para analisis por "especies"?)

  /// \todo add construction of concensus 
  /// \todo make std with vector or list of std::shared_ptr<Csec> and CMultSec> ?
  /// Atention:  It owns and destroy the sequences: Use Remove() or Free() to prevent destruction fo sequences
class CMultSec	 
{	public:
	using LSec  = std::list<std::shared_ptr<CSec    >>;
	using LMSec = std::list<std::shared_ptr<CMultSec>>;

		std::string			_name ;						///< 
        int					_ID       {NewMS_ID()};		///< Unique ID in each programm run
        LonSecPosRang       _SecLim   {1,0};			///< \todo: quitar de aqui?. Pertenece a CSec, o a un objeto "AddFromFile" 
        SecPosRang          _SecLenLim{0,0};            ///< \todo: NumRang<LonSecPos> _SecRange{1,0};// def-unlimited: the range used, not the limits
		float				_MaxTgId  {100};			///< \todo: quitar de aqui?. Pertenece a CSec, o a un objeto "AddFromFile" 
		std::shared_ptr<CSaltCorrNN>	_NNPar ;		///< \todo: quitar de aqui?. Pertenece a CSec, o a un objeto "AddFromFile" 
		CMultSec			*_parentMS {nullptr};		///< std::weak_ptr<CMultSec> _parentMS	;
        CSec				*_Consenso {nullptr};
        bool                 _selected { true };
		std::string			 _orig_file_path ;				///< file path of the original sequence source


     //explicit CMultSec (const std::string &Name  )                 : _name		(trim_string(Name))  {	}

       /// Create a named and free empty group
       explicit CMultSec (std::shared_ptr<CSaltCorrNN> NNpar, const std::string &Name = "")        
                          : _NNPar      (NNpar            ), 
   	  			            _name       (trim_string(Name))            
                      {  }

	   /// Create a named and free empty group using the given group as "template"
       CMultSec(CMultSec	*ms, const std::string &Name = ""): _name       (trim_string(Name)),
                                                                _SecLim     (ms->_SecLim),
                                                                _SecLenLim  (ms->_SecLenLim),
                                                                _MaxTgId    (ms->_MaxTgId), 
                                                                _NNPar      (ms->_NNPar)         
                  {  }

	   /// Read all the sequences from a stream into this group dedducing format and appling filters
       CMultSec (	std::ifstream &	    file	,	 
					std::shared_ptr<CSaltCorrNN>  NNpar	, 
					float		  MaxTgId	= 100,                 ///< Sec. with more % of identity are marked as "filtered" and not selected
					LonSecPosRang  SecLim	= LonSecPosRang {1,0}, ///< Filtre, using only this region. Will take into account alignment coordenates.
                    SecPosRang     SecLenLim= SecPosRang{0,0})     ///< Limit the length. tiny sec: not created, large: get trunkated   
                  : 
	                    _SecLim     (SecLim),
                        _SecLenLim  (SecLenLim),
	                    _MaxTgId    (MaxTgId), 
	                    _NNPar      (NNpar)              
                  { AddFromFile(file); }

         /// The new MSec take the name of the dir, and remember the rest of the path
         CMultSec (	const std::string &path	,                       ///< The name of the file or directory to be loaded 
					std::shared_ptr<CSaltCorrNN>  NNpar	, 
					bool           all_dir  = false,                ///< Load all files and directories recursiverly? 
					float		   MaxTgId	= 100,                  ///< Sec. with more % of identity are marked as "filtered" and not selected
					LonSecPosRang  SecLim	= LonSecPosRang {1,0},  ///< Filtre, using only this region. Will take into account alignment coordenates.
                    SecPosRang     SecLenLim= SecPosRang    {0,0},  ///< Limit the length. tiny sec: not created, large: get trunkated
					bool           loadSec  = true                  ///< Get the sec? False: get only the dir/file structure
				 ) ;


	    bool	Selected(bool select)	{return _selected=select;} 			///< \todo make protected: ??
	    bool	Selected(		) const 
		{return _selected && (_parentMS ? _parentMS->Selected() : true);  }				///< User-editable


        /// Construct a full-current path acording to the current tree
		
		/// can be different from the original path saved in member variable ._Path
		static std::string	Path(CMultSec *ms, const std::string& path_sep="/")
        {
			std::string path  ;			 
			for (CMultSec *parent=ms; parent;parent=parent->_parentMS)
				path = parent->_name + path_sep + path;
			return path;
        }

        /// Construct a filesystem path acording to the current tree, which can be different from the original path saved in member variable ._orig_file_path
		std::string	path( )
		{
			std::string sep(std::string(1, std::filesystem::path::preferred_separator));// ::slash<std::filesystem::path>().value));
            return Path(this, sep);
		}

		/// generate a runing unique ID
		static int			NewMS_ID()
		{
			static int ID(0);
			return ++ID;
		}

		struct CExtremes
		{
			int			_NSec{}, _NMSec{};
		    NumRang<LonSecPos> _Len		;		
		    NumRang<Temperature> _Tm	;
		    // CExtremes():_NSec(0), _NMSec(0)	{    }
		    void Set   (const CSec& s)	{	 
										    _Len.Set(s.Len());		_Tm.Set(s._Tm)	;
									    }
		    void Expand(const CSec& s)	{	
										    _Len.Expand(s.Len());	_Tm.Expand(s._Tm);
									    }
		    void Set   (const CExtremes& e)	{	 
											    _Len.Set(e._Len);		_Tm.Set(e._Tm)	;
										    }
		    bool Expand(const CExtremes& e)	{	
											    bool L =_Len.Expand(e._Len);  
												bool T = _Tm.Expand(e._Tm);
												return ( T || L );
										    }
		    void Clear(){_NSec=0, _NMSec=0;}

		    bool isExtreme(const CSec&		s){return _Tm.isExtrem( s._Tm ) || _Len.isExtrem( s.Len() );}
		    bool isExtreme(const CExtremes& e){return _Tm.isExtrem( e._Tm ) || _Len.isExtrem( e._Len  );}

		} _Local, _Global;

		void Add2LocalExtreme(const CSec& s)
		{ 
			if (_Global._NSec)
			{
				_Global.Expand(s);
				_Local._NSec ? _Local.Expand(s) : _Local.Set(s);
			}else
			{
				_Local.Set(s);
				_Global.Set(s);
			}
			_Local._NSec++;
			_Global._NSec++;
		}
		void Add2GlobalExtreme(const CSec& s)
		{ 
			if (_Global._NSec)
			{
				_Global.Expand(s);
			}else
			{
				_Global.Set(s);
			}
			_Global._NSec++;
		}
		bool Add2LocalExtreme(const CMultSec& ms)
		{
			_Local._NMSec ++ ;
			return Add2GlobalExtreme(ms);
		}
		bool Add2GlobalExtreme(const CMultSec& ms)
		{	
			_Global._NMSec += 1 + ms._Global._NMSec   ;
			if (! ms._Global._NSec )
				return false;
			bool res=true;
			if (_Global._NSec  )
				res=_Global.Expand(ms._Global); 
			else
				_Global.Set(ms._Global) ;
			_Global._NSec  += ms._Global._NSec   ;
			 return res;
		}
		void RecalExtremes()
		{
			_Local.Clear();
			for (  auto &CurSec :  _LSec )		// recorre todos las primeras sec
				Add2LocalExtreme( *CurSec) ; 
		}
		void RecalGlobExtremes()
		{
			_Global=_Local;
            for (auto &CurMSec : _LMSec)  		// recorre todos las primeras sec
				Add2LocalExtreme( *CurMSec) ; 
		}
		bool isGlobExtreme(const CSec     *s ){return _Global.isExtreme( *s				) ;}
		bool isGlobExtreme(const CMultSec *ms){return _Global.isExtreme( ms->_Global	) ;}
		bool isLocExtreme (const CSec     *s ){return  _Local.isExtreme( *s				) ;}
		//bool isLocExtreme (const CMultSec *ms){return _Tm.isExtrem (ms->_TTm) || _Len.isExtrem (ms->_TLen );}
		//void setGloExtreme(const CMultSec *ms){return _Tm.isExtrem (ms->_TTm) || _Len.isExtrem (ms->_TLen );}

		//int			AddFromDir		(const std::string& dir , bool  recurs  /*= false*/)

		/// determine the file format and read sequences into this group

		/// The file format is decided lookind at the first non blanc character: > fasta, < some xml, etc.
		/// @return number of readed sequences 
		/// \todo better exeption mesg and new exeptio classes for format error?
		int		AddFromFile		(const std::string& file);
		int		AddFromFile     (std::ifstream& ifile);
		int		AddFromFileFASTA(std::ifstream &ifileFASTA);
		int		AddFromFileBLAST(std::ifstream &ifileBLAST);
		int		AddFromFileGB	(std::ifstream &ifileGB);
		int		AddFromFileGBtxt(std::ifstream &ifileGB);
		int		AddFromFileODT	(std::ifstream &ifileODT);
		int		AddFromFileODS	(std::ifstream &ifileODS);

		/// Recursively export sequences to new files in fasta format with filters applied

        /// Reproduce the current -in memory- tree, creating directories as need, 
        /// and export the local sequences in files with extention .fasta.
        /// If the file allready exist create a file with a new name
        /// From each child MSec of this MSec Export one file in fasta format  
        /// (will try all child of base to set the base dir and export)
        /// To decide the name of the file it need the path of the base node that do not form part of the file name
        /// That is for example: "all_seq/Primers for Multiplex PCR/"...
        bool    Export_from   ( CMultSec& tree_base, bool only_selected)  ;

		/// Export local sequences to a new file in fasta format with filters applied

        /// Export the local sequences in a file with extention .fasta.
        /// The name is generated acording to the current postion of the group on the tree.
        /// If the file allready exist create a file with a new name
        /// To decide the name of the file it need the path of the base node that do not form part of the file name
        /// That is for example: "all_seq/Primers for Multiplex PCR/"...
        bool    Export_local_seq   ( CMultSec& base, bool only_selected);

        void   Export_as(std::string filename, bool only_selected)  
        {
			std::ofstream ofile( filename );
	        if ( ! ofile ) 
	        {
	            throw std::ios_base::failure(std::string("Could not create the sequence file: ")+ filename );
	        }
            Export( ofile, only_selected);
        }
        void   Export(std::ofstream& ofile, bool only_selected)
        {
        	for ( auto& CurSec : _LSec   )		// recorre todos las sec locales
				if (CurSec->Selected() || !only_selected)
                    CurSec->ExportFASTA(ofile) ; 

            for (auto& CurMSec : _LMSec  )		// recorre todos las msec
			    CurMSec->Export(ofile, only_selected);
        }


		LSec::const_iterator Idem (CSec &sec);  //		CConsParam	_ConsPar ;
		CSec		*AddSec			( CSec *sec );
		CSec		*InsertSec		( LSec::const_iterator pos, CSec *sec ) ;
		CSec		*InsertSecAfter ( LSec::const_iterator pos, CSec *sec ) ;
		CMultSec    *MoveMSec       (LMSec::const_iterator from) // revise design
		{
			CMultSec* p = (*from)->_parentMS;
			if (p == this) return p;    // no-op ; mover al final?
			LMSec::value_type s = std::move(*from);
			if (p) p->_LMSec.erase(from);
			return AddMultiSec(s);
		}
		int			CountSelectedSeq		()
		{
            int count{0};
			for (auto& CurSec : _LSec)		// recorre todos las primeras sec de esta misma ms
				 count += (CurSec->Selected()) ;
			return count;
		}
		int			CountSelectedSeqRec		()
		{
            int count{0};
				 count += CountSelectedSeq() ;
			for (auto& CurMSec : _LMSec)			// recorre todos las primeras sec
				 count += (CurMSec->CountSelectedSeqRec());
			return count;
		}
		int			CountSelectedNDegSeq	(int MaxGrDeg=0)
		{
			int count(0);
			for (auto& CurSec : _LSec)		// recorre todos las primeras sec de esta misma ms
			{
				if( CurSec->Filtered() ) 
				{	if(CurSec->NonDegSet())
					{	
						if (!MaxGrDeg)
							count += CurSec->NonDegSet()->_Global._NSec ;
						else	
							if( MaxGrDeg > CurSec->NonDegSet()->_Global._NSec)
								count += CurSec->NonDegSet()->_Global._NSec ;
					}
					else ++count;
				}
			}
			return count;
		}
		int			CountSelectedNDegSeqRec	(int MaxGrDeg=0)
		{
			int count(0);
				 count += CountSelectedNDegSeq( MaxGrDeg) ;
			for (auto& CurMSec : _LMSec)		// recorre todos las primeras sec
				 count += CurMSec->CountSelectedNDegSeqRec( MaxGrDeg);
			return count;
		}

		void			CreateNonDegSet	( )
		{
			for (auto& CurSec : _LSec)		// recorre todos las primeras sec de esta misma ms
			  CurSec->CreateNonDegSet();
		}
		void			CreateNonDegSetRec	( )
		{
			CreateNonDegSet( ) ;
			for (auto& CurMSec : _LMSec)		// recorre todos las primeras sec
				 CurMSec->CreateNonDegSetRec( );
		}


		CMultSec	*AddMultiSec	(const std::string &Name )
		{
			return AddMultiSec(new CMultSec (this, Name  ));
		}
		CMultSec	*AddMultiSec	(CMultSec *MultSec);
		CMultSec	*AddMultiSec	(LMSec::value_type MultSec);

//		CSec		CalculateConsenso	(double) ;
//		void		release			()	{_LSec.release(); _LMSec.release();}
		void		clear			()	{_LSec.clear(); _LMSec.clear();} ///\todo uncount !!!!!

		virtual ~CMultSec ()  ;	

        long Len()
           {
               if (_Consenso)
                   return _Consenso->Len();
               else return 0;
           }

		CMultSec		*findComParent		( CMultSec	*ms	);
		const LSec &  SecL() const { return  _LSec; }
		const LMSec& MSecL() const { return _LMSec; }
	private:
		LSec  _LSec;
        LMSec _LMSec;
		//std::list<CSec    > _LSec;
		//std::list<CMultSec> _LMSec;   
        //std::list<std::shared_ptr<CSec    >> _LSec;
        //std::list<std::shared_ptr<CMultSec>> _LMSec;
		//void			UpdateTotals		( CSec		*sec ) ;
		void			UpdateTotalsMoving	( CSec		*sec ) ;
		void			UpdateTotalsMoving	( CMultSec	*sec ) ;
		//static void		RefreshExtremes		( CMultSec	*ms	);
};

#endif


