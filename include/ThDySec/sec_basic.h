/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDySec\sec_basic.h
*
* @brief 
*/

#ifndef _SEC_BASIC_H
#define _SEC_BASIC_H

#pragma unmanaged	

#include <fstream>
#include <string>
#include <memory>

#include "cod_deg.h"
#include "sec_rang.h"
#include "common.h" 

class CMultSec	;

class ISec	: public CLink		// Pure virtual class ?
{public:			
	using sequence = std::basic_string<DegCod::Base> ;

	const char*          charSequence()	const 
       {return reinterpret_cast<const char*>(Sequence().c_str()  );} 

	virtual const sequence& Sequence(							  )	const=0 ;

	// std::unique_ptr<ISec> ?
	virtual ISec *Clone    ( DNAstrand   strnd = DNAstrand::direct )	const=0 ; ///< unique_ptr<ISec> crea una copia muy simple. CUIDADO con copias de CSecBLASTHit y otros derivados
	virtual ISec *Clone    ( long        InicBase,
							 long        EndBase, 
							 DNAstrand   strnd = DNAstrand::direct)	    const=0  ;

	virtual std::string Copy_Seq( DNAstrand   strnd = DNAstrand::direct)const=0;
	virtual std::string Copy_Seq( long InicBase,
		                          long  EndBase,
		                          DNAstrand   strnd = DNAstrand::direct)const=0;
	virtual				~ISec			() { }
};

/// Have a letter sequence but no code sequence and no thermodynamics - abstract because of clone ?!-fix?
class CSecBasInfo : public ISec
{ protected:	
	bool			_selected{true}, _filtered{false};
	int				_ID{NewS_ID()} ;			          ///< num de la sec en file original?? en total??, num unico?
	std::string     _name ;			                      ///< nombre unico? FASTA id	// char			*_name ;
	std::string		_description;
	//long			_len ;			                      ///< longitud corregida, sin los '$'
	long			_GrDeg {1} ;		                  ///< Grado de degeneracion. Cantidad total de diferentes molec, dependiendo de deg
	float			_GCp   {} ;		
	long			_Count[DegCod::n_dgba];	              ///< An array of counters for each deg base - inicializar!!
	long			_NDB   {};			                  ///< cantidad de bases deg
	std::string     _Clas ;		                          ///< clasificacion
    sequence	    _c=sequence{ DegCod::basek[DegCod::n_basek-1]};		  ///< sec char, comienzan y terminan con '$'0
	CMultSec		*_NonDegSet{nullptr} ;                /// \todo: std::unique_ptr<>
	static int	NewS_ID     ()	{static int last_ID{};	return ++last_ID;	}

		CSecBasInfo (int id,     const std::string& nam    , const std::string& clas) 						
			           :	_ID	 ( id ), 								
							_name( trim_string	(nam )),	
						    _Clas(clas )/*,
							_c{2,DegCod::basek[DegCod::n_basek]}*/
                       {} 

		CSecBasInfo()  {}
		CSecBasInfo ( long l)    {_c.reserve(l+2);}   ///   ????????????????????
public:
    std::unique_ptr<Aligned_fragment> _aln_fragment;
    void ExportFASTA(std::ofstream& ofile, int line_len=80)
    {
        ofile << "\n"
              << ">" << _name << " " << Description ()   ;
        for (int i=0 ; i< Len() ; ++i )
        {
            if (!(i % line_len))    // after description, i=0 print end line too
                ofile << "\n";

            ofile << charSequence()[i+1]; // we dont export the first letter: $, and the last
        }
        ofile << "\n";
     }

	//CSecBasInfo			*CopyFirstBases	(long pos)	;			 // copia parcialmente hasta la pos
    ~CSecBasInfo() override;
    std::string Name		()const		{return _name;}              ///< User-editable
	int			ID			()const		{return _ID;}	             ///< Run-time-sistem define, non-editable
	bool		Filtered(bool filter)	{return _filtered=filter;}   ///< User-editable ?????
	bool		Filtered(		) const {return _filtered;}
	bool		Selected(bool select)	{return _selected=select;} 		  	 ///< make protected: ??
	bool		Selected(		) const {return _selected;}					 ///< User-editable
	void	    Description (std::string	description)		{ _description.swap(description);}
	virtual std::string	Description ()const	
        {
            return _description.empty() ? Name()  : _description ; 
        }

	const sequence& Sequence	(						)	const override   /// ver esto !!!
		{	 
			return _c ;
		}

	std::string  Copy_Seq(long InicBase,
		                  long  EndBase,
		                  DNAstrand   strnd = DNAstrand::direct)	    const override;
	std::string  Copy_Seq(DNAstrand   strnd = DNAstrand::direct)	    const override
	{
		return Copy_Seq(1, Len(), strnd);
	}

	bool		 NotIdem		(CSecBasInfo *sec) {return false;}
	DegCod::Base operator[]	(int i)const{return _c[i];}  /// i+1 ????

	long		Len			()const		{return static_cast<long>(_c.length())-2;} //
	long		Degeneracy	()const		{return _GrDeg;}
	long		*BaseCount	()			{return _Count;}
	long		BaseCount	(DegCod::Base b) const{ return  DegCod::is_degbase[b] ?  _Count[DegCod::db2nu[b]] : 0;}
	CMultSec	*NonDegSet	()			{return _NonDegSet;}
	float		GCpercent	()const		{return	_GCp ;}		
};


#endif



