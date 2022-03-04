/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDySec\th_dy_align.h
*
* @brief  Thermodynamic Alignment Algorithm
*
* This representation is based on the ideas and code of Kaderali (http://bioinformatics.oxfordjournals.org/content/21/10/2375.abstract)
* but with many modifications, so that the original authors have no responsibility on the many errors,
* simplifications or inconsistencies I have introduced (most file and class names were changed to avoid confusion with originals).
*
* Some additions/modifications are:
*
* - Integrate different variants of the algorithm into a (virtual) class hierarchy
* - use the code (Base) from CSec instead of char* as input sequences
* - open the dynamic program (DP) black-box to output all hits passing a cutoff (it originally report only the best hit)
* - class CHit to managed hits
* - split each of the two DP matrix into three, because the 0,1,2 index was not used as variable, but as a name
* - introduce a step matrix to go back in the DP matrix with minimal calculation
* - use a shared_ptr<CSaltCorrNN> member
* - idiomatic typedefs: Energy, Entropy, Temperature, etc. 
* - the initial maxlen of probes and targets are only a hint, and the DPT is resized as need. \todo implement as std::vector
* - ThDyAlign_TmCand -designed specially to find probes that hybridize in both sequences and in others targets too.
* -
* -
*
* The original source files have the following headers:
*
*     //=============================================================================
*     // Module:        thermalign.h
*     // Project:       Diploma Thesis - Probe Selection for DNA Microarrays
*     // Type:          header file - Thermodynamic Alignment.
*     // Language:      c++
*     // Compiler:      microsoft visual c++ 6.0, unix/linux gcc
*     // System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
*     // Database:      none
*     // Description:   class CThermAlign - Thermodynamic Alignment Algorithm
*     // Author:        kaderali
*     // Date:          9/2000 - 12/2000
*     // Copyright:     (c) L. Kaderali, 9/2000 - 12/2000
*     //
*     // Revision History
*     // $ 00sep04 LK : created
*     // $ 00dec30 LK : changed to do local alignment of probe against
*     //                one entire sequence
*     // $ 01feb07 LK : optimized
*     // #$
*     //=============================================================================
*
*     //=============================================================================
*     // Module:        galign.h
*     // Project:       Cubic Project - Calculation of melting temperature and free
*     //                energy of two DNA strands
*     // Type:          header file - Thermodynamic Alignment.
*     // Language:      c++
*     // Compiler:      microsoft visual c++ 6.0, unix/linux gcc
*     // System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
*     // Database:      none
*     // Description:   class GAlign - Thermodynamic Alignment Algorithm
*     // Author:        leber
*     // Date:          01/2002 - 02/2002
*     // Copyright:     (c) L. Kaderali & M. Leber, 01/2002 - 02/2002
*     //
*     // Revision History
*     // $ 00sep04 LK : created
*     // $ 00dec30 LK : changed to do local alignment of probe against
*     //                one entire sequence
*     // $ 01feb07 LK : optimized
*     // #$
*     //=============================================================================
*
* Which are accessible under GNU GPL at: http://dna.engr.uconn.edu/?page_id=85
*
*/

#pragma unmanaged	
#ifndef _TH_DY_ALIGN_H
#define _TH_DY_ALIGN_H
//#define _CRT_SECURE_NO_WARNINGS
#include <assert.h>
#include "sec_mult.h"
#include "common.h" 


class CHit;
class CHitAligned ;
//extern char sep[];


//  ------------------------------------------------------   ThDyAlign  ----------------------------------

/// base for thermodynamic dynamic program alignment calculation classes

/// "hybrid" the first sequence (the sonde) as simple-strand onto the 
///  complementary (by default) strand of the double-stranded second sequence (the target). 
/// It uses the code (Base) from CSec as input sequence 
/// Construct once, and reuse for each new paar sonde/target, and after each use of Align()
/// use the "OUTPUT" functions to return results
/// The "OUTPUT" functions Get_X(i,j)- return parameter x at (i,j) (dynamic-table-coordinates-position in sequences) for maximized algorithm parametr
/// while the other "OUTPUT" x() functions return parametr x for maximized algorithm parametr
/// \todo: make abstract class ?   do what ThDyAlign_Tm do... simplify?
/// \todo: make hybridization on both strand simultaneously
class ThDyAlign												
{
 public:		
	friend class	CHitAligned ;	 
	friend class	CHit ;

	ThDyAlign	(LonSecPos MaxLenSond, 
				LonSecPos MaxLenTarg, 
				std::shared_ptr<CSaltCorrNN>  NNpar, 
				float InitMax= 0 );
	virtual	~ThDyAlign	() ;
	
	                               /// just return the name of the method
	virtual std::string AlignMeth() 
	{  return "ThAl"; }

	                              /// will be used in the next "experiment"
	Temperature		SetTa	(Temperature Ta)					
	{
		Temperature T= _Ta;		
		_NNpar->_Ta  = _Ta=Ta;	   
		return T;
	}

	  /// return settled Ta, which will be used in the next "experiment"   \todo: review all around Ta 	
	Temperature	     Ta		()const							
	{  return _Ta; } 

	void			SetSig	(Temperature Tm_sig, Energy G_sig)	///< Tm and G cut off values
	{ _Tm_sig=Tm_sig;	_G_sig=G_sig ;}	

	       /// main trigger of alignment calculation of the sonde on the complementary strand of the target \todo: use const references
	virtual void	Align	( CSec *sonde,  CSec *target)		
	{	ClearHits();
		Use	(sonde, target);	//	InitBorder	();
		Run();
		SelectOptParam();	// a la Ta en que se calculo
	}

	        /// Ta used for the last calculation (-1 before any calculation) 
	Temperature		last_Ta()const{return _Ta_lastCalc;}  

	void			force_ResizeTable();		    //	void force_ResizeTable(long LenSond, long LenTarg) ;
	void			ResizeTable(LonSecPos LenSond, LonSecPos LenTarg);

	///  back pointer matrix  a:i+1 -direct sonde,    b:j+1 -dir tg,   d:diag
	enum Step {	st_a=1, st_i=1, st_b=2, st_j=2, st_d=3, NoStep=0,		 
				st_0=0, st_1, st_2, st_3, st_4, st_5, st_6, st_7,
				c_01=01, c_02=2,c_04=4, c_05=5, c_06=6, c_08=8, c_09=9, c_10=10, 
			  } *_pre	;   /*,*_pre0, *_pre1, *_pre2*/ 
protected:
	void			Use			(CSec  *sonde, CSec *target);            ///< \todo: take const ?
	void			Run			(LonSecPos tg_j_StartPos =1);            ///< the actual calculation
	void			ClearHits	()							{	_Hits.clear(); }
	virtual	bool	AddIfHit	(LonSecPos i, LonSecPos j)	{	return false;    }
	void			ResizeTable	()							{	ResizeTable(_sd->Len() , _tg->Len() ) ;}
private:
	void			InitBorder	();
public:
	float	(ThDyAlign::*CalcParam)(Entropy S, Energy H)const ; ///<  \todo: make it virtual ????
			Temperature	CalcParamTm(Entropy S, Energy H) const {return    _NNpar->CalcTM (S, H);} 
			Energy		CalcParamG (Entropy S, Energy H) const {return  - _NNpar->CalcG  (S, H);} 
			Temperature	CalcParamRs(Entropy S, Energy H) const {	Temperature Tm = _NNpar->CalcTM (S +_restS, H +_restH);	//  EXPERIMENTAL
																	if ( _NNpar->CalcTM (S , H )==0 ) Tm= 0;		//if ( Tm < _minTm ) return 0;
																	Tm -= _sd->_Tm.Ave();
																	return 		Tm ;
																}  


			// "OUTPUT" functions to return results. Use only after Align()  !!!

protected:						
	         // --------- Get_X(i,j)	---------- return parameter at (i,j) for maximized algorithm parameter !!!!    ?????????

	Energy		 Get_H_max_para	(LonSecPos i, LonSecPos j)const;	///< "OUTPUT" function
	Entropy		 Get_S_max_para	(LonSecPos i, LonSecPos j)const;	///< "OUTPUT" function
	Temperature	 Get_Tm_max_para(LonSecPos i, LonSecPos j)const;	///< "OUTPUT" function
	Energy		 Get_G_max_para	(LonSecPos i, LonSecPos j)const 	///< "OUTPUT" function. G at the Ta used for the calculation
	                            {return Get_G_max_para(i,j, last_Ta() );}		 
	Energy		 Get_G_max_para	(LonSecPos i, LonSecPos j, Temperature Ta )const;	///< "OUTPUT" function. G at given Ta

			// --------- Getmaxglo_X	---    

	void		 SelectOptParam_max_para(LonSecPos i, LonSecPos j, Temperature Ta );	///< "OUTPUT" function			// eliges a que Ta

			// ---- SelectOptParam ---    segun step !!!!!-------alternativa a   Get_X	

	void		 SelectOptParam(LonSecPos i, LonSecPos j) {return SelectOptParam(i,j, last_Ta() );}	///< "OUTPUT" function // a la Ta en que se calculo
	void		 SelectOptParam(LonSecPos i, LonSecPos j, Temperature Ta );							///< "OUTPUT" function // eliges a que Ta
protected:										
	       
	        // ---------alternative to   Get_X	----------    da param en (i,j) segun step()i,j) !!!!    ?????????

	Energy		Get_H	(LonSecPos i, LonSecPos j, Step st) const; ///< "OUTPUT" function
	Entropy		Get_S	(LonSecPos i, LonSecPos j, Step st) const; ///< "OUTPUT" function
	Step		Get_pre	(LonSecPos i, LonSecPos j, Step st) const; ///< "OUTPUT" function
	inline Step	step	(LonSecPos i, LonSecPos j		  ) const {return pre(i,j);} ///< "OUTPUT" function

public:										
	
	        // ------ Same that : Get_X , pero calculados todos de una vez. "OUTPUT" x() functions return parametr x for maximized algorithm parametr

	void		 SelectOptParam()				{return SelectOptParam(_maxgloi,_maxgloj);}			// a la Ta en que se calculo
	void		 SelectOptParam(Temperature Ta ){return SelectOptParam(_maxgloi,_maxgloj,Ta);}		
	Energy		 H ()const{return _optH ;}          ///< "OUTPUT" function
	Entropy		 S ()const{return _optS ;}          ///< "OUTPUT" function
	Temperature	 Tm()const{return _optTm;}          ///< "OUTPUT" function
	Energy		 G ()const{return _optG ;}          ///< "OUTPUT" function. G depend on how optParam were selected
	Energy		 G (Temperature Ta )const{return +(_optH-Ta*_optS);  }  ///< "OUTPUT" function G NOT depend on how optParam were selected
	CHit		*GetOptHit();

	void		Export_Hits    (std::ofstream &osHits , char *sep);
	void		Export_DPMz_Tm (std::ofstream &osDP_mz, char *sep);
	void		Export_DPMz_H  (std::ofstream &osDP_mz, char *sep);
	void		Export_DPMz_S  (std::ofstream &osDP_mz, char *sep);
	void		Export_DPMz_Pre(std::ofstream &osDP_mz);

	virtual int	IterationNum()const{return 1;}			//  ??? hace falta aqui?

public:		
	static int const  sti[], stj[], sti1[], stj1[];    /// \todo make non public ?
	std::shared_ptr<CSaltCorrNN>  _NNpar ;
	CSec				*_sd , *_tg ;        
	long				_THits,     ///< total number of hits 
		                _HitsOK ;   ///< total of Ok hits: which are shared with all the the other alignments. \todo: make interface?
	LonSecPos			_maxgloi, _maxgloj	;	///< coordenates of the global max in DP matrix	/*, _maxglot*/ 

protected:	
	Temperature			_Tm_sig,    ///< "significant", cut-off Tm
		                _minTm, 
		                _Ta,
		                _Ta_lastCalc;
	Energy				_G_sig ;    ///< "significant", cut-off G		// _Tm_min, _Tm_max ; // _Tm_min, _Tm_max ; aqui se usan????????????????

	LonSecPos			_LenSond, _LenTarg,  _LenSondPlus1 ;
	long				_TableSize;

	               // Dynamic Programming (DP) Tables for Entropy and Enthalpy  - use vectors ?
	Energy				*_dH0, ///< Enthalpy DP matrix for Enthalpy - for steps that will align x(i+1) with y(j+1) - diagonal, match, mismatch
		                *_dH1, ///< Enthalpy DP matrix for Enthalpy - for steps that will align x(i+1) with a gap in the target  
		                *_dH2; ///< Enthalpy DP matrix for Enthalpy - for steps that will align y(i+1) with a gap in the sonde            
	Entropy				*_dS0, *_dS1, *_dS2;         

	float		_InitMax, _maxglo, _max ;   ///< max of rector parameter  :  like G o Tm

	Entropy		_optS;           ///< optimal value
	Energy		_optH, _optG;    ///< optimal value
	Temperature	_optTm ;         ///< optimal value
	std::vector<CHit> _Hits ;    ///< the list of detected hits.  or better std::list<CHit>?


 protected:
	inline Energy  &dH0(LonSecPos i, LonSecPos j)const{return _dH0[i +j*_LenSondPlus1]; }  //  {return _dH0[i +j*_LenSond]; }   {return *(_dH0 + i +j*_LenSond); }
	inline Energy  &dH1(LonSecPos i, LonSecPos j)const{return _dH1[i +j*_LenSondPlus1]; }
	inline Energy  &dH2(LonSecPos i, LonSecPos j)const{return _dH2[i +j*_LenSondPlus1]; }
	inline Entropy &dS0(LonSecPos i, LonSecPos j)const{return _dS0[i +j*_LenSondPlus1]; }
	inline Entropy &dS1(LonSecPos i, LonSecPos j)const{return _dS1[i +j*_LenSondPlus1]; }
	inline Entropy &dS2(LonSecPos i, LonSecPos j)const{return _dS2[i +j*_LenSondPlus1]; }
	inline Step    &pre(LonSecPos i, LonSecPos j)const{return _pre[i +j*_LenSondPlus1]; }

	Entropy		_restS;		
	Energy		_restH ;
	inline void			RestHS(LonSecPos i)	{	_restS = _sd->_SdS[_sd->Len() -1 ] - (i==1?0:_sd->_SdS[i-2]);
												_restH = _sd->_SdH[_sd->Len() -1 ] - (i==1?0:_sd->_SdH[i-2]);
											;}
};
//Energy		 Getmaxglo_H()				const		{return Get_H (_maxgloi,_maxgloj);}
//Entropy		 Getmaxglo_S()				const		{return Get_S (_maxgloi,_maxgloj);}
//Temperature	 Getmaxglo_Tm()				const		{return Get_Tm(_maxgloi,_maxgloj);}
//Energy		 Getmaxglo_G()				const		{return Get_G (_maxgloi,_maxgloj);}			// a la Ta en que se calculo
//Energy		 Getmaxglo_G(Temperature Ta)const		{return Get_G (_maxgloi,_maxgloj,Ta);}		// eliges a que Ta

/// Hits "inside" of the dynamic programming Matrix, when the algorithm parameter pass the X_sig cut off. 
class CHit 
{	
public: 	
	LonSecPos		_i,_j, _i0, _j0, _l;   ///< begin, end and length of the local alignment hit in both strands
	DNAstrand		_strnd;                ///< ?
	ThDyAlign::Step _Step ;                ///< ?
	Energy			_H,  _G ;
	Entropy			_S ;
	Temperature		_Tm ;
	float			_max ;                 ///< ?
	CHit (LonSecPos i, LonSecPos j, Temperature Tm): _i(i), _j(j), _Tm(Tm){}
	CHit (LonSecPos i, LonSecPos j, Energy H, Entropy S,float max, ThDyAlign::Step st): _i(i), _j(j), _H(H), _S(S), _max(max),_Step(st) {};
	CHit (LonSecPos i, LonSecPos j, LonSecPos i0, LonSecPos j0, LonSecPos l,Energy H, Entropy S,float max, ThDyAlign::Step st)
		: _i(i), _j(j),  _i0(i0), _j0(j0), _l(l)      ,_H(H), _S(S), _max(max),_Step(st) {};
	CHit (ThDyAlign &Al);  ///< save the optimal Hit: recalculate optP of AL at temperature Ta set in NNpar
	explicit CHit (){}
};

class CHitAligned  : public CHit
{
public:
	ISec::sequence          _sd, _tg ;  ///< string with the aligned sequences, (including introduced gaps?)
	std::vector<ThDyAlign::Step> _st ;
	long                    _mt, _mm, _sgap, _tgap ;    // count sonde and target - matchs , mismatch, and gaps
	float		            _Hr, _Sr, _Gr, _Tmr ;

	explicit CHitAligned(ThDyAlign &Al) : CHit(Al) {ExtractAligment(Al);}
	explicit CHitAligned(ISec::sequence s, ISec::sequence t, std::shared_ptr<CSaltCorrNN>  NNpar )  
									: _sd(std::move(s)), 
									  _tg(std::move(t))
	                              {ReCalcule( NNpar );}; 

	void ExtractAligment(ThDyAlign &Al);

	virtual ~CHitAligned()  { }
	void ReCalcule( std::shared_ptr<CSaltCorrNN>  NNpar );
};

class AlignedSecPar : public CHitAligned
{
 public: 
	explicit AlignedSecPar( ISec::sequence s,  ISec::sequence t, std::shared_ptr<CSaltCorrNN>  NNpar ):CHitAligned(s, t, NNpar){};
	Temperature Tm(){return _Tmr;}
	float  G(){return _Gr ;}
		virtual ~AlignedSecPar()  {/*_sd = _tg = nullptr;*/ }
};

//
//ofstream &operator<<(ofstream &stream, CHit &TmACHitl) 
//{	stream	<< endl<< "T melting Align";     //<< (ThDyAlign)TmAl;
//	print_ThDyAlign (stream, TmAl);
//	return stream;
//}

/// the rector parameter for DP is Tm (simplistically) = H / S  (G=H-ST, Tm is T for which G=0)
class ThDyAlign_Tm			: public ThDyAlign  // -----------------------------Tm--------------------ThDyAlign_Tm------
{public:	
	ThDyAlign_Tm ( long MaxLenSond, 
                   long MaxLenTarg, 
                   std::shared_ptr<CSaltCorrNN>  NNpar )
		:ThDyAlign(MaxLenSond, MaxLenTarg, NNpar, NNpar->kein_Tm)
	{
		CalcParam = &ThDyAlign::CalcParamTm;
	} 

	std::string   AlignMeth ()  override   {return "Tm"   ;}
	Temperature	  GetMax_Tm ()  const      {return _maxglo;}
};

/// unused ?
class ThDyAlign_TmHits			: public ThDyAlign_Tm  // -----------------------------Tm---------------ThDyAlign_TmHits----- not in use????------
{public:	
	ThDyAlign_TmHits(   long                    MaxLenSec, 
                      std::shared_ptr<CSaltCorrNN>  NNpar, 
                      float               Tm_min =CtoK(57), 
                      float               Tm_max =CtoK(63) )
		:ThDyAlign_Tm(MaxLenSec, MaxLenSec, NNpar),	 _Tm_min(Tm_min), _Tm_max(Tm_max)
     {} 

	std::string   AlignMeth()  override {return "TmHits";}
	bool	     AddIfHit  (LonSecPos i, LonSecPos j) override;
	Temperature  _Tm_min, _Tm_max ; //aqui se usan????????????????
};

/// Designed specially to find probes that hybridize in both sequences and in others targets too.

/// Take CSecCand instead of CSec, the rector parameter for DP is Tm 
class ThDyAlign_TmCand			: public ThDyAlign_Tm  // -----------------------------Tm-----------------ThDyAlign_TmCand---------
{public:	
	ThDyAlign_TmCand ( long                      MaxLenSec, 
                       std::shared_ptr<CSaltCorrNN>  NNpar)               /*, float Tm_min=CtoK(57), float Tm_max=CtoK(63) */
		:ThDyAlign_Tm(MaxLenSec, MaxLenSec, NNpar){}                            /*,	 _Tm_min(Tm_min), _Tm_max(Tm_max)*/

	std::string   AlignMeth()  override {return "TmCand";}

	virtual bool	AddIfHit	(long i, long j) override;

	void	Use			(CSecCand  *cand1, CSecCand *cand2)	
	                   {	_cs=cand1; 
	                        _ct=cand2; 
							ThDyAlign_Tm::Use	( &_cs->_Sec, &_ct->_Sec);
	                    }

	void	FindCommon	(CSecCand  *cand1, CSecCand *cand2, bool colapse=true) 
	                   {	 Use	(cand1, cand2);
							 Align( &_cs->_Sec, &_ct->_Sec); 
							_cs->ColapseRangs(colapse);
							_ct->ColapseRangs(colapse);		
	                    }
	CSecCand *_cs, *_ct;

	//float  _Tm_min, _Tm_max ; //aqui se usan????????????????
};

/// wrapper for CSec adding Rangs of hit in each position to track existing thermodynamic hits with the other sequences
class CMSecCand  	//--------------------------------Tm------ CMSecCand --------------------------------
{public:
	CMSecCand(	SondeLimits sL ,
				float	Tm_sig		, float G_sig ,			// sonde  - target		: hace falta aqui?
				float	MaxSd_nTgTm , float MinSd_nTgG ,	// sonde  - non target
				float	MaxSelfTm   , float MinSelfG  	)	// sonde 
				:	_sL(sL),
					_Tm_sig(Tm_sig) ,			_G_sig(G_sig),
					_MaxSd_nTgTm(MaxSd_nTgTm) , _MinSd_nTgG(MinSd_nTgG), 
					_MaxSelfTm(MaxSelfTm),		_MinSelfG(MinSelfG),

					_TNumCand(0),	_NumPosCand(0),		_NumCand(0),	_NSecCand (0),
					_TNumPosCand(-1)  /*,		// valor imposible, inicial	
					_osPaarComp(0)	  */	
		{} 

	void		Use(CMultSec* MSec);//	void		Set_PaarComparExport(ofstream &osPaarComp){_osPaarComp=osPaarComp;};
	CSecCand	*Add(CSec &sec);
	std::unique_ptr<CSecCand> AddBeging	(CSec &sec) ;
	void		FindCommon	(CSecCand  &cand1, CSecCand &cand2, bool design=true)	;

	/// Return probes with a percent of other-target coverage with is not intern to the range ExtrCovPerc.

	/// That is: probes which hybrid in one target but in not than more than in ExtrCovPerc.Min % of the others, 
	/// and additionally, probes with hybrid in one target and at last in ExtrCovPerc.Max % of the others.
	void		ExportCommonSonden( bool               colpased, 
		                            NumRang<float>     ExtrCovPerc, 
		                            CMultSec           *res = {},    /// nullptr or a valid sequence group to collect the probe candidates
		                            const std::string  &outpup_fileName ="", 
		                            fileFormat         format =fileFormat::fasta
	                             );

    void write_probes(	CMultSec           *res     ,
		    			const std::string &fileName , 
			    		fileFormat         format   = fileFormat::fasta);
	
	virtual ~CMSecCand(){	 
	                     }
	SondeLimits _sL ;

	float	_Tm_sig, _G_sig ;				// sonde  - target
	float	_MaxSd_nTgTm , _MinSd_nTgG ;	// sonde  - non target
	float	_MaxSelfTm , _MinSelfG  ;		// sonde 

	long	_NumPosCand,	_NumCand;  // Solo los "locales"
	long	_TNumPosCand,	_TNumCand;  // Tambien cuenta los de la multiseq
    long    _NSecCand;

	std::ofstream _osPaarComp;

	CMultSec*            _MSec ;
	std::shared_ptr<ThDyAlign_TmCand>	_TDATmC ;  // donde se crea y se borra???   _______________ PROBLEMA !!!!!!!!!!!!!!!!


	std::list<std::shared_ptr<CSecCand > >	_LSecCand;
	std::list<std::shared_ptr<CMSecCand> >	_LMSecCand;
};

///  the rector parameter for DP is G
class ThDyAlign_G			: public ThDyAlign	// ------------------------------G-------------------------
{public:
	ThDyAlign_G(LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar, Temperature Tm)
		:	ThDyAlign(MaxLenSond, MaxLenTarg, NNpar, NNpar->forbidden_freeEnerg)
																	{CalcParam = &ThDyAlign::CalcParamG;
																	 SetTa (Tm);
																	} 
	ThDyAlign_G(LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar)// usar ultima Ta de calculo en NNpar, o "Set" despues
		:	ThDyAlign(MaxLenSond, MaxLenTarg, NNpar, NNpar->forbidden_freeEnerg)
																	{CalcParam = &ThDyAlign::CalcParamG;} 	
	std::string AlignMeth()  override {return "G";}
	float		 GetMax_G()const{return -_maxglo;} // mas o menos lo mismo, pero "comprobado", usa ultima Ta de calculo
};

/// fractional DP
class FracTDAlign	: public ThDyAlign			// --------------------------Frac-----------------------------
{	bool	_finisch ;
	int		_iterations, _maxNumIt, _fixedNumIter; //, _totalIterations 
	Energy	_maxG_der ;  // epsilon 
 public:
	 std::string AlignMeth()  override {return "FractG";}
	 FracTDAlign(LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar)  
		 :	ThDyAlign	(MaxLenSond, MaxLenTarg, NNpar),
																	_iterations(0),	// prohibido comenzar sin BeginAlign
																	_maxG_der(1.0f),// super excesivo
																	_maxNumIt(10),	// super excesivo
																	_finisch(true),	// prohibido comenzar sin BeginAlign
																	_fixedNumIter(0)// NOT fixed !!	
																{}
	 void			SetmaxG_der	(Energy maxG_der) { _maxG_der = maxG_der;}	// max desviacion permisible de la G sobre 0 o de G(Tm)
	 void			SetEpsilum  (Energy maxG_der) {  SetmaxG_der( maxG_der);} // exactamente lo mismo
	 void			SetMaxNumIt	(int maxNumIt)	 { _maxNumIt = maxNumIt;}
	 int			SetFixedNumIt(int fixedNumIter){int fixedNum = _fixedNumIter;_fixedNumIter=fixedNumIter;return fixedNum;}

	 void			BeginAlign	(CSec  *sonde, CSec *target)	;
	 bool			NotFinisch	();
	 void			iterate		(Temperature ta);
	 inline void	iterate		() {iterate ( Tm() );}

	 virtual void	Align (CSec  *sonde, CSec *target)	{for ( BeginAlign(sonde, target) ; NotFinisch(); iterate() );}


//	 float			G			()const{return _optG  ;}// ------ usa la Ta usada para los optParam
	 float			GetMax_G	()const{return  IterationNum()>1 ? -_maxglo : _maxglo;}// mas o menos lo mismo, pero usa Ta de calculo.	CUIDADO : Tm o G . G solo si interations>1  !!!!!
	 virtual int	IterationNum()const{return _iterations;}
};

class ThDyAlign_restTm			: public ThDyAlign  // ---------------------------HACER----------------------------
{public:	
	ThDyAlign_restTm(long MaxLenSond, long MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar, Temperature minTm)
			:		ThDyAlign(MaxLenSond, MaxLenTarg, NNpar, -274)
			{		_minTm = minTm ; /*CalcParam = &ThDyAlign::CalcParamRs;*/} 
	std::string  AlignMeth()  override {return "restTm";}
	virtual inline float CalcParam (float S, float H) const		
	{	
		Temperature Tm = _NNpar->CalcTM (S +_restS, H +_restH);

		if ( _NNpar->CalcTM (S , H )==0 ) Tm= 0;
		//if ( Tm < _minTm ) return 0;
		Tm -= _sd->_Tm.Ave();
		return 		Tm ;
	}  
	float	GetMax_Tm()const{return _maxglo;}
};





void		print_ThDyAlign (std::ofstream &osTm,ThDyAlign &Al);
//void		print_Tm		(ofstream &osTm, CMultSec	&pr, int MaxGrDeg=-1, char sep[]=";" );
void print_Tm (std::ofstream &osTm, CMultSec	&pr, int MaxGrDeg, char sep[]);
std::ofstream	&operator<<(std::ofstream &stream,	ThDyAlign_Tm	&TmAl) ;
//std::ofstream	&operator<<(std::ofstream &osTm,	ThDyAlign		&Al) ;
std::ofstream	&operator<<(std::ofstream &stream,	ThDyAlign_G		&G_Al) ;
std::ofstream	&operator<<(std::ofstream &stream,	FracTDAlign		&FrAl) ;
#endif

