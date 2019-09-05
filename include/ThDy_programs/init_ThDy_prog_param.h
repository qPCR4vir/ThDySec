/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @author Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2019
*
* @file  ThDySec/include/ThDy_programs/init_ThDy_prog_param.h
*
* @brief 
*
*/

#ifndef _INIT_ThDy_PROG_PARAM_H
#define _INIT_ThDy_PROG_PARAM_H

#include <memory>
#include <set>

#  include <nana/filesystem/filesystem_ext.hpp>   // std filesystem ?

#include <init_prog_param.h>       // from ../ProgParam/include
#include <Units.hpp>
#include "ThDySec/sec_mult.h"


const SecPos MAX_SEQ_LEN_ALIGN{ 12501 };     // ?

class CParamSondeLimits: public IBParam  /// \todo Use delegating constructor to limit code duplication
{   SondeLimits sL;
 public: 
    CParamNumMinMax<Energy>      G;
    CParamNumMinMax<Temperature> T; 
    CParamNumMinMax<SecPos>      L;
	 /// Acepta un parametro y por tanto no usa _v. Por compatibilidad.
    CParamSondeLimits (IProg *pp, const std::string& titel, SondeLimits &parRef, 
		             const std::string& etiqGmin, Energy minGmin,        Energy maxGmin,         Energy defValueGmin,
		             const std::string& etiqGmax, Energy minGmax,        Energy maxGmax,         Energy defValueGmax,
		             const std::string& etiqTmin, Temperature minTmin,   Temperature maxTmin,    Temperature defValueTmin,
		             const std::string& etiqTmax, Temperature minTmax,   Temperature maxTmax,    Temperature defValueTmax,
		             const std::string& etiqLmin, SecPos minLmin,        SecPos maxLmin,         SecPos defValueLmin,
		             const std::string& etiqLmax, SecPos minLmax,        SecPos maxLmax,         SecPos defValueLmax,
                           std::string  UnitG="kcal/mol", std::string UnitT= RTunits::grC, std::string UnitL="nt"

					) : IBParam ( titel), 
					    G(pp, titel+". Perf.Match dG", parRef._G,
						    " Min", etiqGmin,  minGmin, maxGmin, defValueGmin, 
						    " Max", etiqGmax,  minGmax, maxGmax, defValueGmax, 
						    UnitG),
					    T(pp, titel+". Perf.Match Tm", parRef._Tm,
						    " Min", etiqTmin,  minTmin, maxTmin, defValueTmin, 
						    " Max", etiqTmax,  minTmax, maxTmax, defValueTmax, 
						    UnitT),
					    L(pp, titel+". Length", parRef._L,
						    " Min", etiqLmin,  minLmin, maxLmin, defValueLmin, 
						    " Max", etiqLmax,  minLmax, maxLmax, defValueLmax, 
						    UnitL)
	          { 
	          } 
	        /// Num &parRef,   usa _v y por tanto no necesita un parametro externo
    CParamSondeLimits (IProg *pp, const std::string& titel, //SondeLimits &parRef, 
		             const std::string& etiqGmin, Energy minGmin,        Energy maxGmin,         Energy defValueGmin,
		             const std::string& etiqGmax, Energy minGmax,        Energy maxGmax,         Energy defValueGmax,
		             const std::string& etiqTmin, Temperature minTmin,   Temperature maxTmin,    Temperature defValueTmin,
		             const std::string& etiqTmax, Temperature minTmax,   Temperature maxTmax,    Temperature defValueTmax,
		             const std::string& etiqLmin, SecPos minLmin,        SecPos maxLmin,         SecPos defValueLmin,
		             const std::string& etiqLmax, SecPos minLmax,        SecPos maxLmax,         SecPos defValueLmax,
                     const std::string& UnitG="kcal/mol", const std::string& UnitT= RTunits::grC, const std::string& UnitL="nt"

					) : IBParam ( titel), 
					    G(pp, titel+". Perf.Match dG", sL._G,
						    " Min", etiqGmin,  minGmin, maxGmin, defValueGmin, 
						    " Max", etiqGmax,  minGmax, maxGmax, defValueGmax, 
						    UnitG),
					    T(pp, titel+". Perf.Match Tm", sL._Tm,
						    " Min", etiqTmin,  minTmin, maxTmin, defValueTmin, 
						    " Max", etiqTmax,  minTmax, maxTmax, defValueTmax, 
						    UnitT),
					    L(pp, titel+". Length", sL._L,
						    " Min", etiqLmin,  minLmin, maxLmin, defValueLmin, 
						    " Max", etiqLmax,  minLmax, maxLmax, defValueLmax, 
						    UnitL)
	          { 
	          } 
};

class CMultSec;
class CSaltCorrNN;

/// concreta los parametros comunes. Mantiene lista de los prog Espec que los usan
class ThDyCommProgParam : public CCommProgParam 
{public:	
    CParamString     _OutputFile      {this, "Results output file",                                  "OutputFile",     "" } ;///\todo make it a std::path

    CParamString     _InputTargetFile {this, "Input file for Targets",				                 "TargetFile",     "" } ;
    CParamBool       _TRecurDir       {this, "Recursively add all Target seq-files from all dir",    "TRecursDir",    false} ;
    CParamBool       _TDirStrOnly     {this, "Reproduce only the dir struct in targets"         ,    "TOlyDirStr",    true } ;

    CParamString     _NonTargetFile   {this, "Input file for non-Targets",			                  "NonTargetF",    "" } ;
    CParamBool       _nTRecurDir      {this, "Recursively add all non-Target seq-files from all dir",   "nTRecurDir",   false} ;
    CParamBool       _nTDirStrOnly    {this, "Reproduce only the dir struct in non-targets"         ,   "nTlyDirStr",   true } ;

    CParamString     _PCRfiltrPrFile  {this, "Input file with primers for filtering",                        "PCRftrFile", "" } ;
    CParamBool       _FiltrRecuDir    {this, "Recursively add all filter PCR primer seq-files from all dir", "FltrRecDir", false} ;
    CParamBool       _FiltrStrOnly    {this, "Reproduce only the dir struct in filter PCR primer"     ,        "FltrStrOly", true } ;
   
    CParamString     _InputNNFile     {this, "Input file with NN parameters",		                       "iNNParFile",   ""	 } ;

	SaltCorrection                 _SaltCorr {StLucia }  ;			
    CParamEnumRange<SaltCorrection>	SaltCorr {this, "Salt correction method",		  "SaltCorrMt", _SaltCorr, StLucia, StLucia, StLucia }  ;

	float					_ConcSd,	_ConcTg,	_ConcSalt ;
	CParamNumRange<float>	 ConcSd  {this,    "Concentration of the probes" ,				"ConcSonden", _ConcSd,  0.0f,50e-3f,  0.8e-6f,"M" } ,	 
                             ConcTg  {this,    "Concentration of the targets",				"ConcTarget", _ConcTg,  0.0f,50e-3f,  1e-9f  ,"M" } ,	 
                             ConcSalt{this,    "Concentration of salt",						"ConcenSalt", _ConcSalt,0.0f,1.0f,    50e-3f ,"M" } ;

	Temperature		               _Ta  {55.0f }  ;				
    CParamNumRange <Temperature>	Ta  {this, "Temp annealing expected in exp",		"TempAnnelg", _Ta,  20.0f, 90.0f, 55.0f, RTunits::grC  }  ;
	AlignMeth		           _TAMeth  {TAMeth_Tm }  ;			
    CParamEnumRange<AlignMeth>	TAMeth  {this, "Optimized parameter during DynProg",	"AlignMethd", _TAMeth, TAMeth_Tm, TAMeth_Fract, TAMeth_Tm }  ;	

    ///  Sequence import filtres
	float			           _MaxTgId {99.0f }  ;			
    CParamNumRange <float>	    MaxTgId {this, "Filter Seq Target with more %id"  ,	"MaxTgIdent",  _MaxTgId,  0.0f,100.0f,  99.0f  ,"%" }  ;	
	LonSecPosRang	           _SecLim  {1, 0 } ;			
    CParamNumMinMax<LonSecPos>  SecLim {this, "Fragments of seq. to use",	_SecLim, 
													"Begin",						"SecBeginng", 1, std::numeric_limits<LonSecPos>::max(), 1 ,
													"End",							"SecuencEnd", 0, std::numeric_limits<LonSecPos>::max(), 0 ,		"nt" } ;	
	SecPosRang			       _SecLenLim {1,0 } ;			
    CParamNumMinMax<SecPos>     SecLenLim {this, "Limits of the length of the seq. to be useful",_SecLenLim ,
													"Minimal length",				"MinSecuLen", 15, MAX_SEQ_LEN_ALIGN, 20 ,
													"Maximal length",				"MaxSecuLen", 15, MAX_SEQ_LEN_ALIGN, MAX_SEQ_LEN_ALIGN ,		"nt" } ;	
	///  Program setup bool parameters
    CParamBool       loadNNPar  {this, "Programs option- Load NN parameter",		"LoadNNPara",    false } ,	     
                     saveNNPar  {this, "Programs option- save NN parameter",		"SaveNNPara",    false }  ; 
	CParamBool		 st_savTm   {this, "Programs option- Save Tm Table",		"SavTmTable",    false },
		             st_savPos  {this, "Program option- Save Position Table",	"SavPosTabl",    false },
					 st_savG    {this, "Program option- Save free Energy Table","SaveGTable",   false },
					 st_savAlign{this,"Program option- Save Table of aligns",   "SavAlignTb",   false },
					 st_savProj {this,"Program option- Save Project file",      "SavProject",   true  }, 
					 st_savG_Plasm{this, "Program option- Gtable for Plasmid",  "SavG_Plasm",   false }, 
					 st_savTm_Plasm{this,"Program option- Tm table for Plasmid","SavTmPlasm",   false }, 
					 st_savLog  {this, "Program option- save a log",			 "Save_Logfi",   false }, 
					 st_Exp_sond{this, "Program option- re-export probes",	     "Exp_probes",   false }, 
					 st_ExpTarg {this, "Program option- re-export targets",	 "Exp_target",   false };

   	std::set<CMultSec*>           _primersGr  ;  /// the real CMultiSec are in the main tree. Here just for convenience.
    std::shared_ptr<CSaltCorrNN>  _pSaltCorrNNp     {Create_NNpar( )};
    /// In this tree will be ALL other CMultiSec. The next pointers are just for convenience.
	std::shared_ptr<CMultSec>     _pSeqTree      = CreateRoot()  ;  ///<  The super-root of the sequences tree
	CMultSec*     _pSeqNoUsed       {(AddSeqGroup(*_pSeqTree, "Don t use"  ))->get() } ; ///< the root for deleted
	CMultSec*     _pSeqTargets      {(AddSeqGroup(*_pSeqTree, "Target seq"))->get() } ;  ///< the root for Viral seq
	CMultSec*     _pSeqNonTargets   {(AddSeqGroup(*_pSeqTree, "Non Target seq"))->get()  } ; ///<
	CMultSec*     _pPCRfiltrePrimers{(AddPrimerGroup(*_pSeqTree, "PCR Primers to <filtre> sequences"))->get()} ;

	///
	void Actualize_All_NNp()
	{
        if ( ! _pSaltCorrNNp->NeedActualization(_ConcSd, _ConcTg, _ConcSalt, _SaltCorr)  )
			return;
        _pSaltCorrNNp=Create_NNpar ();
		Actualize_NNp_recur(_pSeqTree.get());
	}
	
	///
	void Actualize_NNp_recur(CMultSec *ms)   // ?
		{
			Actualize_NNp(ms);
            for (auto &CurMSec : ms->MSecL() ) 
                Actualize_NNp_recur(CurMSec.get());
		}
    
	///
	void Actualize_NNp(CMultSec*ms)
    {
        ms->_NNPar=_pSaltCorrNNp;
		//for (auto &CurMSec : ms->MSecL())
		//	Actualize_NNp_recur(CurMSec.get());
	}

	///
    ThDyCommProgParam(const std::string& titel,   CProject *proj)
		:	CCommProgParam(titel,proj)
		{	
			TAMeth.AddStrValues("TAMeth_Tm",	TAMeth_Tm);
			TAMeth.AddStrValues("TAMeth_G",		TAMeth_G);
			TAMeth.AddStrValues("TAMeth_Fract",	TAMeth_Fract);
		
			SaltCorr.AddStrValues("NoSelect",	NoSelect);
			SaltCorr.AddStrValues("StLucia",	StLucia);
			SaltCorr.AddStrValues("Owczarzy",	Owczarzy);
		} 

	///< Create a new set of NN parameters based on current concentration. and Ta set in com-par.
    std::unique_ptr<CSaltCorrNN> Create_NNpar        ( )
    {
        std::unique_ptr<CSaltCorrNN> NNpar ( new  	CSaltCorrNN	(	 _ConcSd,  _ConcTg,  _ConcSalt, _SaltCorr, loadNNPar.get()?_InputNNFile.get():""));
        //auto NNpar = std::make_unique<CSaltCorrNN>();
        NNpar->SetTa(	CtoK(	 _Ta));
        return NNpar ;
    }
    std::unique_ptr<CSaltCorrNN> Init_NNpar          ()   //< Initialize the set of NeirN parametrs  in com-par, loading/saving if necesary
    {                                      /// Depend on  _ConcSd,  _ConcTg,  _ConcSalt , _Ta, _loadNNPar,  _InputNNFile,  _saveNNPar,_OutputFile
        std::unique_ptr<CSaltCorrNN> NNpar=Create_NNpar();
        if (loadNNPar.get())
        {
            std::ifstream isTDP(_InputNNFile.get());	assert(isTDP);
            NNpar->LoadNNParam(isTDP) ;
        }
        if ( saveNNPar.get())
        {
            //std::string OutputTDP( _OutputFile.Get()) ; OutputTDP += ".ThDyParam.csv";

            std::ofstream osTDP	(_OutputFile.get() + ".ThDyParam.csv");				assert(osTDP);
            osTDP << *NNpar ;
        }
        return NNpar ;
    }

    std::shared_ptr<CSaltCorrNN> Get_Actualiced_NNpar(const std::shared_ptr<CSaltCorrNN>& currNNpar )
    {
        if ( currNNpar)
        {
            if ( currNNpar->NeedActualization(_ConcSd, _ConcTg, _ConcSalt, _SaltCorr)  )
                return Create_NNpar ();

            currNNpar->SetTa(	CtoK(_Ta));
            return currNNpar ;
        }
        if ( _pSaltCorrNNp )
            return Get_Actualiced_NNpar(_pSaltCorrNNp);
        return Init_NNpar ();
    }

    void                    Actualice_NNp       ()
    {
        _pSaltCorrNNp=Get_Actualiced_NNpar(_pSaltCorrNNp);
    }

    void Check_NNp_Targets (/*ThDyCommProgParam& cp*/)
    {
        if (! _pSaltCorrNNp )
            Actualice_NNp();

        assert( _pSeqTargets );
        if ( !  _pSeqTargets->_Global._NSec)
                AddSeqFromFile (    *_pSeqTargets,
                                    _InputTargetFile.get()  );
    }

    CMultSec::pMSec CreateRoot	();

    CMultSec::MSecIt AddSeqGroup(CMultSec &parentGr, const std::string& Name);

    /// used to create only base, sub-root primer groups.
    CMultSec::MSecIt AddPrimerGroup(CMultSec &parentGr, const std::string& Name)
    {
        CMultSec::MSecIt pgr = AddSeqGroup(parentGr, Name);
        _primersGr.insert((pgr)->get());
        return pgr;
    }


	/// take parameters from the parent if possible
    CMultSec::MSecIt AddSeqFromFile(CMultSec &parentGr,
                                    const std::filesystem::path &FileName,
                                    bool recursive = false,
                                    bool onlyStructure = false);

    CMultSec::MSecIt CopyStructFromDir(CMultSec &parentGr, const std::string& FileName)
    {
        return AddSeqFromFile(parentGr, FileName, true, true);
    }

    void      LoadSequences ()
    {
        if (!_InputTargetFile.get().empty())   AddSeqFromFile(*_pSeqTargets       , _InputTargetFile.get(), _TRecurDir   .get() , _TDirStrOnly .get() );
        if (!_NonTargetFile  .get().empty())   AddSeqFromFile(*_pSeqNonTargets    , _NonTargetFile  .get(), _nTRecurDir  .get() , _nTDirStrOnly.get() );
        if (!_PCRfiltrPrFile .get().empty())   AddSeqFromFile(*_pPCRfiltrePrimers , _PCRfiltrPrFile .get(), _FiltrRecuDir.get() , _FiltrStrOnly.get() );
    }
    void      CopyStructFromDir ()
    {
        if (!_InputTargetFile.get().empty())   CopyStructFromDir(*_pSeqTargets      , _InputTargetFile.get());
        if (!_NonTargetFile  .get().empty())   CopyStructFromDir(*_pSeqNonTargets   , _NonTargetFile  .get());
        if (!_PCRfiltrPrFile .get().empty())   CopyStructFromDir(*_pPCRfiltrePrimers, _PCRfiltrPrFile .get());
    }

    CMultSec *AddTargetFromFile(const std::string& FileName)                      ///< deprecated   !!!
	{
		return AddSeqFromFile(*_pSeqTargets,FileName)->get();
	}
	CMultSec *AddNoTargetFromFile(const std::string& FileName)                    ///< deprecated   !!!
	{
		return AddSeqFromFile(*_pSeqNonTargets,FileName)->get();
	}
	CMultSec *AddPCRfiltreFromFile(const std::string& FileName)                   ///< deprecated   !!!
	{
		return AddSeqFromFile(*_pPCRfiltrePrimers,FileName)->get();
	}


	// convertirlas en funciones "previas a la paralelizacion", que hacen copias propias de los parametros en serie, no en paralelo ??
	void    OutputFile(const std::string &OutputFile)     {	_OutputFile     .set( trim_string(OutputFile     ))   ;	}
	void    TargetFile(const std::string &InputTargetFile){	_InputTargetFile.set( trim_string(InputTargetFile));	}
	void    NNParaFile(const std::string &InputNNFile)	  {	_InputNNFile    .set( trim_string(InputNNFile    ))	;	}
	
	~ThDyCommProgParam(void) override//;	
    {	/*delete []_ProgList;*/
       //_pSeqTree->Free();     // ?????????
    }
};

class CProgParam_microArray ;
int microArrayProg   ( CProgParam_microArray   *IPrgPar_uArr )  ;

template <typename Num, typename RowTitle> class CTable ;          //typedef uArr_RT CTable<TmGPos> ;
using Table = CTable<TmGPos, std::shared_ptr<CSec    >>;

class CEspThDyProgParam : public CEspProg
{public:
	ThDyCommProgParam& _cp;
	CEspThDyProgParam(const std::string& titel, ThDyCommProgParam &commThDyParam)
		: CEspProg(titel, commThDyParam),
		  _cp(commThDyParam)
	{}
	~CEspThDyProgParam()override{}
};


class CProgParam_microArray : public CEspThDyProgParam
{public:	
    CMultSec::pMSec _probesMS{*_cp.AddPrimerGroup(*_cp._pSeqTree, "Probes for Virtual uArr")};
    CParamString	            _InputSondeFile{ this, "Input file for probes", "iSonde_uAr", "" };
    CParamBool       _PrRecurDir    {this, "Recursively add all probe seq-files from all dir", "ProbRecDir", false} ;
    CParamBool       _PrDirStrOnly  {this, "Reproduce only the dir struct in probe"          , "ProbDirStr", true } ;
      
    std::shared_ptr<Table> _rtbl ;		                //uArr_RT *_rtbl;

	//bool			    _I, _G;			// Outpu table of I, G. 
	//CParamBool		 I,  G;			// Outpu table of I, G. 

    CProgParam_microArray(const std::string& titel, ThDyCommProgParam &commThDyParam) 
		    :	CEspThDyProgParam(titel, commThDyParam)
	        {		
			} 

	void      SondeFile(const std::string &InputSondeFile)	  {	_InputSondeFile .set( trim_string(InputSondeFile    ))	;	}
    void RenameSondesMS(const std::string& name);

    void      LoadSequences ()
    {
        if (!_InputSondeFile.get().empty())
            _cp.AddSeqFromFile(*_probesMS, _InputSondeFile.get(), _PrRecurDir   .get() , _PrDirStrOnly .get() );
    }

    void CopyStructFromDir ()                                    ///< deprecated   !!!
    {
        if (!_InputSondeFile.get().empty())   
			_cp.CopyStructFromDir(*_probesMS, _InputSondeFile.get());
    }

    CMultSec *AdduArrFromFile(const std::string& FileName)       ///< deprecated   !!!
	{
		return _cp.AddSeqFromFile(*_probesMS,FileName)->get();
    }

    void Check_NNp_Targets_probes      (CMultSec *probes) 
{
	_cp.Check_NNp_Targets ();
    assert(("Trying to load probe seq into inexistent MultiSec", probes));
    if (! probes->_Global._NSec)
	 if (!_InputSondeFile.get().empty())
		probes->AddFromFile ( _InputSondeFile.get() );	
     else assert((std::cerr<<"No seq. in the probes MS in uArr routine.", true)); //false??
}
    

	virtual int Run (){	return microArrayProg ( this )  ;}

	// cuando se corre un proceso paralelo ver donde es mejor hacer estos delete.
	~CProgParam_microArray()	override	;
};
 
//,_tlG(nullptr),_tlPos(nullptr), UpDate(nullptr)
//				_I(false),   I				(this, "Programm option- Save Tm Table",			"SavTmTable", _I,   false), 
//				  _G(true),    G				(this, "Programm option- Save Tm Table",			"SavTmTable", _G,   true),  
//class	CProgParam_uArrExp;
//int		microArrayProgTest ( CProgParam_uArrExp *IPrgPar_uArr)  ;
//class CProgParam_uArrExp  : public CProgParam_microArray
//{public:	
//	C_str		_Input_uArrExpFile ; 
//	bool		_exclSd;					// No analizar (excluir) determinadas sondas (otros org - micobact, contaminadas, Post contr, Neg contr, etc)
//	C_str		_exclSdMark;				// solo las sondas con esta marca son excl. Si ="" todas las marcadas. Si _exclSd -  se ignora, y nada se excl.
//	bool		_IxI,		_IxI_d;			// output these table (metodo de comparacion - solo similitudes o similitudes menos diferencias?)
//	bool		_Normalize;					// las tablas anteriores : puede suponer perdida de informacion para calcular signif estadistica??
//	Energy		_Isat, _Isen, _Gsat, _Gsen;
//
//    explicit 	CProgParam_uArrExp (const std::string& titel, ThDyCommProgParam &commThDyParam);
//
//	int		Run		(	){	
//                            Check_NNp_Targets (/*IPrgPar_SdDes->*/_cp);
//                            Check_probes      (/*IPrgPar_SdDes->*/_probesMS);
//                            return microArrayProg ( this )  ;
//                         }
//	void    ExpFile (const char *Input_uArrExpFile )	{	_Input_uArrExpFile.CopyTrim(Input_uArrExpFile) ;	}
//
//	// OJO !!!!!!!!!   la sig funcion se aduena del pointer, y luego lo deletea   !!!!!!!!
//	void SetExpFile (      char *Input_uArrExpFile )	{	ExpFile( Input_uArrExpFile) ;	delete []Input_uArrExpFile   ;	}
//
//	// cuando se corre un proceso paralelo ver donde es mejor hacer estos delete.
//	virtual ~CProgParam_uArrExp()		{}
//};

class CProgParam_MultiplexPCR ;
int MultiplexPCRProg ( CProgParam_MultiplexPCR *IPrgPar_uArr )  ;

class CProgParam_MultiplexPCR : public CProgParam_microArray
{public:
	CProgParam_MultiplexPCR(const std::string& titel, ThDyCommProgParam &commThDyParam) ;
	int		Run		()   {	
                            return MultiplexPCRProg ( this )  ;
                         }
    std::shared_ptr<Table> _rtbl_self ;		                //uArr_RT *_rtbl;
};

class CProgParam_SondeDesign ;
int SondeDesignProg  ( CProgParam_SondeDesign  &IPrgPar_SdDes)  ;

class CProgParam_SondeDesign : public CEspThDyProgParam			//  .------------------------	CProgParam_SondeDesign	----------------
{public:

	CProgParam_SondeDesign(const std::string& titel,ThDyCommProgParam &commThDyParam)
      : CEspThDyProgParam (titel, commThDyParam),
		_design (true), design(this, "Make only design or full comp",	        "DesigVsCmp", _design,  true),
		                unique(this, "Find probes unique for a few sequences",	"FindUnique",           true),
		                common(this, "Find probes common for most sequences",	"FindCommon",           true),
        sL(this, "Probe design limits", _sL,							 
		        "MinSonde_G",  -100.0f,   100.0f,  -5.0f,
		        "MaxSonde_G",  -100.0f,   100.0f,  -1.0f,
		        "MinSondeTm",  -270.0f,    90.0f,	57.0f,
		        "MaxSondeTm",  -270.0f,    99.0f,	63.0f,
		        "MinSondeLn",  3,         200,		20,
		        "MaxSondeLn",  4,         200,		35,
                /*UnitG=*/"kcal/mol",  /*UnitT=*/RTunits::grC,  /*UnitL=*/"nt"
		   ),
		_G_sig (10),		  G_sig  (this, "Significant G probe-target",		"MaxSd_Tg_G",    _G_sig,  -15.0f, 30.0f,	10.0f, "kcal/mol" ),	
		_Tm_sig (30)	, 	  Tm_sig (this, "Significant Tm probe-target",	"MinSd_TgTm",    _Tm_sig, -0.0f,  80.0f,	30.0f, RTunits::grC ),

		_MinSd_nTgG (15) , MinSd_nTgG(this, "Significant G probe-non target",	"MinSdnTg_G", _MinSd_nTgG,  0.0f, 30.0f,	15.0f, "kcal/mol" ),		
		_MaxSd_nTgTm (10),MaxSd_nTgTm(this, "Significant Tm probe-non target",	"MaxSdnTgTm",_MaxSd_nTgTm, -0.0f,  70.0f,	10.0f, RTunits::grC ),

		_MinSelfG (10),		 MinSelfG(this, "Significant self probe G",		"MinSdSlf_G", _MinSelfG,  0.0f, 30.0f,		10.0f, "kcal/mol" ),		
		_MaxSelfTm (10),    MaxSelfTm(this, "Significant self probe Tm",		"MaxSdSlfTm",_MaxSelfTm, -0.0f,  70.0f,		10.0f, RTunits::grC )
		//_MinTgCov (100),	MinTgCov (this, "Find probes with more % coverage",	"Min_Tg_Cov",  _MinTgCov,  0.0f,100.0f,		99.0f  ,"%")							
        {
	    }  // revisar cuales deben ser estos valores !!!!	

    bool		_design ;  // realizar solo diseno de sondas o solo comparacion de sec????
	CParamBool	design, unique, common ;  // realizar solo diseno de sondas o solo comparacion de sec????

	SondeLimits _sL ;	
	CParamSondeLimits sL;                                    

	Temperature	               _Tm_sig ;	    Energy		            _G_sig ;				// sonde  - target
	Temperature	               _MaxSd_nTgTm ;   Energy		            _MinSd_nTgG ;			// sonde  - non target
	Temperature	               _MaxSelfTm ;	    Energy		            _MinSelfG  ;			// sonde 
	CParamNumRange<Temperature>	Tm_sig ;	    CParamNumRange<Energy>   G_sig ;				// sonde  - target
	CParamNumRange<Temperature>	MaxSd_nTgTm ;	CParamNumRange<Energy>   MinSd_nTgG ;			// sonde  - non target
	CParamNumRange<Temperature>	MaxSelfTm ;		CParamNumRange<Energy>   MinSelfG  ;			// sonde 

	//float		           _MinTgCov ;
	//CParamNumRange<float>	MinTgCov ;

    CParamNumMinMax<float> Coverage {this, "Find probes with % of target coverage",  
											" less than of the others (unique probes)", "PercUnique",  0.0f,  100.0f,   0.0f, 
											" at last of the others (common probes)"  , "PercCommon",  0.0f,  100.0f, 100.0f, 
											"%"}; 


	int		Run		(){	
                        //Check_NNp_Targets (/*IPrgPar_SdDes->*/_cp);
                        return  SondeDesignProg( *this )  ;
                      }
	/// results
	CMultSec::pMSec probes;

	struct targets_comp
	{
		//std::string title;
		static std::vector<std::string> headers;
		long         iteration_num, target_num ;
		std::string  target_1_name,  target_2_name;

		struct comp
		{
			long         t_n_pos, t_n_cand;
			long         target_1_n_cand_pos, target_1_n_cand;
			long         target_2_n_cand_pos, target_2_n_cand;
		} before, after;

		long THits, HitsOK;
	};

	std::vector<targets_comp> targets_comparitions;
};


class CProgParam_TmCalc ;
int DegTmCalc ( CProgParam_TmCalc *IPrgPar)  ;  

class CProgParam_TmCalc : public CProgParam_MultiplexPCR
{public:
	CParamBool		 save {this, "Save result to file"         , "TmCsaveRes",  false},  
		             align{this, "Align primers before Tm calc", "TmUseAlign",  true } ; //  save results as mPCR, using the TargetSecFile and ResultFile comunes. Use ThDyAlign or just calc Tm?
	CParamString    _Sec       {this, "Primer",						"TmCalc_Sec", ""} ,			
		            _Sec2Align {this, "Primer to align",			"TmC_Sec2Al", ""};

	TemperatureRang _TmS, _Tm2A , _TmHy ;	// Para recoger los reslutados y display it back. 
	EnergyRang		_GS , _G2A  , _GHy  ; 
	std::string     _AlignedSec,	_AlignedSec2Align ;

	CProgParam_TmCalc (const std::string& titel,ThDyCommProgParam &commThDyParam) :	
						CProgParam_MultiplexPCR (titel,commThDyParam)
			{
                _InputSondeFile.SetTitel("Imput oligos for TmCalc"); 
		        _InputSondeFile.SetEtiq ("iOligo_TmC", this); 
				RenameSondesMS("Tm calulator sondes??");

                _PrRecurDir.SetTitel("Recursively add all oligos seq-files from all dir"); 
		        _PrRecurDir.SetEtiq("OligRecDir", this); 

                _PrDirStrOnly.SetTitel("Reproduce only the dir struct in Oligos"); 
		        _PrDirStrOnly.SetEtiq("OligDirStr", this); 
	        }

	void	Set_Sec				 (const std::string& Sec){_Sec      .set (Sec)	;	  }
	void	Set_Sec2Align		 (const std::string& Sec){_Sec2Align.set (Sec)	;	  }
	void	Set_AlignedSec		 (const std::string& Sec){_AlignedSec      =Sec	;	  }
	void	Set_AlignedSec2Align (const std::string& Sec){_AlignedSec2Align=Sec ;     }

	void	Update_Sec			(bool rev, bool complem)	{ Set_Sec      ( Generate_DegSec_char( _Sec.get().c_str(),	     rev, complem)  ); }
	void	Update_Sec_Sec2Align(bool rev, bool complem)	{ Set_Sec2Align( Generate_DegSec_char( _Sec.get().c_str(),	     rev, complem)  ); }
	void	Update_Sec2Align	(bool rev, bool complem)	{ Set_Sec2Align( Generate_DegSec_char( _Sec2Align.get().c_str(), rev, complem)  ); }
	void	Update_Sec2Align_Sec(bool rev, bool complem)	{ Set_Sec	   ( Generate_DegSec_char( _Sec2Align.get().c_str(), rev, complem)  ); }

	~CProgParam_TmCalc(){}
	int		Run		()
    {	
        //if (_save)   Check_NNp_Targets ( _cp);
        return DegTmCalc ( this )  ;
    }
};


class ThDyProject : public CProject /// Permite manejar todo el projecto: con un miembro para los parametros comunes y otro para los de cada programa
{public:
		ThDyCommProgParam		_cp   {"Common parameters for all functions",this}  ;
		CProgParam_microArray   _uArr {"Virtual microarray experiment"     ,_cp }  ;
		CProgParam_MultiplexPCR _mPCR {"Check multiplex PCR"               ,_cp }  ;
		CProgParam_SondeDesign	_SdDes{"Find probes"                       ,_cp }  ;
		CProgParam_TmCalc		_TmCal{"Tm calculator"                     ,_cp }  ;

        void LoadSequences()
        {
            _cp  .LoadSequences();
            _uArr.LoadSequences();
            _mPCR.LoadSequences();
            _TmCal.LoadSequences();
        }
        void ExportFASTA(CMultSec* ms, bool only_selected)
        {
            if (ms) ms->Export_from   (*_cp._pSeqTree,  only_selected);

                    //nana::filebox  fb{ *this, true };
                    //fb .add_filter ( SetupPage::FastaFiltre( )                   )
                    //   .title(("Directory load: Add a tree of groups of sequences from a directory"));
                    //if (fb()) 
                    //    AddMSeqFiles(nana::charset(fb.file()), true);

        }
	    int	Run	(IProg &prog ) override
        {
           //try
           //{
               std::filesystem::path dir(_cp._OutputFile.get());
               dir.remove_filename();
               std::filesystem::create_directories(dir);

               return CProject::Run(prog);
        //   }
        //catch(std::exception& e)
        //{
             //(nana::msgbox(*_EdWd, ("std::exception during EditableWidget ReCollocation: "))
             //       .icon(nana::msgbox::icon_error)
             //                    <<("\n   in form: ") << nana::API::window_caption(*_EdWd)
             //                    <<("\n   exception : ") << e.what() 
             //).show();
  //      }
		//catch(...)
		//{
             //(nana::msgbox(*_EdWd, ("An uncaptured exception during EditableWidget ReCollocation: "))
             //       .icon(nana::msgbox::icon_error)
             //                    <<("\n   in form: ") << nana::API::window_caption(*_EdWd)
             //).show();
	    //}
     //      return -1;
        }       //  ???????


 explicit	ThDyProject():	CProject("ThDy DNA Hybrid Project.","Def.ThDy.txt","Def.ThDy.txt")
					{}
 	virtual void saveTMP() const     override       // Se me habia olvidado redefinir esta funcion para usar _cp.OutF
	{	  saveToFile(   ( _cp._OutputFile.get() + ".ThDy.txt" ).c_str()   );	}
};  




#endif


