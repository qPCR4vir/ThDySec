/**
* Copyright (C) 2009-2019, Ariel Vina-Rodriguez ( arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
*
* @file  ThDySec\src\ThDy_programs\init_thdy_prog_param.cpp
*
* @brief
*
*/

#include <cassert>

#include "ThDy_programs/init_ThDy_prog_param.h"
#include "ThDy_programs/prog_comm_functions.h"
#include "ThDySec/sec.h"

namespace fs = std::filesystem;


//ThDyCommProgParam::~ThDyCommProgParam(void)        {/*delete []_ProgList;*/}

CProgParam_microArray::~CProgParam_microArray()		{ /*delete _tlTm;*/}

std::shared_ptr<CMultSec>  ThDyCommProgParam::CreateRoot	()
{
	return std::make_shared<CMultSec> (_pSaltCorrNNp, "All seq") ;
}

CMultSec::MSecIt ThDyCommProgParam::AddSeqGroup(CMultSec &parentGr, const std::string& Name)
{
	return parentGr.AddMultiSec(
		std::make_shared<CMultSec>(parentGr._NNPar ? parentGr._NNPar : _pSaltCorrNNp, Name)   );
}


CMultSec::MSecIt ThDyCommProgParam::AddSeqFromFile(CMultSec &parentGr,
                                                   const std::filesystem::path &FileName,
                                                   bool recursive/*=false*/,
                                                   bool onlyStructure/*=false*/)
{
    CMultSec::pMSec sG;
    // Primers and normal sequences are added differently: primers are not filtered.
    // We need to know if we are adding primers. For that we will test if one of the parents of the given
    // parent is in the list (set) of primer trees.
	for (CMultSec* pgr : _primersGr)       // scan primers groups
	{
		if (pgr == pgr->findComParent(&parentGr))   // adding in a primer group: dont use filters
		{   
				sG=std::make_shared<CMultSec>(FileName ,
                                              parentGr._NNPar ? parentGr._NNPar : _pSaltCorrNNp,
                                              recursive,
							                  parentGr._MaxTgId  ,   // use parent filters?
							                  parentGr._SecLim ,
                                              parentGr._SecLenLim,
                                              !onlyStructure);
				sG->CreateNonDegSetRec();
				break;
		}
	}
	if (!sG)    // not primers, just some target sequences. Will be filtered.
	    sG=std::make_shared<CMultSec>(FileName ,
                                      parentGr._NNPar ? parentGr._NNPar : _pSaltCorrNNp,
                                      recursive,
							          _MaxTgId  ,
							          _SecLim ,
                                      _SecLenLim,
                                      !onlyStructure);
    CMultSec::MSecIt it = parentGr.AddMultiSec(sG);
    auto parent_path = sG->_orig_file.parent_path();

    if (parentGr._Local._NMSec == 1)
        parentGr._orig_file = parent_path ;
    else if (parentGr._orig_file != parent_path)
            parentGr._orig_file.clear();
	return it;
}

void CProgParam_microArray::RenameSondesMS(const std::string& name)
{
    _probesMS->_name=name;
}


//CProgParam_uArrExp::CProgParam_uArrExp (const string& titel, ThDyCommProgParam &commThDyParam):	
//	                _exclSd(false),		
//					_IxI(true),		_IxI_d(true),			
//					_Normalize(true), 
//					_Isat(Energy(0.87f)),_Isen(Energy(0.01f)),	_Gsat(Energy(-2.0f)),	_Gsen(Energy(2.0f)),
//					CProgParam_microArray (titel,commThDyParam) 
//					{  _probesMS->_name="Probes of Exp uArr";
//					} 
CProgParam_MultiplexPCR::CProgParam_MultiplexPCR(const std::string& titel, ThDyCommProgParam &commThDyParam)
	: CProgParam_microArray(titel,commThDyParam), _rtbl_self(nullptr)
	{	
        _InputSondeFile.SetTitel("Imput file for primers"); 
		_InputSondeFile.SetEtiq("iSonde_PCR", this); 
        _probesMS->_name="Primers for Multiplex PCR";

        _PrRecurDir.SetTitel("Recursively add all primers seq-files from all dir"); 
		_PrRecurDir.SetEtiq("PrimRecDir", this); 

        _PrDirStrOnly.SetTitel("Reproduce only the dir struct in primers"); 
		_PrDirStrOnly.SetEtiq("PrimDirStr", this); 

 
	}


//#include "ThDySec\sec.h"

