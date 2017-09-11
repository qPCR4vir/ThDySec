/**
* Copyright (C) 2009-2015, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
*
* @file  ThDySec\src\ThDy_programs\init_thdy_prog_param.cpp
*
* @brief
*
*/



//#include "StdAfx.h"
#pragma unmanaged
#include "thdy_programs\init_thdy_prog_param.h"
#include "ThDy_programs\prog_comm_functions.h"
#include "ThDySec/sec.h"
#include <assert.h>
using namespace std;


//ThDyCommProgParam::~ThDyCommProgParam(void)        {/*delete []_ProgList;*/}

CProgParam_microArray::~CProgParam_microArray()		{ /*delete _tlTm;*/}

std::shared_ptr<CMultSec>  ThDyCommProgParam::CreateRoot	()
{
	return std::make_shared<CMultSec> (_pSaltCorrNNp, "All seq") ;
}

CMultSec* ThDyCommProgParam::AddSeqGroup		(CMultSec   *parentGr, const std::string&     Name)
{
	if(parentGr)
	    if (  parentGr->_NNPar ) 
            return parentGr->AddMultiSec(new CMultSec(parentGr->_NNPar, Name));
		else
            return parentGr->AddMultiSec(new CMultSec(_pSaltCorrNNp, Name));

    return new CMultSec(_pSaltCorrNNp, Name);
}


CMultSec* ThDyCommProgParam::AddSeqFromFile( CMultSec           *parentGr, 
	                                         const std::string&  FileName, 
	                                         bool                recursive/*=false*/, 
	                                         bool                onlyStructure/*=false*/)
{
	if(parentGr)
    {
		CMultSec *sG{};

		for (CMultSec* pgr : _primersGr)
		{
			if (pgr == pgr->findComParent(parentGr) )   // adding in a primer group: dont use filtres, use parent filtres?
			{   
				  sG=new CMultSec ( FileName , 
                                    parentGr->_NNPar ? parentGr->_NNPar : _pSaltCorrNNp,
                                    recursive,
							        parentGr->_MaxTgId  ,
							        parentGr->_SecLim ,
                                    parentGr->_SecLenLim, 
                                   !onlyStructure);
				  sG->CreateNonDegSetRec();
				  break;
		    }
		}
	    if (!sG)
		          sG=new CMultSec ( FileName , 
                                    parentGr->_NNPar ? parentGr->_NNPar : _pSaltCorrNNp,
                                    recursive,
							        _MaxTgId  ,
							        _SecLim ,
                                    _SecLenLim, 
                                   !onlyStructure);

		parentGr->AddMultiSec(sG);
        std::string parent_path = filesystem::path(sG->_Path).remove_filename().string();

        if (parentGr->_Local._NMSec == 1)
            parentGr->_Path= parent_path ;
        else if (parentGr->_Path != parent_path)
                parentGr->_Path.clear();

		return sG;
    }
	else
	return new  CMultSec  (     FileName , 
                                _pSaltCorrNNp, 
                                recursive,
							   _MaxTgId,
							   _SecLim ,
                               _SecLenLim, 
                               !onlyStructure);

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
CProgParam_MultiplexPCR::CProgParam_MultiplexPCR(const string& titel, ThDyCommProgParam &commThDyParam) 
	: CProgParam_microArray(titel,commThDyParam), _rtbl_self(nullptr)
	{	
        _InputSondeFile.SetTitel("Imput file for primers"); 
		_InputSondeFile.SetEtiq("iSonde_PCR", this); 
        _probesMS->_name="Primers of Multiplex PCR";

        _PrRecurDir.SetTitel("Recursively add all primers seq-files from all dir"); 
		_PrRecurDir.SetEtiq("PrimRecDir", this); 

        _PrDirStrOnly.SetTitel("Reproduce only the dir struct in primers"); 
		_PrDirStrOnly.SetEtiq("PrimDirStr", this); 

 
	}


//#include "ThDySec\sec.h"

