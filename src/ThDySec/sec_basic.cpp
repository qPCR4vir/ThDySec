/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\src\ThDySec\sec_basic.cpp
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

using namespace std ; 

#include "ThDySec/sec_mult.h"
#include "ThDySec/common.h" 
using namespace DegCod;

//#define SEQUENCES_MAX_SIZE 100000
char *DNAstrandName[]=	{""		, "(c)", ""		, "(r)"	, "(i)", "(c)"		} ;
// enum DNAstrand		{plus	, minus, direct	, rev	, compl, rev_compl	} ;


CSecBasInfo::~CSecBasInfo()
{
	if (_NonDegSet) 
	     if (!_NonDegSet->Prev() && !_NonDegSet->Next()) 
			 delete _NonDegSet;
	// en otro caso, donde borrar _NonDegSet ????. Lo borra la lista en la que esta insertado
}
	 
std::string& CSecBasInfo::Copy_Seq  	(std::string &SecHier,  long InicBase, long EndBase, DNAstrand strnd) const
{	
	SecHier.clear();
	if ( EndBase< 1 || Len() <EndBase ) EndBase= Len(); 
	long l=EndBase-InicBase+1 ;  
	if (l>=0) 	//assert(l>=0);
	{
		SecHier.reserve(l); /// \todo resize ?? and use [] directly instead of push_back 
		switch (strnd)
		{	case DNAstrand::plus :
			case DNAstrand::direct:		for(long p=InicBase;  p<=EndBase;    p++) 	SecHier.push_back ( _c[p] );  				    break;
			case DNAstrand::compl:		for(long p=InicBase;  p<=EndBase;    p++) 	SecHier.push_back ( c_degbase[_c[p]]);      	break;
			case DNAstrand::rev:		for(long p=EndBase ;  p>=InicBase;   p--)	SecHier.push_back ( _c[p] );  				    break;
			case DNAstrand::minus:
			case DNAstrand::rev_compl:	for(long p=EndBase ;  p>=InicBase;   p--)	SecHier.push_back ( c_degbase[_c[p]]);      	break;
 		}
	}
	return  SecHier ;
}

/// Best just return sequence
//Base  *	CSecBasInfo::Copy_charSec(Base *charSecHier,long InicBase, long EndBase, DNAstrand strnd)//DNAstrand strnd=direct)
//{	if ( EndBase< 1 || Len() <EndBase ) EndBase= Len(); 
//	long l=EndBase-InicBase+1 ; charSecHier[l]=0 ;
//	if (l>=0) 
//	//assert(l>=0);
//	switch (strnd)
//	{	case DNAstrand::plus :
//		case DNAstrand::direct:		for(long i=0,   p=InicBase;  p<=EndBase;    i++, p++) 	charSecHier[i]=_c[p];				break;
//		case DNAstrand::compl:		for(long i=0,   p=InicBase;  p<=EndBase;    i++, p++) 	charSecHier[i]=c_degbase[_c[p]];	break;
//		case DNAstrand::rev:		for(long i=l-1, p=InicBase;  p<=EndBase;    i--, p++)	charSecHier[i]=_c[p];				break;
//		case DNAstrand::minus:
//		case DNAstrand::rev_compl:	for(long i=l-1, p=InicBase;  p<=EndBase;    i--, p++)	charSecHier[i]=c_degbase[_c[p]];	break;
//
//		default : return 0;
//	}
//	return charSecHier ;
//}
//Base  *	CSecBasInfo::GetCopy_charSec(DNAstrand strnd)
//{	return GetCopy_charSec(1, Len(), strnd);
//}
//Base  *	CSecBasInfo::GetCopy_charSec(long InicBase, long EndBase, DNAstrand strnd)       // recuerde los $...$, aqui se cuentan, 
//{	if ( InicBase< 1 )					InicBase=1;
//	if ( EndBase< 1 || Len() <EndBase )	EndBase=Len();
//	long l=EndBase-InicBase+1 ;
//	assert(l>=0);
//	Base *charSecHier=new Base[l+1];						 // asi como InicBase y EndBase inclusive!!
//	assert(charSecHier);	
//
//	return Copy_charSec(charSecHier, InicBase, EndBase, strnd)  ;
//}



