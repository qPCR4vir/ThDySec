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

// enum DNAstrand		{plus	, minus, direct	, rev	, compl, rev_compl	} ;

CSecBasInfo::~CSecBasInfo()
{
	//if (_NonDegSet && !_NonDegSet->_parentMS)
	     //if (!_NonDegSet->Prev() && !_NonDegSet->Next()) 
	//		 delete _NonDegSet;
	/// \todo en otro caso, donde borrar _NonDegSet ????. 
	// Lo borra la lista en la que esta insertado
	// Pero como se que la otra lista lo borro?
	// poner todo como shared_ptr ???
	// por ahora entregar la "propiedad"
}

std::string CSecBasInfo::Copy_Seq  	(long InicBase, 
	                                 long EndBase, 
	                                 DNAstrand strnd) const
{	
	std::string SecHier;
	if ( EndBase< 1 || Len() <EndBase ) EndBase= Len(); 
	long l=EndBase-InicBase+1 ;  
	if (l>=0) 	//assert(l>=0);
	{
		SecHier.reserve(l); /// \todo resize ?? and use [] directly instead of += 
		switch (strnd)
		{	case DNAstrand::plus :
			case DNAstrand::direct:		for(long p=InicBase;  p<=EndBase;    p++) 	SecHier += _c[p] ;  		    break;
			case DNAstrand::complem:	for(long p=InicBase;  p<=EndBase;    p++) 	SecHier += c_degbase[_c[p]];   	break;
			case DNAstrand::rev:		for(long p=EndBase ;  p>=InicBase;   p--)	SecHier += _c[p] ;  		    break;
			case DNAstrand::minus:
			case DNAstrand::rev_compl:	for(long p=EndBase ;  p>=InicBase;   p--)	SecHier += c_degbase[_c[p]];   	break;
 		}
	}
	return SecHier;
}




