/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\src\ThDySec\link.cpp
*
* @brief 
*/

#ifdef WINDOWS_FORM_GUI
#include "stdafx.h"
#pragma unmanaged
#endif

#include "ThDySec/link.h"



void CLink::Insert ( CLink *p, CLink *n ) // 		void Insert ( CLink *p, CLink *n=nullptr ) ; 
{	
	// CUIDADO con las ramificaciones, use Remove() or Move() para evitarlas ;
	if( p ) 
		if (p->_next == n ) // pone a this como next de p SOLO si ya antes n era el next de p !!!???
			p->_next = this ;

	if( n ) 
		if (n->_prev == p )
			n->_prev = this ;
	
	_next = n ;
	_prev = p ;
}

void CLink::Remove()
{	if( _next ) 
		if ( _next->_prev == this) 
			 _next->_prev = _prev ;

	if( _prev ) 
		if ( _prev->_next == this) 
			 _prev->_next = _next ;
					
	_prev = _next = nullptr ;
}

void CList::Insert (CLink *link)   // Asume que Cur es correcto !!!. insert 'link' antes de 'Cur()'. Si no existe Cur() lo anade al final, antes de last
{	if (!link) return ;			//CLink &c= *this->Cur();
	if (this->Cur())
			link->MoveBefore (this->Cur()) ;	//	link->Insert (this->Cur()->_prev, this->Cur()) ;   // insert 'link' antes de 'Cur'.
	else	Add (link)	;					//	Ahora si no existe 'Cur' Add el link AL FINAL de la CList   !!???
}

void CList::Destroy ()				// funciona bien solo si la lista es "lineal". Destruye todos los miembros de la lista. Usar cuando la lista es "duena" de ellos
{	while (_First.Next() != &_Last)
			delete _First.Next() ;	// si la lista es lineal como ultimo se borra "last._prev" tambien
}
void CList::free ()				// funciona bien solo si la lista es "lineal". Solo libera los miembros de la lista. Usar cuando la lista no es "duena" de ellos
{	while (_First.Next() != &_Last)
			 _First.Next()->Remove() ;	// si la lista es lineal como ultimo se borra "last._prev" tambien
}
//void CList::Add (CLink *link) // Anade al final, antes de last, que se queda siempre de last. independiente de CList::Insert (CLink *link), usa CLink::Insert ( CLink *p, CLink *n=nullptr  )
//{	if (!link) return ;
//	
//	link->InsertBefore (&_Last) ;
//
//	//if (!_first) 
//	//	_first	= link ;	
//
//	//_last = link ;
//}
//void CList::Insert (CLink *link)   // Asume que Cur es correcto !!!. insert 'link' antes de 'Cur()'. Si no existe Cur() lo anade al final, antes de last
//{	if (!link) return ;
//	//CLink &c= *this->Cur();
//	if (this->Cur())
//		link->InsertBefore (this->Cur()) ;	//	link->Insert (this->Cur()->_prev, this->Cur()) ;   // insert 'link' antes de 'Cur'.
//	else	Add (link)	;					//	Ahora si no existe 'Cur' Add el link AL FINAL de la CList   !!???
//
//	if (_first == this->Cur() || !_first)   // _first == this->Cur() :  por tanto this->Cur()->_prev seria =0. Innecesario-> || !_first
//		_first	= link ;					// y _first == this->Cur() ahora "queda" detras de link, es su next, y link queda de primero
//
//	if (!_last)
//		_last   = link ;
//}
//void CList::Destroy ()				// funciona bien solo si la lista es "lineal"
//{	
//	if(_first)
//		while (_first->_next)
//			delete _first->_next ;	// si la lista es lineal como ultimo se borra "last" tambien
//	//if(_last)
//	//	while (_last->_prev)
//	//		delete _last->_prev ;		
//		
//	//	for (_cur = _first->_next ; _cur ; _cur = _first->_next ) // va borrando siempre el segundo elemento de la lista, last incluso
//	//	{	delete _cur ; /*_first->_next=0;*/   }	// Como el destructor de CLink es virtual primero se invoca el destructor redefinido del tipo concreto
//													// y luego se invoca Remove() que quita el link de la lista
//	delete _first ;
////	delete _last  ;
//	_first=_last=_cur=nullptr;
//}