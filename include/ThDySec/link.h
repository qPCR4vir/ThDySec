/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\include\ThDySec\th_dy_param.h
*
* @brief
*/

#pragma unmanaged	
#ifndef _LINK_H
#define _LINK_H
#include <cassert>
// #include "cod_deg.h"



class CList  ;
class CLink     // base para los clases que quieran ser elementos de listas, tal vez demaciado listo
				// si no se redefine adecuadamente el destructor del link y del dueno de la lista entonces la lista no se "aduena" de sus elementos
{	// friend CList ;
	public:
		CLink (CLink *p, CLink *n=nullptr) : _next (n), _prev(p) {Insert (p,n);} ;
		CLink () : _next (nullptr), _prev(nullptr) {} ;


		void Remove ();								// y si pertenecia a una Clist y coincidia con  cur, last or first ??

		void Move ( CLink *p, CLink *n=nullptr ) { if(this==p || this==n) return; Remove ();    Insert (p,n);};
		void MoveAfter (CLink *p)				 {assert(p);	Move(p,p->_next)  ;}
		void MoveBefore(CLink *n)				 {assert(n);	Move(n->_prev,n)  ;}

		CLink *Next(){return _next;}
		CLink *Prev(){return _prev;}

		virtual ~CLink() {if (this) Remove();}	// Cada derivado debe decidir si implementar destructor "adicional"
												// Este solo desconecta el link de donde estaba, "cerrando" el hueco que dejaria
												// y dejando el link "huerfano"
		void InsertBefore(CLink *n)				{assert(n);	  Insert(n->_prev,n)  ;}		// Poner PRIVATE ??????????				
		void InsertAfter (CLink *p)				{assert(p);	  Insert(p,p->_next)  ;}
	private:
		CLink *_next, *_prev ;
		void Insert ( CLink *p, CLink *n=nullptr ) ;  // CUIDADO con las ramificaciones !!!

};

class CList    // No tiene destructor. No se aduena de los elementos, pero el dueno puedo usar  Destroy () para delete'rlos
{	public:
		CList () : _First(nullptr, &_Last), _Last( &_First, nullptr), _cur(_First.Next()) {} 
		void Add	(CLink *link) {	if (link) link->MoveBefore   (&_Last) ;}    // adiciona al final, "detras" de _last
		void Acople	(CLink *link) {	if (link) link->InsertBefore (&_Last) ;}    // adiciona al final, "detras" de _last
		void Insert (CLink *link) ;   // insert 'link' antes de 'this'

		void Destroy () ;				// "desbarata" la lista, y si el link tiene redefinido el destructor tambien los borra
		void free() ;
		void SetCur(CLink *cur){_cur =cur;} // CUIDADO, no se "salga" de la lista !!!!!
		CLink *goBeging()	{return _cur=_First.Next(); }
		CLink *goLast()		{return _cur=_Last.Prev(); }
		CLink *goNext()		{return _cur= _cur ? (_cur->Next()) : nullptr ; }  // devolver Last o nullptr ???
		CLink *goPrev()		{return _cur= _cur ? (_cur->Prev()) : nullptr ; }
		CLink *Cur()		{return _cur; }
		CLink *Last()		{return _Last.Prev(); }
		CLink *First()		{return _First.Next(); }

		bool NotEnd()		{return _cur!= &_First && _cur!= &_Last && _cur!= nullptr; }		/*   || _cur== _last     ;*/
		bool    End()		{return _cur== &_First || _cur== &_Last || _cur== nullptr; }		/*   || _cur== _last     ;*/
	private:
		 CLink _First, _Last, *_cur ;// first, last and current --- solo un "cur" cada vez en contraposicion con un "iterator" independiente
};

//class CList    // No tiene destructor. No se aduena de los elementos, pero el dueno puedo usar  Destroy () para delete'rlos
//{	public:
//		CList () : _first(nullptr), _last(nullptr), _cur(nullptr) {} 
//		void Add	(CLink *link) ;   // adiciona al final, "detras" de _last
//		void Insert (CLink *link) ;   // insert 'link' antes de 'this'
//
//		void Destroy () ;				// "desbarata" la lista, y si el link tiene redefinido el destructor tambien los borra
//
//		void SetCur(CLink *cur){_cur =cur;} // CUIDADO, no se "salga" de la lista !!!!!
//		CLink *goBeging()	{return _cur=_first; }
//		CLink *goLast()		{return _cur=_last; }
//		CLink *goNext()		{return _cur= _cur ? _cur->Next() : nullptr ; }
//		CLink *goPrev()		{return _cur= _cur ? _cur->Prev() : nullptr ; }
//		CLink *Cur()		{return _cur; }
//		CLink *Last()		{return _last; }
//		CLink *First()		{return _first; }
//
//		bool NotEnd()		{return _cur!= nullptr ; }		/*   || _cur== _last     ;*/
//	private:
//		 CLink *_first, *_last, *_cur ;// first, last and current --- solo un "cur" cada vez en contraposicion con un "iterator" independiente
//};

#endif

