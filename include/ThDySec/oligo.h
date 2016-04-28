/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file  ThDySec\include\ThDySec\oligo.h
*
* @brief ?? not in use
*/
#pragma unmanaged	
#ifndef _OLIGO_H
#define _OLIGO_H
#include "sec.h"

class Oligo
{
	LonSecPos _i;
	LonSecPos _f;		// or _l  unsigned char  ?? unsigned short
public:
	CSecCand &_SecCand;
	Oligo (CSecCand &sec,	LonSecPos i, LonSecPos f):_SecCand(sec),_i(i),		_f(f)		{}
	Oligo (CSecCand &sec,	const CRangBase& r)		 :_SecCand(sec),_i(r.Min()),_f(r.Max())	{}
												//char *begin(){return _pSec->}
	bool operator	==		(const Oligo &o	)const
					{	if(_f-_i != o._f-o._i) 
							return false;
						for(LonSecPos p1=_i, p2=o._i ; p1<=_f ; ++p1, ++p2)
							if (_SecCand._Sec[p1]!=o._SecCand._Sec[p2])
								return false;
						return true;
					}
	float			G		(				)const 
					{
						return _SecCand._Sec.G(_i,_f);
					}
	size_t			HashG	(				) const 
					{	double y;
							return size_t(  modf((double(G())-40)/(20-(-40)), &y) * std::numeric_limits<size_t>::max() ); 
					}
	static size_t	HashG	(Energy		g	)  
					{	double y;
							return size_t(  modf((double(g)-40)/(20-(-40)), &y) * std::numeric_limits<size_t>::max() ); 
					}
	virtual float	Gcroos	(				)const 
					{
						return G();
					}
	virtual float	Gcroos	(Temperature Ta	)const 
					{
						return _SecCand._Sec.G(Ta);
					}
};
									// http://en.cppreference.com/w/cpp/utility/hash
namespace std 
{		template<>
	class hash<Oligo> 
	{ public:
		size_t operator()(const Oligo &o) const 
		{	
			return o.HashG(); 
		}
	};
}
class CrossHy : public Oligo
{   Energy			_H, _G;
	Entropy			_S;
	Temperature		_Tm;
	float			_I;

public:
	CrossHy (Oligo &o)	:	Oligo (o)	{}

	CrossHy (CSecCand &sec,	LonSecPos i, LonSecPos f, Energy H, Entropy S, Energy G, Temperature Tm)
			:	Oligo (sec, i,  f),
				_H(H), _S(S), _G(G), _Tm(Tm) {} 
	CrossHy (CSecCand &sec,	const CRangBase& r,		  Energy H, Entropy S)	:	Oligo (sec,	 r)	,	_H(H), _S(S)	{}
	static Temperature	Ta;
	static Energy		Gsens, Gsat;
	static float		Isens, Isat;				
};

//#include <list>
#include <unordered_map>
class CGlobalOligoList
{
	unordered_map     <Oligo, unordered_multimap<CSecCand , CrossHy>> _oli;
	//unordered_multimap<Oligo, unordered_multimap<CSecCand , CrossHy>> _oli;
public:
	CGlobalOligoList( Temperature	Ta, Energy	Gsens, Energy Gsat, float Isens, float Isat)
	{
		CrossHy::Ta		=Ta;
		CrossHy::Gsens	=Gsens;
		CrossHy::Gsat	=Gsat;
		CrossHy::Isens	=Isens;
		CrossHy::Isat	=Isat;
	}
	void AddOligo(Oligo &o)
	{
		_oli[o].insert(std::pair<CSecCand , CrossHy>(o._SecCand,CrossHy(o))); ;
	}
	void AddCrossHy(Oligo &o, CrossHy& h)
	{
		_oli[o].insert(std::pair<CSecCand , CrossHy>(o._SecCand,h)); ;
	}

};


			//size_t h1 = std::hash<std::string>()(s.first_name);
			//size_t h2 = std::hash<std::string>()(s.last_name);
			// h1 ^ ( h2 << 1 );
//					 h								s	
//"	ndH(C,G,G,C)="	-10.0  	6f*1000; ndS(C,G,G,C)=	-27.2  	f;   // CG/GC 04
//"	ndH(G,C,C,G)="	-9.8  	f*1000;  ndS(G,C,C,G)=	-24.4  	f;   // GC/CG 04
//"	ndH(C,A,G,T)="	-8.5  	f*1000;  ndS(C,A,G,T)=	-22.7  	f;   // CA/GT 04
//"	ndH(T,G,A,C)="	-8.5  	f*1000;  ndS(T,G,A,C)=	-22.7  	f;   // TG/AC adapted CA/GT
//"	ndH(A,C,T,G)="	-8.4  	f*1000;  ndS(A,C,T,G)=	-22.4  	f;   // AC/TG adapted GT/CA
//"	ndH(G,T,C,A)="	-8.4  	f*1000;  ndS(G,T,C,A)=	-22.4  	f;   // GT/CA 04
//"	ndH(G,A,C,T)="	-8.2  	f*1000;  ndS(G,A,C,T)=	-22.2  	f;   // GA/CT 04
//"	ndH(T,C,A,G)="	-8.2  	f*1000;  ndS(T,C,A,G)=	-22.2  	f;   // TC/AG adapted GA/CT
//"	ndH(C,C,G,G)="	-8.0  	f*1000;  ndS(C,C,G,G)=	-19.9  	f;   // CC/GG adapted GG/CC
//"	ndH(G,G,C,C)="	-8.0  	f*1000;  ndS(G,G,C,C)=	-19.9  	f;   // GG/CC 04
//"	ndH(A,G,T,C)="	-7.8  	f*1000;  ndS(A,G,T,C)=	-21.0  	f;   // AG/TC adapted CT/GA
//"	ndH(C,T,G,A)="	-7.8  	f*1000;  ndS(C,T,G,A)=	-21.0  	f;   // CT/GA 04
//"	ndH(A,A,T,T)="	-7.6  	f*1000;  ndS(A,A,T,T)=	-21.3  	f;   // AA/TT 04
//"	ndH(T,T,A,A)="	-7.6  	f*1000;  ndS(T,T,A,A)=	-21.3  	f;
//"	ndH(A,T,T,A)="	-7.2  	f*1000;  ndS(A,T,T,A)=	-20.4  	"f;   // AT/TA 04 //	basek[]="".ACGT$"" :.-0, A-1, C-2, G-3, T-4, $-5 , g=0 /* gap */ , e=5 /* extremo, end */"
//"	ndH(T,A,A,T)="	-7.2  	f*1000;  ndS(T,A,A,T)=	-21.3  	f;   // TA/AT 04


#endif 