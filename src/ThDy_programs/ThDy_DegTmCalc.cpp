/**
* Copyright (C) 2009-2015, Ariel Vina Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2015
*
* @file ThDySec\src\ThDy_programs\ThDy_DegTmCalc.cpp
*
* @brief
*/

//#include "StdAfx.h"
#pragma unmanaged
#include "ThDy_programs/prog_comm_functions.h"


int MultiplexPCRProg ( CProgParam_MultiplexPCR *IPrgPar_uArr, 	CMultSec		&primers	)  ;

int DegTmCalc ( CProgParam_TmCalc *IPrgPar_Calc)  
{
	const int MaxGrDeg=300 ;			// crear NonDegSet para las sondas con menos de este gr de deg. Poner como ProgParam??

    if (IPrgPar_Calc->save.get())  
        IPrgPar_Calc->_cp.Check_NNp_Targets ();
    else 
        IPrgPar_Calc->_cp.Actualice_NNp ( );

	Temperature Ta=  IPrgPar_Calc->_cp._pSaltCorrNNp->Ta() ; 

	auto Sec = std::make_shared<CSec>(	IPrgPar_Calc->_Sec.get(), 0, "Sec", IPrgPar_Calc->_cp._pSaltCorrNNp);

	if (Sec->Len() < 1)  return 0 ; // Error :  no sec !!!!!!
	Sec->CreateNonDegSet();	

	if (CountDegBases(			IPrgPar_Calc->_Sec2Align.get().c_str())		< 1)				
								IPrgPar_Calc->Update_Sec_Sec2Align(true,true);	

	auto Sec2Align = std::make_shared<CSec>(IPrgPar_Calc->_Sec2Align.get(), 0, "Sec2Align",	IPrgPar_Calc->_cp._pSaltCorrNNp);
	Sec2Align->CreateNonDegSet();

	std::shared_ptr<CMultSec>	pr,	tg;			// Esto se puede hacer mejor

	if   (      Sec->NonDegSet()) pr=      Sec->NonDegSet() ; else {pr.reset(new CMultSec(IPrgPar_Calc->_cp._pSaltCorrNNp)); pr->AddSec(      Sec);}
	if   (Sec2Align->NonDegSet()) tg=Sec2Align->NonDegSet() ; else {tg.reset(new CMultSec(IPrgPar_Calc->_cp._pSaltCorrNNp)); tg->AddSec(Sec2Align);}

		
	IPrgPar_Calc->_TmS  = KtoC(pr->_Local._Tm) ;// (    KtoC(pr->_minTm)   ,   KtoC(pr->_maxTm)   ) ; 
	IPrgPar_Calc->_Tm2A = KtoC(tg->_Local._Tm) ;

	CSec *pr_maxTmH= pr->SecL().begin()->get()  ;	IPrgPar_Calc->_GS.Set ( pr_maxTmH->G()/1000  )  ; 
	CSec *tg_maxTmH= tg->SecL().begin()->get()  ;   IPrgPar_Calc->_G2A.Set( tg_maxTmH->G()/1000  )  ;
														
	std::unique_ptr<ThDyAlign> apAl; 

	//LonSecPos TgMaxLen= (tg->_TMaxLen > pr->_TMaxLen) ? tg->_TMaxLen : pr->_TMaxLen ;

	if ( IPrgPar_Calc->align.get())	
	{	apAl= Create_ThDyAlign(	IPrgPar_Calc->_cp, pr->_Global._Len.Max(), 
		                                           tg->_Global._Len.Max(), IPrgPar_Calc->_cp._pSaltCorrNNp);
		
	    apAl->Align( pr_maxTmH, tg_maxTmH);					
	    apAl->SelectOptParam( Ta);	//  virtual !!! Si G la Ta pudo cambiar, por eso aqui explicita
								
		IPrgPar_Calc->_TmHy.Set ( KtoC( apAl->Tm() ) );		//	FrAl.GetOptHit();
		IPrgPar_Calc->_GHy.Set  ( apAl->G ()/1000    );		//		print_ThDyAlign (osAl, Al);	//Al.Export_DPMz_Pre(osAl);

	}
	else 
	{	AlignedSecPar al( pr_maxTmH->Sequence(), 
		                  tg_maxTmH->Sequence(), IPrgPar_Calc->_cp._pSaltCorrNNp ); // la Ta en NNpar no cambio

		IPrgPar_Calc->_TmHy.Set ( KtoC( al.Tm() ) );
		IPrgPar_Calc->_GHy.Set  (	al.G ()/1000  );
	}

	for ( auto &CurSec : 	pr->SecL() )     // recorre todos las var no deg de la sonda
	{	CSec &s = *CurSec ; 					 Energy  g= s.G (Ta)/1000;			Temperature tm ;
		IPrgPar_Calc->_GS.Expand(g);
		
		for (auto &tg_CurSec : tg->SecL())   // recorre todos las var no deg de la sonda
		{	CSec &t = *tg_CurSec ;					     g= t.G (Ta)/1000 ;
			IPrgPar_Calc->_G2A.Expand(g) ;

			if ( IPrgPar_Calc->align.get())	
			{	ThDyAlign	&Al=*apAl.get();
				Al.Align( &(s), &(t));				Al.SelectOptParam(Ta);			//	FrAl.GetOptHit();					
															 g= Al.G ()/1000 ;			tm=  KtoC( Al.Tm() ) ;
				IPrgPar_Calc->_GHy.Expand(g) ;
				if (IPrgPar_Calc->_TmHy.Max() <=  tm  ) 
				{	IPrgPar_Calc->_TmHy.Max()  =  tm; 
					CHitAligned Hit (Al);
					IPrgPar_Calc->Set_AlignedSec      ( (char*)(Hit._sd.c_str() ) /*)*/);
					IPrgPar_Calc->Set_AlignedSec2Align( (char*)(Hit._tg.c_str() ) /*)*/);
				} 
				else if  (IPrgPar_Calc->_TmHy.Min() >  tm  ) {IPrgPar_Calc->_TmHy.Min() = tm; }

			} else 
			{
				AlignedSecPar al( (s.Sequence())  , (t.Sequence()), IPrgPar_Calc->_cp._pSaltCorrNNp ); 		
															g= al.G ()/1000 ;	float tm=KtoC( al.Tm() ) ; 
				IPrgPar_Calc->_GHy.Expand(g) ;
				if       (IPrgPar_Calc->_TmHy.Max() <  tm  ) 
				{	IPrgPar_Calc->_TmHy.Max() = tm; 
					pr_maxTmH=&s; 
					tg_maxTmH=&t;
				} 
				else if  (IPrgPar_Calc->_TmHy.Min() >  tm  ) {IPrgPar_Calc->_TmHy.Min() = tm; }
			}

		}// recorre todos las var no deg
	}
	if ( ! IPrgPar_Calc->align.get())	
	{
		IPrgPar_Calc->Set_AlignedSec      ( pr_maxTmH->charSequence());
		IPrgPar_Calc->Set_AlignedSec2Align( (char*)tg_maxTmH->Sequence().c_str()  );
	}

	//delete pAl;	
	if (IPrgPar_Calc->save.get())	
	{
		CMultSec primers(IPrgPar_Calc->_cp._pSaltCorrNNp); 
		primers.AddMultiSec( pr);
		primers.AddMultiSec( tg);

		int t=MultiplexPCRProg ( IPrgPar_Calc, primers		)  ;
		return t;
	}

	return 1;

}

	//IPrgPar_uArr->_GS.Set();
	//IPrgPar_uArr->_G2A.Set();
	//IPrgPar_uArr->_GHy.Set();
	//IPrgPar_uArr->_TmHy.Set();
	//CSaltCorrNN		_NNpar	(			IPrgPar_uArr->_cp._ConcSd, 
	//									IPrgPar_uArr->_cp._ConcTg, 
	//									IPrgPar_uArr->_cp._ConcSalt   ); 		CSaltCorrNN *NNpar = &_NNpar;
	//_NNpar.SetTa(				CtoK(	IPrgPar_uArr->_cp._Ta));			// Aqui por si acaso. Revisar.
	//if (IPrgPar_uArr->_cp._loadNNPar){ifstream isTDP(IPrgPar_uArr->_cp._InputNNFile.Get());	assert(isTDP);	NNpar->LoadNNParam(isTDP) ;	}
	//if (IPrgPar_uArr->_cp._saveNNPar)
	//{	string OutputTDP(IPrgPar_uArr->_cp._OutputFile.Get()) ; 		OutputTDP += ".ThDyParam.csv";
	//	ofstream osTDP	(OutputTDP.c_str());				assert(osTDP);	
	//	osTDP << _NNpar ;
	//}

		//switch (	IPrgPar_uArr->_cp._TAMeth )
		//{	case TAMeth_Tm: default:	pAl=new ThDyAlign_Tm( pr->_TMaxLen ,  TgMaxLen, *NNpar);   break;
		//	case TAMeth_Fract:			pAl=new FracTDAlign ( pr->_TMaxLen ,  TgMaxLen, *NNpar);   break;
		//	case TAMeth_G:				pAl=new ThDyAlign_G ( pr->_TMaxLen ,  TgMaxLen, *NNpar);   break;
		//}
		//ThDyAlign	&Al=*pAl;
		//Al.SetTa		 (			CtoK(	IPrgPar_uArr->_cp._Ta	    ));		// OK

//AlignedSecPar al( reinterpret_cast <char *> (pr_maxTmH->_c)  , reinterpret_cast <char *> (tg_maxTmH->_c), *NNpar );