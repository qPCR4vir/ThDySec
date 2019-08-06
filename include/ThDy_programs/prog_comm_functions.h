/**
* Copyright (C) 2009-2016, Ariel Vina-Rodriguez ( ariel.rodriguez@fli.bund.de , arielvina@yahoo.es )
*  https://www.fli.de/en/institutes/institut-fuer-neue-und-neuartige-tierseuchenerreger/wissenschaftlerinnen/prof-dr-m-h-groschup/
*  distributed under the GNU General Public License, see <http://www.gnu.org/licenses/>.
*
* @autor Ariel Vina-Rodriguez (qPCR4vir)
* 2012-2016
*
* @file  ThDySec\include\ThDy_programs\prog_comm_functions.h
*
* @brief 
*
*/

#ifndef _PROG_COMM_FUNTIONS_H
#define _PROG_COMM_FUNTIONS_H
#pragma unmanaged
#include <memory>

#include "init_ThDy_prog_param.h"
#include "..\ThDySec\th_dy_align.h"
#include "..\..\ProgParam\include\matrix.h"


//unique_ptr<CSaltCorrNN> Create_NNpar    (ThDyCommProgParam& _cp);
//
//void Check_NNp_Targets (ThDyCommProgParam& cp);


std::unique_ptr<ThDyAlign>   Create_ThDyAlign(const ThDyCommProgParam& _cp, LonSecPos MaxLenSond, LonSecPos MaxLenTarg, std::shared_ptr<CSaltCorrNN>  NNpar);

class OutStr
{public:
	std::ofstream &Tm, &G, &Pos, &Pl_Tm, &Pl_G, &Al;
	OutStr( std::ofstream &osTm,
			std::ofstream &osG,
			std::ofstream &osPos,
			std::ofstream &osPl_Tm,
			std::ofstream &osPl_G,
			std::ofstream &osAl
		  )
			:Tm(osTm), G(osG),Pos(osPos),Pl_Tm(osPl_Tm),Pl_G(osPl_G),Al(osAl)
			{}
};


inline void Hybrid(CSec &s, CSec &t, 	ThDyAlign &Al,  std::ofstream &osTm,
														std::ofstream &osG,
														std::ofstream &osPos,
														std::ofstream &osPl_Tm,
														std::ofstream &osPl_G,
														std::ofstream &osAl,
														CTable<TmGPos, std::string> *_rtbl = nullptr	/*,
														CTable<Temperature> *_tlTm = nullptr	,
														CTable<Energy>	*tlG = nullptr,
														CTable<SecPos> *tlPos  = nullptr*/);


void HybridPr(CMultSec &pr, CSec &t, 	ThDyAlign &Al,  std::ofstream &osTm,
														std::ofstream &osG,
														std::ofstream &osPos,
														std::ofstream &osPl_Tm,
														std::ofstream &osPl_G,
														std::ofstream &osAl,
														CTable<TmGPos, std::string> *_rtbl = nullptr	/*,
														CTable<Temperature> *_tlTm = nullptr	,
														CTable<Energy>	*tlG = nullptr,
														CTable<SecPos> *tlPos  = nullptr*/);


#endif