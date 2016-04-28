#ifndef _INIT_PROG_PARAM_IMPL_H
#define _INIT_PROG_PARAM_IMPL_H
#define _CRT_SECURE_NO_WARNINGS
#pragma unmanaged
#pragma warning( disable : 4996 )
#include <iostream>
#include <string>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <map>
#include <vector>
#include <functional>

#include "..\ThDySec\common.h" 

/// \brief Organize a "software" or Project into Specialized Programns and manage the input/config parametrs for each of the programs. 
/// Make ease to write a program interfase with a command-line, a text "project file" or a GUI, or all of then together. 
/// Definitions and declarations for the user interface. Also to be used directly by the programs. Almost primary, depend only on Common basics.
///  \todo:  PROBLEMA : como organizar estos parametros si usamos procesos? Hacer copia de ellos !!!!!!!!?
namespace Programs{
class IProg ;
/// Non abstract basic Interface of programms parametr which define default save, load and title use.
class IBParam
{
   std::string _Titel;
 public:

   IBParam (  const std::string &titel)      //  pp->?????????????
		    :_Titel(titel) 
	        {         
	        }    
	 virtual ~IBParam(){}

	 std::string Titel   ()const{return _Titel;}               ///<  ???  Human redeable
	 void        SetTitel(std::string titel){ _Titel=titel;}   ///<  ???  Human redeable


 virtual std::ostream &save	(std::ostream	&osPr) const        ///< The default. To be change in derivate classes
			            {   osPr<< ".\t"<<Titel()<<std::endl; 
							return osPr;
			            } 
 virtual bool   load	(std::istream   &isPr)  ///< throw( std::out_of_range)                
			            {   return false;}   
 virtual bool   load	(std::string &etiq, std::istream &isPr) ///<  throw( std::out_of_range) 
			            {   return false;}   

};
   // ifstream& operator >>(ifstream& ifs,IParam& p);                             //    ?????????????????
   // ofstream& operator <<(ofstream& ofs,const IParam& p){ p.save(ofs); return ofs;};  //    ?????????????????

///  \brief Base clase to manage a parametr of "some" type (defined in derived classes).  
///
/// It will call a IParam::ValueChanged callback function and set IParam::Unit and IParam::Etiq. 
/// Used to save to and load from an stream, without implementing the load and IParam::saveValue functions. 
class IParam : public IBParam
{    
    std::string _etiq, _unit;
 protected: 
	void changed()
	        {   
                if(ValueChanged) 
				    ValueChanged(/**this*/);
	        }

 public:
     std::function<void(void/*IBParam& param*/)> ValueChanged ; ///< A callback each time the value change

	 IParam (  IProg *pp,                     ///< the programm who own the parametr (tipicaly "this", the object for which the parametr is a member) 
               const std::string& titel, 
               const std::string& etiquete, 
               const std::string& unit="" ) ;

	 std::string Etiq(           )const{return _etiq;}      ///< semiHuman readable and unique. I use with length 10
	 void     SetEtiq(std::string etiq, IProg *prog);
	
	 std::string Unit()const{return _unit;}                 ///< Human redeable and optional

	 std::ostream	&save	(std::ostream	&osPr) const override  ///< Implemented using IParam::saveValuee plus a call to IBParam::save 
			            {    
                            osPr<< _etiq << ":\t"; 
			                saveValue(osPr)<<"\t"<<_unit<<"\t";
							IBParam::save(osPr)	;              ///< << ".\t"<<Titel()<<std::endl; 
							return osPr;
			            } 
     bool       load	(std::istream   &isPr)    override   ///< Asume etiquete allready tested OK !!!!   throw( std::out_of_range).  
			            {   
                            return loadValue(isPr);}   
	 bool       load	(std::string		&etiq, std::istream &isPr)  override  ///<   throw( std::out_of_range).  
			            {   
                            if (etiq!=_etiq) 
						        return false;
			                return load(isPr);
			            }  
                                            /// Default behavior, not yet a real implementation 
	virtual std::ostream	    &saveValue	(std::ostream	&osPr	) const  ///<  =0;   ??No salva nada, no tiene "value" todavia
	                                {return osPr;} 
                                            /// Default behavior, not yet a real implementation 
    virtual bool        loadValue	(std::istream   &isPr)  ///<   throw( std::out_of_range).  
	                                {return false;}         // =0;   ??    ?????????????????

	  ~IParam()override{}
};

           /// \brief Only partialy manage a parametr of type Num (a "numeric" type) for with the value have to be in a range defined by min and max. 
           ///
           /// Implement get and set (with check if value is in range and throw) but NOT loadValue and saveValue. Do not check the DefValue. 
  template <typename Num>
class CParamBNRange: public IParam, public NumRang<Num>
{
    Num  _v, &_value;
 public:
								/// It accepts a parameter and therefore does not use _v. For compatibility.
    CParamBNRange (IProg *pp, const std::string& titel, const std::string& etiq, Num &parRef, 
						Num min, Num max, Num defValue,
						const std::string& unit=""
					) : IParam (pp, titel, etiq, unit), 
					    NumRang<Num>(min,max), 
						_value(parRef)
	          { /*if (!inRang(defValue)) 
			        throw ParamOutOfNumRange(string("Default Value out of Range while trying to construct: ")+Titel() ); */
	             _value=defValue ;
	          }

								/// Num &parRef,   _v used and therefore does not need an external parametr
    CParamBNRange (IProg *pp, const std::string& titel, const std::string& etiq, 
						Num min, Num max, Num defValue,
						const std::string& unit=""
					) : CParamBNRange ( pp,  titel,  etiq, _v, min,  max,  defValue, unit )  
	          {  
	          }
	void set(Num value){ if (!inRang(value)) ///  \todo: Incluir value y rang que no concuerdan en mensaje to throw
		                    throw OutOfNumRange(std::string("Value out of Range while trying to set: ")+Titel(), value, *this );
	                     if (_value==value) return; 
						 _value=value ; 
						 changed();
	                    }
	Num get()const{return _value;
    }
	//virtual ostream	    &saveValue	(ostream	&osPr	) const  // =0;   ??No salva nada, no tiene "value" todavia
};


class CProject;
//class CEspProg ;

/**  \brief   Abstract "Interface" for parametrs of programs which only save/load project and run programs (all virtual)
 *   \todo Para crear y anadir un nuevo programa:
 *		- General plan: crear (plan) interfase de usuario para tener idea de los parametros a usar
 *		- Common Programm Parametrs: crear class tomando como base Programs::CCommProgParam o uno de sus derivados para contener todos los parametros comunes
 *		- Specific programs: Para cada tarea crear class tomando como base CEspProg o uno de sus derivados para contener los parametros especificos
 *      - Project: crear class tomando como base Programs::CProject y anadirle todos los Programs::IProg anteriores
 *      - Parametrs: Anadir a cada uno de los anteriores Programs::IProg los parametr especif (como se vio en la UI) e inicializarlos 
 *      - New Parametrs: Si es necesario implementar derivado de Programs::IParam:
 *		    - implementar funciones load/save del project. Load: con etiqueta para cada param, "	<< boolalpha " para bool, 
 *		    - implementar funciones de actualizar parametros <=> interfase de usuario : UpdateThDyP() & UpdateThDyForm()
*		- crear funcion del programa en su propio .cpp que toma como argumento un puntero a class C, call it in a overloaded Programs::IProg::Run 
 */
class IProg : public IBParam // -------	  Clase    ----------
{   
 public:
    std::map<std::string,IParam*> _parametrs;
    void insertParam(IParam *pp){_parametrs[pp->Etiq ()]=pp;}

	IProg (const std::string& titel, CProject *proj=nullptr); /*:_Titel(titel){ if (proj) proj->_ProgList.push_back(this);}*/
	std::ofstream	&save		(std::ofstream	&osPr				 )  const
	                            {   
                                    osPr << std::endl <<"\t------\t"<<Titel()<<" "<<std::endl ;
									for (auto &par : _parametrs) 
								        par.second->save(osPr); 
	                                return osPr;
	                            }
	         bool	load		(std::string		&etiq, std::ifstream &isPr) 
	                            {   
                                    auto p=_parametrs.find(etiq); 
	                                if (p==_parametrs.end()) 
										return false;
	                                return p->second->load(isPr); //throw execption if false ????
	                            }
	virtual	int		Run			(IProg &prog				){return prog.Run();}       //  ???????
	virtual int		Run			(		void					)
	                            {   
                                    for(int WorkToDo=Initialize(); WorkToDo>0 ; WorkToDo=Continue()) 
										CallBack(WorkToDo); 
									return Finalize();			 
	                            } 
	virtual int		Initialize	(		void					){ return 0;} 
	virtual int		Continue	(		void					){ return false;}
	virtual int		Finalize	(		void					){ return 0;} 
	virtual void	CallBack	(		int WorkToDo			){}
	virtual ~IProg() override{}
};	

//typedef CEspProg *pCEspProgParam ;
/// Clase base para los parametros "Especificos" de programas "Especificos".
/// derivar para concretar parametros comunes. Mantiene link a proj de los prog Espec que los usan.
/// 
//How to use? 
// Each parameter have an unique identificator or etiquette. 
// While loading, the text between the beginning of a line and the first : will be taken as
// an etiquette (discarding surrounding but not internal spaces). 
//If the etiquette is known (valid), the rest of the line will be use to deduce the value of the parameter. 
//Some parameter (like names of files) will assume this rest-line-text entirely as his valid value. 
//For such parameter please, add any comment in the next line. 
//Others parameter (like numeric or bool parameters) will only use the beginning of this rest-line-text 
//and will ignore the end. 
//Any line without a valid etiquette will be ignore (they are comments!).” 
//Only the last valid value of each parameter will be used.
//For not defined parameters, the previous value (from the previously active project 
//or from the program´s default) will be use.
//Send questions please to: ArielVina.Rodriguez@fli.bund.de

class CProject : public IProg
{
		std::string		    _defPr ;
		std::string		    _ProjetFileName ;
	    std::vector<IProg*> _ProgList;
public:
	CProject(std::string titel, std::string	prFname="", std::string	defProFN="Def.Proj.txt")
		:   IProg           (titel), 
            _defPr          (defProFN)  ,
		    _ProjetFileName (prFname.empty () ? _defPr :  prFname )  
	{  
	} 

   ~CProject()override { }

	void	    ProjetFile	(const std::string &ProjetFileName){	_ProjetFileName=trim_string(ProjetFileName);	}
	std::string ProjetFile	(                            )const{	return _ProjetFileName;	}

    bool	         load		(); 
	bool	         load_defPr	()                                 {  ProjetFile(_defPr);         return load();	}
	bool	         load	    (const std::string &ProjetFileName){  ProjetFile(ProjetFileName); return load();	}

    ///  \todo  Este es el verdadero save !!! El que abre el fichero y lo salva todo.
	std::ofstream	&saveToFile	(const char *ProjetFileName) const{	std::ofstream osPr(ProjetFileName);			return save_all(osPr);}

	std::ofstream	&save		()		const		            {	return saveToFile(_ProjetFileName.c_str())	;   }
	std::ofstream	&save_defPr	()                              {   ProjetFile(_defPr);         
                                                                    return save();	    }
	std::ofstream	&save_asDefPr()		const		            {	return saveToFile(_defPr.c_str())	        ;   }
	std::ofstream   &save	 (const std::string &ProjetFileName){	ProjetFile(ProjetFileName); 
                                                                    return save();	    }

	virtual std::ofstream &saveTMP() const            /// Reescribe el projecto actual. Pensar algo mejor? Preguntar al user? usar # conscuti?
	                            {	return save();	}

	   std::ofstream&	save_all(std::ofstream &osPr)	const 			
	{   
        for(auto p : _ProgList) 
			p->save(osPr) ;		
	        IProg::save(osPr) ;   
	   
	   osPr<< std::endl<<std::endl<<
			 "How to use? \n Each program´s parameter have an unique identificator or etiquette. \n "
			 "While loading, the text between the beginning of a line and the first : will be taken as\n "
			 "an etiquette (discarding surrounding but not internal spaces). \n"
			 "IF the etiquette is know (valid), the rest of the line will be use to deduce the value of the parameter. \n"
			 "Some parameter (like file´s names) will assume this rest-line-text entirely as his valid value. \n"
			 "For such parameter please, add any comment in the next line. \n"
			 "Other parameter (like numeric or bool parameters) will only use the beginning of this rest-line-text and ignore the end. \n"
			 "Any line without a valid etiquette will be ignore (they are comments!).” \n"
			 "Only the last valid value of each parameter will be used\n"
			 "For not defined parameters, the previous value (from the previously active project or from the program´s default) will be use.\n"
			 "Direct questions please to ArielVina.Rodriguez@fli.bund.de\n"
			;

	   return (osPr) ;
	   
	   }   // por que solo funciona con el IProg:: ???
	bool	            load_all(std::string &etiq, std::ifstream &isPr)	//override
	{   
        for(auto p : _ProgList)	
			if ( p->load(etiq, isPr)) 
				return true ;
		return IProg::load(etiq, isPr);					 
    }

    int		Run (IProg &prog)	override                    //   ??????
	{	
    saveTMP( ) ; 
	return prog.Run();
	}

	void AddProg (IProg* par) {_ProgList.push_back(par);}

};

class CCommProgParam : public IProg 
{	CProject *_proj;
 public:	
	CCommProgParam  (const std::string& titel,       CProject *proj=nullptr)
		            : IProg(titel,proj),   _proj(proj) {}
	~CCommProgParam() override	{}

	std::ofstream	&save_all(std::ofstream	&osPr 				 ) const
	                    {   
							assert(("Atemt to use an unitialized project pointer in save_all",_proj));
							return _proj->save_all(osPr);
	                    }
	bool	    load_all(std::string     &etiq, std::ifstream &isPr)
	                    {
							assert(("Atemt to use an unitialized project pointer in load_all",_proj));
							return _proj->load_all(etiq,isPr);
	                    } 
	void        AddProgToProject(IProg *p)
	                    {
							assert(("Atemt to use an unitialized project pointer in AddProgToProject",_proj));
							_proj->AddProg(p);
	                    }
    virtual std::string  MakeRuningName()const {return "";}

};
class	CEspProg  : public IProg 
{public:										// Permite no duplicar los parametros comunes en los parametros especificos
	explicit CEspProg(const std::string& titel, CCommProgParam &commParam ) 
		                    : _cp(commParam), IProg(titel ) 
	                        { _cp.AddProgToProject(this);}
	CCommProgParam &_cp;
	std::ofstream	&save_all(std::ofstream &osPr)	const		   // Save all needed parametrs for this programm, not only the spesific ones
	                     {          _cp.save(osPr) ;			
	                         return     save(osPr) ;    }
	bool	    load_all(std::string &var, std::ifstream &isPr)	  // Usar estas dos funciones solo si se quiere save or load olny this program
	                    { 	if ( _cp.load(var, isPr))       // PAra salvar el projecto completo use el save_all del projecto
						        return true ;
							return load(var, isPr);					 }
} ;

}
//using namespace Programs ;
#endif


