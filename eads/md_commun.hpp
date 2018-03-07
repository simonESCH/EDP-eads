#ifndef __MD_COMMUN_HPP__
#define __MD_COMMUN_HPP__ 1

//! file md_commun.hpp

// use in the child file
#include <feel/feel.hpp>


// use in this file
//#include <feel/option.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace Feel;
using namespace vf;

enum Modele_type{modele0=0, modele1, modele2, modele3};
typedef Simplex<FEELPP_DIM> MyMesh_type;
#define TEST_PROJET_HPP 1


/// \fn makeOptions
/// \brief declaration of the option 
inline
    po::options_description
makeOptions()
{
    po::options_description myappOptions( "My app options" );
    myappOptions.add_options()
        ("Air.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of air")
        ("Proc.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of processor")
        ("PCB.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of the print circuit board")

        ("Air.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of air")
        ("Proc.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of processor")
        ("PCB.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of the print circuit board")

        ("Air.mu", po::value<double>()->default_value( 1.8e-5), 
         "viscosity of Air")
        ("Air.rho", po::value<double>()->default_value( 1.184), 
         "density of Air")

        ("Modele.Tamb", po::value< double>()->default_value(293.15), 
         "temperature of reference")
        ("Proc.Q1", po::value< std::string>()->default_value("1e6"), 
         "Quantity of heat of the processor IC2")
        ("Proc.Q2", po::value< std::string>()->default_value("1e6"), 
         "Quantity of heat of the processor IC1")
        ("Proc.Q", po::value< std::string>()->default_value("1e6"), 
         "Quantity of heat of the processor IC1")
        ("Modele.flux", po::value< std::string>()->default_value("0."), 
         "profile du flux d'air a l'entree")
        ("Modele.modele", po::value< std::string>()->default_value("modele0"), 
         "choix du Modele")

        ("Air.D", po::value<std::string>()->default_value("1e3"), 
         "power of the air flow")
        ("Modele.epsilon", po::value<double>()->default_value(1), 
         "epsilon de GaLS")

        //("Time.Tfinal", po::value< double>()->default_value(0), 
        // "temps maximal")
        //("Modele.dt", po::value< double>()->default_value(.1), 
        // "pas de temps")
        ("Modele.GaLS", po::value<bool>()->default_value( false ), 
         "activation de la stabilisation GaLS")
#if 1
        ("Geo.ePCB", po::value<double>()->default_value(2e-3), 
         "thickness of the print circuit board")
        ("Geo.hPCB", po::value<double>()->default_value(13e-2), 
         "height of the print circuit board")
        ("Geo.eAIR", po::value<double>()->default_value(4e-3), 
         "thickness of the areation")
        ("Geo.eIC", po::value<double>()->default_value(2e-3), 
         "thickness of the processors")
        ("Geo.hIC", po::value<double>()->default_value(2e-2), 
         "height of the processors")
        ("Geo.h1", po::value<double>()->default_value(2e-2), 
         "position of the processor IC1")
        ("Geo.h2", po::value<double>()->default_value(7e-2), 
         "position of the processor IC2")
        ("Time.time", po::value<bool>()->default_value(false), 
         "choisis entre avec dependance en temps ou non")
        ("Time.dt", po::value<double>()->default_value(.1), 
         "increment de temps")
        ("Time.save", po::value<double>()->default_value(-1), 
         "increment de temps into two save")
        ("Time.Tfinal", po::value<double>()->default_value(1), 
         "temps de la simulation")

#endif
        ("Exporter.save", po::value< std::string>()->default_value("Solve"), 
         "nom du fichier .case")
        //("Exporter.load", po::value< std::string>()->default_value(""), 
        // "chemin de l'approximation fine")
        ;
    myappOptions.add(backend_options("backend_heat"));
    myappOptions.add(backend_options("backend_fluid"));
    return myappOptions.add( feel_options() );
}





/// \fn affichage
/// \brief affichage of the option in the file of config
//void affichage()
//{
//
//    Feel::cout<< "\n\t========================"
//        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/m^3"
//        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/m^3"
//        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/m^3"
//        << "\n\t|======================="
//        << "\n\t|Air.k : " << doption("Air.k") << " W/K/m"
//        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/m"
//        << "\n\t|PCB.k : " << doption("PCB.k") << " W/K/m"
//        << "\n\t|======================="
//        << "\n\t|Tamb  : " << doption("Modele.Tamb") << " K"
//        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/m^3"
//        << "\n\t|ventil : " << "..."
//        << "\n\t|modele : " << soption("Modele.modele");
//
//    //if(soption("Exporter.load").compare(""))
//    //  Feel::cout<< "\n\t|load  : "<< soption("Exporter.load");
//    //if(soption("Exporter.save").compare(""))
//    //  Feel::cout<< "\n\t|save  : "<< soption("Exporter.save");
//
//    Feel::cout<< "\n\t|======================="
//        << "\n\t|hsize : " << doption("gmsh.hsize")*1e3<< " mm"
//        << "\n\t|Tmax  : " << doption("Time.Tfinal")<< " sec"
//        << "\n\t|dt   : " << doption("Time.dt")<< " sec"    
//        << "\n\t========================\n\n";
//}



/// \fn init_choice
/// \brief set the using model
Modele_type init_modele()
{
    std::string choice=soption("Modele.modele");
    if(!choice.compare("modele1"))
        return modele1;
    else if(!choice.compare("modele2"))
        return modele2;
    else if(!choice.compare("modele3"))
        return modele3;
    else if(!choice.compare("modele0"))
        return modele0;
    else
    {
        Feel::cout<< "attention : le modele choisi n'est pas repertorie\n\tle modele choisi sera par defaut\n";
        return modele0;
    }
}






std::string init_edge_in()
{
    // recupere la force du debit d'air
    std::string D= soption(_name="Air.D");
    size_t n_d= D.find(":t");

    // recuperation des longueurs
    double eIC= doption("Geo.eIC");
    double eAir= doption("Geo.eAIR");
    double ePCB= doption("Geo.ePCB");
    //modele choisi
    Modele_type modele= init_modele();

    // parametre du profile de poiseuille
    double epaisseur, milieu;
    if ((modele == modele1) || (modele == modele0))
    {
        epaisseur= eAir-eIC;
        milieu= (eAir+eIC)/2. + ePCB;
    }
    else
    {
        epaisseur= eAir;
        milieu= eAir/2. + ePCB;
    }

    // creation de la chiane de caractere
    std::ostringstream ostr;
    ostr
        <<std::scientific
        << "{0, 3/2*"<< D.substr(0, n_d)<< "/("<< epaisseur<<")"
        << "*(1-("
        << "(x-"<< milieu<< ")/("<< epaisseur<< "/2)"
        << ")^2)*"
        << "(x>"<< milieu-epaisseur/2<< ")*"
        << "(x<"<< milieu+epaisseur/2<< ")";
    ostr<< "}:x:t";
    if((n_d <D.size())&&(! boption("Time.time")))
        Feel::cerr<<"\n\t  WARNING!!! : D is dependent of the time\n";
    return ostr.str();
}


    gmsh_ptrtype
createGMSH()
{
    gmsh_ptrtype desc(new Gmsh);

    std::ostringstream ostr_desc;
    std::string filename_path= Environment::findfile(soption("gmsh.filename"));
    Feel::cout << "description 1\n" << desc->getDescriptionFromFile( filename_path ) << "\n";
    
    

#if(FEELPP_DIM==2)
    ostr_desc
#if 1// veritable test
        << "//create by md_commun.hpp\n"
        << "//parameter of the mesh\n"
        << "\n"
        << "h    = "<< doption("gmsh.hsize")<< ";//m\n"
        << "\n"
        << "// dimension of the motherboard( PCB)\n"
        << "hPCB = "<< doption("Geo.hPCB")<< ";//m\n"
        << "ePCB = "<< doption("Geo.ePCB")<< ";//m\n"
        << "\n"
        << "// dimension of the processor ( IC)\n"
        << "hIC = "<< doption("Geo.hIC")<< ";//m\n"
        << "eIC = "<< doption("Geo.eIC")<< ";//m\n"
        << "\n"
        << "// placement of the processor\n"
        << "h1   = "<< doption("Geo.h1")<< ";//m\n"
        << "h2   = "<< doption("Geo.h2")<< ";//hPCB-hIC-h1;//m\n"
        << "\n"
        << "//dimension of the areation( A)\n"
        << "eAIR = "<< doption("Geo.eAIR")<< ";//m\n\n"
        //<< "Include \"eads.geo\";\n";

        << desc->preamble()<< "\n"





        << "//______point__________________________________________\n"
        << "// corner of the motherboard\n"
        << "Point(1)    = {0, 0, 0, h};\n"
        << "Point(2)    = {ePCB, 0, 0, h};\n"
        << "Point(3)    = {ePCB, h1, 0, h};\n"
        << "Point(4)    = {ePCB, h1+hIC, 0, h};\n"
        << "Point(5)    = {ePCB, h2, 0, h};\n"
        << "Point(6)    = {ePCB, h2+hIC, 0, h};\n"
        << "Point(7)    = {ePCB, hPCB, 0, h/2};\n"
        << "Point(8)    = {0, hPCB, 0, h};\n"
        << "\n"
        << "// other corner of processors\n"
        << "Point(9)    = {ePCB+eIC, h1, 0, h};\n"
        << "Point(10)   = {ePCB+eIC, h1+hIC, 0, h};\n"
        << "Point(11)   = {ePCB+eIC, h2, 0, h};\n"
        << "Point(12)   = {ePCB+eIC, h2+hIC, 0, h};\n"
        << "\n"
        << "// other corner of the areation\n"
        << "Point(13)   = {ePCB+eAIR, 0, 0, h};\n"
        << "Point(14)   = {ePCB+eAIR, hPCB, 0, h};\n"
        << "Point(15)   = {ePCB+eIC, 0, 0, h};\n"
        << "\n"
        << "//______line __________________________________________\n"
        << "Line(1) = {1, 2};\n"
        << "Line(2) = {2, 3};\n"
        << "Line(3) = {3, 4};\n"
        << "Line(4) = {4, 5};\n"
        << "Line(5) = {5, 6};\n"
        << "Line(6) = {6, 7};\n"
        << "Line(7) = {7, 8};\n"
        << "Line(8) = {8, 1};\n"//PBC's cycle

        << "Line(9) = {3, 9};\n"
        << "Line(10) = {9, 10};\n"
        << "Line(11) = {10, 4};\n"//IC1's cycle

        << "Line(12) = {5, 11};\n"
        << "Line(13) = {11, 12};\n"
        << "Line(14) = {12, 6};\n"//IC2's cycle

        << "Line(15) = {2, 15};\n"
        << "Line(16) = {15, 13};\n"
        << "Line(17) = {13, 14};\n"
        << "Line(18) = {14, 7};\n"
        << "\n"
        << "//______surface________________________________________\n"
        << "Line Loop(20) = {1, 2, 3, 4, 5, 6, 7, 8};\n"
        << "Plane Surface(100) = {20};//motherBoard\n"
        << "Line Loop(21) = {9, 10, 11, -3};\n"
        << "Plane Surface(101) = {21};//processor IC1\n"
        << "Line Loop(22) = {12, 13, 14, -5};\n"
        << "Plane Surface(102) = {22};//processor IC2\n"
        << "Line Loop(23) = {15, 16, 17, 18, -6, -14, -13, -12, -4, -11, -10, -9, -2};\n"
        << "Plane Surface(103) = {23};//the conduct of areation\n"
        << "\n"
        << "//______physical-line__________________________________\n"
        << "Physical Line(\"in1\") = {15};//condition of dirichlet: T=T0, u_a=f_entre\n"
        << "Physical Line(\"in2\") = {16};//condition of dirichlet: T=T0, u_a=f_entre\n"
        << "Physical Line(\"out\") = {18};//condition of Robin\n"
        << "Physical Line(\"wall\") = {17};// paroi du conduit d'areation\n"
        << "Physical Line(\"borderFluid\") = {2,9,10,11,4,12,13,14,6};// paroi du conduit d'areation interieur\n"
        << "\n"
        << "//______physical-surface_______________________________\n"
        << "Physical Surface(\"PCB\") = {100};//motherboard PCB\n"
        << "Physical Surface(\"IC1\") = {101};//processor IC1\n"
        << "Physical Surface(\"IC2\") = {102};//processor IC2\n"
        << "Physical Surface(\"AIR\") = {103};//conduct of areation\n";
#else
 << "Point(1) = {0, 0, 0, h};\n"
 << "Point(2) = {2, 0, 0, h};\n"
 << "Point(3) = {2, 1, 0, h};\n"
 << "Point(4) = {2, 2, 0, h};\n"
 << "Point(5) = {0, 2, 0, h};\n"
 << "Point(6) = {3, 0, 0, h};\n"
 << "Point(7) = {3, 1, 0, h};\n"
 << "Point(8) = {3, 2, 0, h};\n"
 << "Point(9) = {4, 0, 0, h};\n"
 << "Point(10) = {4, 2, 0, h};\n"
 << "Line(1) = {1, 2};\n"
 << "Line(2) = {2, 3};\n"
 << "Line(3) = {3, 4};\n"
 << "Line(4) = {4, 5};\n"
 << "Line(5) = {5, 1};\n"
 << "Line(6) = {2, 6};\n"
 << "Line(7) = {6, 7};\n"
 << "Line(8) = {7, 3};\n"
 << "Line(9) = {7, 8};\n"
 << "Line(10) = {8, 4};\n"
 << "Line(11) = {6, 9};\n"
 << "Line(12) = {9, 10};\n"
 << "Line(13) = {10, 8};\n"
 << "Line Loop(15) = {1, 2, 3, 4, 5};\n"
 << "Plane Surface(15) = {15};\n"
 << "Line Loop(17) = {6, 7, 8, -2};\n"
 << "Plane Surface(17) = {17};\n"
 << "Line Loop(19) = {8, 3, -10, -9};\n"
 << "Plane Surface(19) = {19};\n"
 << "Line Loop(21) = {9, -13, -12, -11, 7};\n"
 << "Plane Surface(21) = {21};\n"

         << "//______physical-line__________________________________\n"
        << "Physical Line(\"in1\") = {2};//condition of dirichlet: T=T0, u_a=f_entre\n"
        << "Physical Line(\"in2\") = {3};//condition of dirichlet: T=T0, u_a=f_entre\n"
        << "Physical Line(\"out\") = {12};//condition of Robin\n"
        << "Physical Line(\"wall\") = {1,4};// paroi du conduit d'areation\n"
        << "Physical Line(\"borderFluid\") = {5};// paroi du conduit d'areation interieur\n"
        << "\n"
        << "//______physical-surface_______________________________\n"
        << "Physical Surface(\"PCB\") = {21};//motherboard PCB\n"
        << "Physical Surface(\"IC1\") = {19};//processor IC1\n"
        << "Physical Surface(\"IC2\") = {17};//processor IC2\n"
        << "Physical Surface(\"AIR\") = {15};//conduct of areation\n";

#endif

    //Feel::cout<< "dimension de la geometrie:\n"<< ostr_desc.str()<< "\n\n";
#else
    Feel::cerr<< "le programme ne gere pas les dimensions autre que 2\n";
    exit();
#endif

    std::ostringstream nameStr;
    nameStr<< "geo_eads";

    desc->setPrefix(nameStr.str());
    desc->setCharacteristicLength(doption("gmsh.hsize"));
    desc->setDescription(ostr_desc.str());

    return desc;


}







#endif
