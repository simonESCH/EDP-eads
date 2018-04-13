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
        
        // capacite calorifique des differents composant 
        ("Air.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of air")
        ("Proc.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of processor")
        ("PCB.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of the print circuit board")

        // conductivite thermique des differents composants
        ("Air.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of air")
        ("Proc.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of processor")
        ("PCB.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of the print circuit board")

        // quantite de chaleur emit par les processeur
        ("Proc.Q", po::value< std::string>()->default_value("1e6"), 
         "Quantity of heat of the processor IC1")
        
        // parametre du fluide de refroidissement
        ("Fluid.mu", po::value<double>()->default_value( 1.8e-5), 
         "viscosity of Air")
        ("Fluid.rho", po::value<double>()->default_value( 1.184), 
         "density of Air")
        ("Fluid.D", po::value<std::string>()->default_value("1e-3"), 
         "power of the air flow")

        // parametre global du modele
        ("Modele.Tamb", po::value< double>()->default_value(293.15), 
         "temperature of reference")
        ("Modele.modele", po::value< std::string>()->default_value("modele0"), 
         "choix du Modele")
        
        // stabilisation
        ("Modele.GaLS", po::value<bool>()->default_value( false ), 
         "activation de la stabilisation GaLS")
        ("Modele.epsilon", po::value<double>()->default_value(1), 
         "epsilon de GaLS")
        
        // nom du fichier ou l'on stocke les valeurs des parametres de controle durant le temps
        ("Modele.fileHeat", po::value<std::string>()->default_value("HeatEvolution"),
        "name of the file with the value of the control parameter over time")


        // parametre temporel
        ("Time.time", po::value<bool>()->default_value(false), 
         "choisis entre avec dependance en temps ou non")
        ("Time.dt", po::value<double>()->default_value(.1), 
         "increment de temps")
        ("Time.save", po::value<double>()->default_value(-1), 
         "increment de temps into two save")
        ("Time.Tfinal", po::value<double>()->default_value(1), 
         "temps de la simulation")
        
        // parametre de dimension
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

        // nom du repertoire contenant les images de la simulation
        ("Exporter.save", po::value< std::string>()->default_value(""), 
         "nom du fichier .case")
        ;

    myappOptions.add(backend_options("backend_heat"));
    myappOptions.add(backend_options("backend_fluid"));

    return myappOptions.add( feel_options() );
}








/// \fn init_choice
/// \brief set the using model
Modele_type init_modele(std::string choice)
{
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
        Feel::cout<< "attention : le modele choisi n'est pas repertorie\n"
            <<"\tle modele choisi sera par defaut modele0\n";
        return modele0;
    }
}





std::string init_edge_in(Modele_type modele)
{
    // recupere la force du debit d'air
    std::string D= soption(_name="Fluid.D");
    size_t n_d= D.find(":t");

    // recuperation des longueurs
    double eIC= doption("Geo.eIC");
    double eAir= doption("Geo.eAIR");
    double ePCB= doption("Geo.ePCB");

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

    // creation de la chaine de caractere
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







void test_gmsh()
{
    double hPCB = doption("Geo.hPCB"),
           ePCB = doption("Geo.ePCB"),
           hIC  = doption("Geo.hIC"),
           eIC  = doption("Geo.eIC"),
           h1   = doption("Geo.h1"),
           h2   = doption("Geo.h2"),
           eAIR = doption("Geo.eAIR"),
           h    = doption("gmsh.hsize");

    CHECK(h<=hPCB) << "le pas est trop grand:" << h << ">" << hPCB << "\n";
    CHECK(h<=ePCB) << "le pas est trop grand:" << h << ">" << ePCB << "\n";
    CHECK(h<=hIC) << "le pas est trop grand:" << h << ">" << hIC << "\n";
    CHECK(h<=eIC) << "le pas est trop grand:" << h << ">" << eIC << "\n";
    CHECK(h<=h1) << "le pas est trop grand:" << h << ">" << h1 << "\n";
    CHECK(h<=h2) << "le pas est trop grand:" << h << ">" << h2 << "\n";
    CHECK(h<=eAIR) << "le pas est trop grand:" << h << ">" << eAIR << "\n";
    
    CHECK(eIC<eAIR) << "le conduit est trop petit:" << eIC << ">=" << eAIR << "\n";
#if FEELPP_DIM==3

#endif
}


    gmsh_ptrtype
createGMSH()
{
    test_gmsh();
    gmsh_ptrtype desc(new Gmsh);
    desc->setCharacteristicLength( doption("gmsh.hsize") );

    std::ostringstream ostr_desc;
    std::string filename= soption("gmsh.filename");
    
    std::string filenameExpand = Environment::expand(filename);
    fs::path mesh_name=fs::path(Environment::findFile(filenameExpand));
    std::string filename_path= Environment::findFile(mesh_name.string());

    Feel::cout << "file: " << filename_path << "\t issu :" << filename << "\n";
    std::string desc_str = desc->getDescriptionFromFile( filename_path);




#if(FEELPP_DIM==2)
    ostr_desc
        << "//create by md_commun.hpp\n"
        << "//parameter of the mesh\n"
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
        << desc_str;


#else
    Feel::cerr<< "le programme ne gere pas les dimensions autre que 2\n";
    exit();
#endif

    //Feel::cout << "description 1\n" << ostr_desc.str() << "\n";
    std::ostringstream nameStr;
    nameStr<< "geo_eads";

    desc->setPrefix(nameStr.str());
    desc->setDescription(ostr_desc.str());

    return desc;


}


#endif
