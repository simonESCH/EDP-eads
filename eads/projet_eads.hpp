#ifndef __MD_PROJET_HPP__
#define __MD_PROJET_HPP__ 1


/**
 *\file md_projet.hpp
 */

#include <feel/feel.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 
//#include "md_ns.hpp"
//#include "md_heat.hpp"
#include "chaleur.hpp"
#include "fluide.hpp"

using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> simplex_type;
//typedef Feel::Expr



class Projet
{
    //le but de cette classe sera d'initialiser les deux autres classes
    // faire plusieur mode en fonction du m_modele 
    typedef Mesh<simplex_type> mesh_type;
    typedef boost::shared_ptr<Exporter<mesh_type>> exporter_ptrtype;

    Heat m_heat;
    NavierStokes m_ns;

    boost::shared_ptr<mesh_type> m_mesh;
    exporter_ptrtype m_export;
    Modele_type m_modele;
    std::string m_poiseuille;
    std::ostringstream m_sortie;
    int m_nbSave=0;
    bool is_save(double t,double t_save, double dt);


    public:
    Projet();

    void affiche_parametre();
    void affiche_sortie();
    void add_sortie(double t);
    void add_export(double t);

    void run();
    void run_static();
    void run_dynamic();
    //void run_dynamic_fluid();
    //void run_dynamic_no_fluid();


};// classe Projet



/// \fn name_export
// petite fonction donnant un nom a l'export
std::string name_export()
{
    std::ostringstream ostr_export;
    if(! soption("Exporter.save").compare(""))
    {
        ostr_export << "solve_" << soption("Modele.modele");
        if(boption("Time.time"))
            ostr_export << "_dynamic";
        else
            ostr_export << "_fixed"; 
    }
    else
        ostr_export << soption("Exporter.save");
    return ostr_export.str();
}


// constructeur
Projet::Projet():m_heat(), m_ns(), m_sortie()
{
    Feel::cout << "initialisation\n";

    tic();
    m_modele= init_modele(soption("Modele.modele"));
    m_poiseuille= init_edge_in(m_modele);

    //if( (m_modele == modele3) && !boption("Time.time") )
    //{
    //    Feel::cerr << "le modele3 est forcement dépendant du temps\n";
    //    exit(1);
    //}


    m_mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>, 
            _desc= createGMSH()
            );

    m_export= exporter(_mesh= m_mesh, _name= name_export());

    // si le m_modele n'est pas aere ou sans l'equation de Stokes 
    // m_ns n'est pas initialise
    if((m_modele == modele2)||(m_modele == modele3))
        m_ns.init(m_mesh, m_modele);
    m_heat.init(m_mesh, m_modele);


    affiche_parametre();//mis ici car, si erreur, pas attendre de faire m_mesh
    m_sortie
        << std::setw(15) << "time"
        << std::setw(15) << "temp_IC1"
        << std::setw(15) << "temp_IC2"
        << std::setw(15) << "temp_out"
        << "\n";



    toc("init system");

}








void Projet::run()
{
    if(boption("Time.time"))
        run_dynamic();
    else
        run_static();
}








void Projet::run_static()
{
    auto Q= expr(soption("Proc.Q"));
    auto beta= expr<FEELPP_DIM,1>(m_poiseuille);

    // mise en place d'un flux de refroidissement
    switch(m_modele)
    {
        case modele0:
            break;
        case modele1:
            m_heat.beta.on(
                    _range= elements(m_heat.m_mesh),
                    _expr= beta
                    );
            break;
        default:
            // operateur permettant de transformer le flux d'air en flux de convection pour la chaleur
            tic();
            auto fluidToConv=opInterpolation(
                    _domainSpace= m_ns.m_fluidt.element<0>().functionSpace(),
                    _imageSpace= m_heat.beta.functionSpace()
                    );
            m_ns.run(beta);
            fluidToConv->apply(
                    m_ns.m_fluid.element<0>(), 
                    m_heat.beta
                    );
            toc("run fluid");
    }

    tic();
    m_heat.run(Q);
    toc("run heat");

    add_export(0);

    // affichage des parametres de controle
    Feel::cout
        << "\n"
        << "|> temperature of IC1 : " << m_heat.get_heat_IC1() << " K\n"
        << "|> temperature of IC2 : " << m_heat.get_heat_IC2() << " K\n"
        << "|> temperature of out : " << m_heat.get_heat_out() << " K\n"
        << "|>      L(u_sol)      : " << m_heat.get_error_Lu() << "  \n";

}




std::string timeT(double t)
{
    std::ostringstream ostr;
    ostr << "time: " << std::setw(6) << t << "s";
    return ostr.str();
}

bool Projet::is_save(double t,double t_save, double dt)
{
    if(t+dt>t_save*m_nbSave)
        {
            m_nbSave++;
            return true;
        }
        return false;
}



void Projet::run_dynamic()
{
    auto Q= expr(soption("Proc.Q"));
    double dt= doption("Time.dt");
    double Tfinal= doption("Time.Tfinal");


    switch(m_modele)
    {
        case modele0://pas de refroidissement

            for(double t= 0; t<Tfinal; t+=dt)
            {
                tic();

                // reactualisation des parametres dépendant du temps
                Q.setParameterValues({{"t", t}});

                // resolution a l'instant t
                m_heat.run(Q);


                // recuperation des donnée
                if(is_save(t,doption("Time.save"),dt)) 
                    add_export(t);
                add_sortie(t);

                //std::ostringstream ostr;
                //ostr << "time : " <<std::setw(5) << t;
                toc(timeT(t));//ostr.str());
                Feel::cout << "\n";
            }
            break;

        case modele1:// flux d'air fixe avec un profile de Poisseuille
            {
                auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);
                for(double t=0; t<Tfinal;t+=dt)
                {
                    tic();

                    // reactualisation des parametres dépendant du temps
                    Q.setParameterValues({{"t", t}});
                    beta.setParameterValues({{"t", t}});
                    m_heat.betaUpdate(beta);

                    // resolution a l'instant t
                    m_heat.run(Q);

                    // recuperation des donnée
                if(is_save(t,doption("Time.save"),dt)) 
                        add_export(t);
                    add_sortie(t);

                    toc(timeT(t));
                    Feel::cout << "\n";

                }
                break;
            }
        default:// flux d'air suivant equation de Stokes et Navier-Stokes
            {
                auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);

                auto fluidToConv= opInterpolation( 
                        _domainSpace= m_ns.m_fluid.element<0>().functionSpace(), 
                        _imageSpace= m_heat.beta.functionSpace()
                        );

                for(double t=0; t<Tfinal; t+=dt)
                {
                    tic();

                    // reactualisation des parametre dépendant du temps
                    Q.setParameterValues({{"t",t}});
                    beta.setParameterValues({{"t",t}});

                    // mise en place du flux de refroidissement
                    m_ns.run(beta);
                    fluidToConv->apply(m_ns.m_fluid.element<0>(),m_heat.beta);

                    // resolution de l'équation de la chaleur
                    m_heat.run(Q);

                    //recuperation des données
                    if(is_save(t,doption("Time.save"),dt))
                        add_export(t);
                    add_sortie(t);


                    toc(timeT(t));
                    Feel::cout << "\n";
                }
            }
    }

    Feel::cout
        << "|> temperature of IC1 : " << m_heat.get_heat_IC1() << "\n"
        << "|> temperature of IC2 : " << m_heat.get_heat_IC2() << "\n"
        << "|> temperature of out : " << m_heat.get_heat_out() << "\n"
        << "|>      L(u_sol)      : " << m_heat.get_error_Lu() << "\n";

    affiche_sortie();

}







void Projet::add_sortie(double t)
{
    m_sortie << std::scientific
        << std::setw(15) << t
        << std::setw(15) << m_heat.get_heat_IC1()
        << std::setw(15) << m_heat.get_heat_IC2()
        << std::setw(15) << m_heat.get_heat_out()
        << "\n";
}




void Projet::add_export(double t)
{
    tic();
    m_export->step(t)->add("temperature",m_heat.u);
    switch(m_modele)
    {
        case modele0:
            break;
        case modele1:
            m_export->step(t)->add("convection",m_heat.beta);
            break;
        default:
            m_export->step(t)->add("fluidVelocity",m_ns.m_fluid.element<0>());
            m_export->step(t)->add("fluidPressure",m_ns.m_fluid.element<1>());
    }
    m_export->save();
    toc("export");
}


void Projet::affiche_sortie()
{
    if(Environment::isMasterRank())
    {
        std::cout << "contenue du fichier ...\n" << m_sortie.str() << "\n";

        std::ostringstream ostr_sortie;
        ostr_sortie << soption("Modele.fileHeat") << "_" << soption("Modele.modele");
        std::fstream file(ostr_sortie.str(), std::ios::out|std::ios::trunc);
        if(file)
        {
            std::cout << "the evolution of the heat is in the file :"
                << "heat_evolution.dat" << "\n";
            file << m_sortie.str();
            file.close();
        }
        else
            std::cerr << "le fichier \"" 
                << "heat_evolution.dat"
                << "\" n'a pas pu s'ouvrir\n";
    }
}



void Projet::affiche_parametre()
{
    Feel::cout
        << "\n\t|== CAPACITY THERMIC ================= "
        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/m^3"
        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/m^3"
        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/m^3"
        << "\n\t|== CONDUCTIVITY THERMIC ============= "
        << "\n\t|Air.k  : " << doption("Air.k") << " W/K/m"
        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/m"
        << "\n\t|PCB.k  : " << doption("PCB.k") << " W/K/m";
    if( (m_modele == modele2) || (m_modele == modele3) )
        Feel::cout
            << "\n\t|== FLUID PARAMETER ================== "
            << "\n\t|rho    : " << doption("Fluid.rho") << " kg/m^3"
            << "\n\t|mu     : " << doption("Fluid.mu") << " kg/m/s";
    Feel::cout
        << "\n\t|===================================== "
        << "\n\t|Tamb   : " << doption("Modele.Tamb") << " K"
        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/m^3";
    if(m_modele != modele0)
        Feel::cout
            << "\n\t|ventil : " << m_poiseuille;
    Feel::cout
        << "\n\t|m_modele : modele" << m_modele
        << "\n\t|===================================== "
        << "\n\t|hsize  : " << doption("gmsh.hsize")*1e3 << " mm";
    if(boption("Time.time"))
        Feel::cout
            << "\n\t|Tmax   : " << doption("Time.Tfinal") << " sec"
            << "\n\t|dt     : " << doption("Time.dt") << " sec"
            << "\n\t|Tsave  : " << doption("Time.save") << " sec\t soit "
            << (int)(doption("Time.Tfinal")/doption("Time.save")) << "images";
            
    Feel::cout
        << "\n\t|===================================== "
        << "\n\t|stab   : " << std::boolalpha << boption("Modele.GaLS")
        << "\n\t|time   : " << boption("Time.time")
        << "\n\t|save   : " << m_export->path() << "\n\n";

}
#endif
