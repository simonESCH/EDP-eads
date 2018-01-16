
#define __MD_PROJET_HPP__


/**
 *\file md_projet.hpp
 */

#include <feel/feel.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include "md_ns.hpp"
#include "md_heat.hpp"

using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> M_type;
//typedef Feel::Expr

class Projet:
{
    //le but de cette classe sera d'initialiser les deux autres classes
    // faire plusieur mode en fonction du modele 
    // PENSER A DEPLACER L'AFFICHAGE ICI

    Heat m_heat;
    NavierStokes m_ns;
    
    boost::shared_ptr<Mesh<M_type>> m_mesh;
    Modele_type m_modele;
    std::string m_poiseuille;
    
    public:
    Projet();

    void affichage();
    void run();
    void run_static();
    void run_dynamic();

};

Projet::Projet()
{
    Feel::cout<<"initialisation\n";

    tic();
    modele=init_modele();
    m_poiseuille= init_edges_in();
    m_mesh= loadMesh(_mesh= new M_type);
    
    // si le modele n'est pas aere ou sans l'equation de Stokes 
    // m_ns n'est pas initialise
    if((modele==modele2)||(modele==modele3))
        m_ns= NavierStokes(mesh,modele);
    m_heat= Heat(mesh,modele);
    toc("init system");

    affichage();
    // choix du modele
    // initialisation du profile de poiseuille
    // creation du maillage
    // initialisation de m_ns si le modele est ...
    // initialisation de heat
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
    auto expCase= exporter(
            _mesh=m_mesh,
            _name=doption("Exporter.save")
            );
    
    switch(modele)
    {
        case modele0:
            m_heat.run();
            expCase->add("heat",m_heat.ut);
            break;

        case modele1:
            m_heat.beta.on(
                    _range=markedelements(m_heat.m_mesh,"AIR"),
                    _expr=expr<FEELPP_DIM,1>(m_poiseuille)
                    );
            m_heat.run();
            expCase->add("heat",m_heat.ut);
            expCase->add("convection",m_heat.beta);
            break;

        default:
            auto flowToConv= Opinterpolation(
                    _domainSpace=,
                    _imageSpace=,
                    _range=,
                    );

            m_ns.run();
            flowToConv.apply(
                    m_ns.fluidt.element<0>(),
                    m_heat.beta
                    );
            m_heat.run();

            expCase->add("heat",m_heat.ut);
            expCase->add("fluid_pressure",m_ns.fluid.element<1>());
            expCase->add("fluid_velocity",m_ns.fluid.element<0>());
    }
    expCase->save();
    
    Feel::cout
        <<"temperature of IC1 : "<<m_heat.get_heat_IC1()<<"\n"
        <<"temperature of IC2 : "<<m_heat.get_heat_IC2()<<"\n"
        <<"temperature of out : "<<m_heat.get_heat_out()<<"\n"
        <<"     L(u_sol)      : "<<m_heat.get_error_Lu()<<"\n";


}


void Projet::run_dynamic()
{

}

void Projet::affichage()
{
    Feel::cout
        << "\n\t|== CAPACITY THERMIC ================="
        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/m^3"
        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/m^3"
        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/m^3"
        << "\n\t|== CONDUCTIVITY THERMIC ============="
        << "\n\t|Air.k  : " << doption("Air.k") << " W/K/m"
        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/m"
        << "\n\t|PCB.k  : " << doption("PCB.k") << " W/K/m"
        << "\n\t|== FLUID PARAMETER =================="
        << "\n\t|rho    : " << doption("Air.rho") << " kg/m^3"
        << "\n\t|mu     : " << doption("Air.mu") << " kg/m/s"
        << "\n\t|====================================="
        << "\n\t|Tamb   : " << doption("Modele.Tamb") << " K"
        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/m^3"
        << "\n\t|ventil : " << m_poiseuille
        << "\n\t|modele : " << m_modele
        << "\n\t|====================================="
        << "\n\t|hsize  : " << doption("gmsh.hsize")*1e3<<" mm"
        << "\n\t|Tmax   : " << doption("Time.Tfinal")<<" sec"
        << "\n\t|dt     : " << doption("Time.dt")<<" sec"
        << "\n\t|====================================="
        << "\n\t|stab   : " << boption("Modele.GaLS")
        << "\n\t|save   : " << soption("Exporter.save");

}

