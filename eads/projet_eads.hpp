
#define __MD_PROJET_HPP__


/**
 *\file md_projet.hpp
 */

#include <feel/feel.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 
#include "md_ns.hpp"
#include "md_heat.hpp"

using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> M_type;
//typedef Feel::Expr

class Projet
{
    //le but de cette classe sera d'initialiser les deux autres classes
    // faire plusieur mode en fonction du m_modele 
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






Projet::Projet():m_heat(),m_ns()
{
    Feel::cout << "initialisation\n";

    tic();
    m_modele= init_modele();
    m_poiseuille= init_edge_in();
    
    affichage();//mis ici car, si erreur, pas attendre de faire m_mesh

    m_mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>,
            _desc=createGMSH()
            );

    // si le m_modele n'est pas aere ou sans l'equation de Stokes 
    // m_ns n'est pas initialise
    if((m_modele == modele2)||(m_modele == modele3))
        //m_ns= NavierStokes(m_mesh, m_modele);
        m_ns.init(m_mesh, m_modele);
    //m_heat= Heat(m_mesh, m_modele);
    m_heat.init(m_mesh, m_modele);
    toc("init system");

}

void Projet::run()
{
    if(boption("Time.time"))
        //run_dynamic();
        Feel::cout<<"Work In Progress\n";
    else
        run_static();
}

void Projet::run_static()
{
    auto expCase= exporter(
            _mesh= m_mesh, 
            _name= soption("Exporter.save")
            );

    auto Q= expr(soption("Proc.Q"));

    // estimation de la convection
    if(m_modele != modele0)
    {
        auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);
        if(m_modele == modele1)
            m_heat.beta.on(
                    _range= markedelements(m_heat.m_mesh, "AIR"), 
                    _expr= beta
                    ); 
        else
        {//avec equation des m_fluides
            auto flowToConv= opInterpolation(
                    _domainSpace= m_ns.m_fluidt.element<0>().functionSpace(), 
                    _imageSpace= m_heat.beta.functionSpace()
                    );

            m_ns.run(beta);
            flowToConv->apply(
                    m_ns.m_fluidt.element<0>(), 
                    m_heat.beta
                    );
        }
    }

    m_heat.run(Q);

    expCase->add("heat", m_heat.ut);

    // ajout de la convection dans le .case
    if(m_modele == modele1)
        expCase->add("convection", m_heat.beta);
    else if(m_modele != modele0)
    {
        expCase->add("m_fluid_pressure", m_ns.m_fluid.element<1>());
        expCase->add("m_fluid_velocity", m_ns.m_fluid.element<0>());
    }

    expCase->save();
    Feel::cout
        << "|> temperature of IC1 : " << m_heat.get_heat_IC1() << "\n"
        << "|> temperature of IC2 : " << m_heat.get_heat_IC2() << "\n"
        << "|> temperature of out : " << m_heat.get_heat_out() << "\n"
        << "|>      L(u_sol)      : " << m_heat.get_error_Lu() << "\n";
}












#if 0
void Projet::run_dynamic()
{
    // stockage des temperature au cour du temps
    std::ostringstream ostr_heat;
    ostr_heat
        << std::setw(15) << "time"
        << std::setw(15) << "temp_IC1"
        << std::setw(15) << "temp_IC2"
        << std::setw(15) << "temp_out";

    auto Q= expr(soption("Proc.Q"));
    double dt= doption("Time.dt");
    double Tfinal= doption("Time.Tfinal");


    // cas du m_modele sans refroidissement
    if(m_modele == modele0)
    {
        for(double t= 0; t<Tfinal; t+= dt)
        {
            Q.setParameterValues({{"t", t}});
            m_heat.run_step(Q, dt);

            ostr_heat << std::scientific
                << std::setw(15) << t
                << std::setw(15) << m_heat.get_heat_IC1()
                << std::setw(15) << m_heat.get_heat_IC2()
                << std::setw(15) << m_heat.get_heat_out();
        }
    }
    else
    {
        auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);

        // cas du m_modele de refroidissement par un profil de poisseuille constant
        // sur la longueur
        if(m_modele == modele1)
        {
            for(double t = 0; t<Tfinal; t+=dt)
            {
                Q.setParameterValues({{"t", t}});
                beta.setParameterValues({{"t", t}});

                m_heat.beta.on(
                        _range= markedelements( m_heat.m_mesh, "AIR"), 
                        _expr= beta
                        );
                m_heat.run_step(Q, dt);

                ostr_heat << std::scientific
                    << std::setw(15) << t
                    << std::setw(15) << m_heat.get_heat_IC1()
                    << std::setw(15) << m_heat.get_heat_IC2()
                    << std::setw(15) << m_heat.get_heat_out();
            }
        }
        else
        {
            // les derniers modeles avec les equations de m_fluide
            auto flowToConv= opInterpolation( 
                    _domainSpace= m_ns.m_fluidt.element<0>().functionSpace(), 
                    _imageSpace= m_heat.beta.functionSpace()
                    );
            for(double t = 0; t<Tfinal; t+=dt)
            {
                Q.setParameterValues({{"t", t}});
                beta.setParameterValues({{"t", t}});
                m_ns.run_step(beta, dt);
                flowToConv->apply(
                        m_ns.m_fluidt.element<0>(), 
                        m_heat.beta
                        );
                m_heat.run_step(Q, dt);

                ostr_heat << std::scientific
                    << std::setw(15) << t
                    << std::setw(15) << m_heat.get_heat_IC1()
                    << std::setw(15) << m_heat.get_heat_IC2()
                    << std::setw(15) << m_heat.get_heat_out();
            }
        }
    }

    Feel::cout
        << "|> temperature of IC1 : " << m_heat.get_heat_IC1() << "\n"
        << "|> temperature of IC2 : " << m_heat.get_heat_IC2() << "\n"
        << "|> temperature of out : " << m_heat.get_heat_out() << "\n"
        << "|>      L(u_sol)      : " << m_heat.get_error_Lu() << "\n";


    std::fstream file("heat_evolution.dat", std::ios::out|std::ios::trunc);
    if(file)
    {
        Feel::cout << "the evolution of the heat is in the file :"
            << "heat_evolution.dat" << "\n";
        file << ostr_heat.str();
        file.close();
    }
    else
        Feel::cerr << "le fichier n'a pas pu s'ouvrir\n";

}
#endif

void Projet::affichage()
{
    Feel::cout
        << "\n\t|== CAPACITY THERMIC ================= "
        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/m^3"
        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/m^3"
        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/m^3"
        << "\n\t|== CONDUCTIVITY THERMIC ============= "
        << "\n\t|Air.k  : " << doption("Air.k") << " W/K/m"
        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/m"
        << "\n\t|PCB.k  : " << doption("PCB.k") << " W/K/m"
        << "\n\t|== FLUID PARAMETER ================== "
        << "\n\t|rho    : " << doption("Air.rho") << " kg/m^3"
        << "\n\t|mu     : " << doption("Air.mu") << " kg/m/s"
        << "\n\t|===================================== "
        << "\n\t|Tamb   : " << doption("Modele.Tamb") << " K"
        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/m^3"
        << "\n\t|ventil : " << m_poiseuille
        << "\n\t|m_modele : " << m_modele
        << "\n\t|===================================== "
        << "\n\t|hsize  : " << doption("gmsh.hsize")*1e3 << " mm"
        << "\n\t|Tmax   : " << doption("Time.Tfinal") << " sec"
        << "\n\t|dt     : " << doption("Time.dt") << " sec"
        << "\n\t|===================================== "
        << "\n\t|stab   : " << std::boolalpha << boption("Modele.GaLS")
        << "\n\t|save   : " << soption("Exporter.save");

}

