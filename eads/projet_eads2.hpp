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
#include "md_ns2.hpp"
#include "md_heat2.hpp"

using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> M_type;
//typedef Feel::Expr



#if TEST_PROJET_HPP
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
    std::ostringstream m_sortie;

    public:
    Projet();

    void affiche_parametre();
    void affiche_sortie();
    void add_sortie(double t);

    void run();
    void run_static();
    void run_dynamic_fluid();
    void run_dynamic_no_fluid();
    

};






Projet::Projet():m_heat(), m_ns(), m_sortie()
{
    Feel::cout << "initialisation\n";

    tic();
    m_modele= init_modele();
    m_poiseuille= init_edge_in();

    affiche_parametre();//mis ici car, si erreur, pas attendre de faire m_mesh

    m_mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>, 
            _desc= createGMSH()
            );

    Feel::cout << "marqueur mesh\n";

    m_sortie
        << std::setw(15) << "time"
        << std::setw(15) << "temp_IC1"
        << std::setw(15) << "temp_IC2"
        << std::setw(15) << "temp_out"
        << "\n";


        m_heat.init(m_mesh, m_modele);
        if( (m_modele== modele2) || (m_modele== modele3))
            m_ns.init(m_mesh, m_modele);

    toc("init system");

}








void Projet::run()
{
    if(boption("Time.time"))
        //if( (m_modele == modele0) || (m_modele == modele1) )
        //    run_dynamic_no_fluid();
        //else
        //    run_dynamic_fluid();
        Feel::cout<< "on ne fait pas dans le temps\n";
    else
    {
        run_static();
    }

    Feel::cout
        << "\n"
        << "|> temperature of IC1 : " << m_heat.get_heat_IC1() << " K\n"
        << "|> temperature of IC2 : " << m_heat.get_heat_IC2() << " K\n"
        << "|> temperature of out : " << m_heat.get_heat_out() << " K\n"
        << "|>      L(u_sol)      : " << m_heat.get_error_Lu() << "  \n";
}








void Projet::run_static()
{
    auto Q= expr(soption("Proc.Q"));
    // estimation de la convection
    
    tic();
    auto bilinear_heat= form2(_test= m_heat.Th, _trial= m_heat.Th,_matrix= m_heat.matrix);
    auto linear_heat= form1(_test= m_heat.Th,_vector= m_heat.vector);
    m_heat.init_matrix(bilinear_heat,linear_heat);
    

    if(m_modele != modele0)
    { // the implementation of the flow is the heat equation
        auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);

        if(m_modele == modele1)
        {
            m_heat.beta.on(
                    _range= elements(m_heat.m_mesh), //markedelements(m_heat.m_mesh, "AIR"), 
                    _expr= beta
                    ); 
            toc("create flow");
        }
        else
        {//avec equation des m_fluides
            

            
            auto flowToConv= opInterpolation(
                    _domainSpace= m_ns.m_fluidt.element<0>().functionSpace(), 
                    _imageSpace= m_heat.beta.functionSpace()
                    );
            toc("init flow");


            tic();
            if( m_modele== modele2)
            {
                auto bilinear_fluid= form2(_test= m_ns.Vph, _trial= m_ns.Vph,_matrix= m_ns.matrix); 
                auto linear_fluid= form1(_test= m_ns.Vph, _vector= m_ns.vector); 
                m_ns.init_matrix(bilinear_fluid, linear_fluid);
                
                m_ns.run(beta);
            }
            else
            {
                auto bilinear_fluid_tmp= form2(_test= m_ns.Vph, _trial= m_ns.Vph,_matrix= m_ns.matrix); 
                auto linear_fluid_tmp= form1(_test= m_ns.Vph, _vector= m_ns.vector);

                auto bilinear_fluid= form2(_test= m_ns.Vph, _trial= m_ns.Vph); 
                auto linear_fluid= form1(_test= m_ns.Vph); 
                m_ns.init_matrix(bilinear_fluid, linear_fluid);

                auto tmp= m_ns.Vph->element();
                tmp.zero();
                
                bool continu=true;
                for(int i= 0;(i<2) && continu;i++)
                {
                    bilinear_fluid_tmp= bilinear_fluid;
                    linear_fluid_tmp= linear_fluid;
                    
                    double error= m_ns.run_picard(tmp, beta); 
                    Feel::cout<< "error picard : " << error << "\n";
                    if(error < 1e-7)
                        continu=false;
                }
                for(int i =0;(i<7) && continu;i++)
                {
                    bilinear_fluid_tmp= bilinear_fluid;
                    linear_fluid_tmp= linear_fluid;
                    
                    double error= m_ns.run_newton(tmp, beta); 
                    Feel::cout<< "error newton : " << error << "\n";
                    if(error < 1e-7)
                        continu=false;
                }
            }
            flowToConv->apply(
                    m_ns.m_fluidt.element<0>(), 
                    m_heat.beta
                    );
            Feel::cout<< "convection modele2\n";
            toc("create flow");
        }
    }
    else toc("init flow");


    tic();
    m_heat.run(Q);
    toc("Heat");


    tic();
    auto expCase= exporter(
            _mesh= m_mesh, 
            _name= soption("Exporter.save")
            );
    expCase->add("heat", m_heat.ut);

    // ajout de la convection dans le .case
    if(m_modele == modele1)
        expCase->add("convection", m_heat.beta);
    else if(m_modele != modele0)
    {
        expCase->add("fluid_pressure", m_ns.m_fluidt.element<1>());
        expCase->add("fluid_velocity", m_ns.m_fluidt.element<0>());
    }
    expCase->save();
    toc("exporter");


}










void Projet::run_dynamic_no_fluid()
{
    // stockage des temperature au cour du temps

    auto Q= expr(soption("Proc.Q"));
    double dt= doption("Time.dt");
    double Tfinal= doption("Time.Tfinal");

    auto expCase= exporter(
            _mesh= m_heat.m_mesh, 
            _name= soption("Exporter.save")
            );

    auto bilinear_static= form2(_test= m_heat.Th, _trial= m_heat.Th);
    auto linear_static= form1(_test= m_heat.Th);
    m_heat.init_matrix(bilinear_static, linear_static);

    auto bilinear= form2(_test= m_heat.Th, _trial= m_heat.Th, _matrix= m_heat.matrix);
    auto linear= form1(_test= m_heat.Th, _vector= m_heat.vector);

    // cas du m_modele sans refroidissement
    if(m_modele == modele0)
    {
        for(double t= 0; t<Tfinal; t+= dt)
        {
            linear= linear_static;
            bilinear= bilinear_static;

            Q.setParameterValues({{"t", t}});
            m_heat.run(Q, dt);

            expCase->step(t)->add("heat", m_heat.uPrec);
            expCase->save();

            m_heat.uPrec= m_heat.ut;
            add_sortie(t);
        }
    }
    else
    {
        auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);

        // cas du m_modele de refroidissement par un profil de poisseuille constant
        // sur la longueur
        if(m_modele == modele1)
        {
            for(double t = 0; t<Tfinal; t+= dt)
            {
                linear= linear_static;
                bilinear= bilinear_static;

                Q.setParameterValues({{"t", t}});
                beta.setParameterValues({{"t", t}});

                m_heat.beta.on(
                        _range= markedelements( m_heat.m_mesh, "AIR"), 
                        _expr= beta
                        );
                m_heat.run(Q, dt);

                expCase->step(t)->add("heat", m_heat.ut);
                expCase->save();

                m_heat.uPrec= m_heat.ut;
                add_sortie(t);
            }
        }
    }
    affiche_sortie();
}

#if 0
void Projet::run_dynamic_fluid()
{

    auto Q= expr(soption("Proc.Q"));
    double dt= doption("Time.dt");
    double Tfinal= doption("Time.Tfinal");
    double error;

    tic();
    auto lin_fluid= form1(_test= m_ns.Vph,_vector= m_ns.vector);
    auto bil_fluid= form2(_test= m_ns.Vph, _trial= m_ns.Vph, _matrix= m_ns.matrix);

    auto lin_heat= form1(_test= m_heat.Th, _vector= m_heat.vector);
    auto bil_heat= form2(_test= m_heat.Th, _trial= m_heat.Th, _matrix= m_heat.matrix);

    auto lin_fluid_static= form1(_test= m_ns.Vph);
    auto bil_fluid_static= form2(_test= m_ns.Vph, _trial= m_ns.Vph);
    m_ns.init_matrix(bil_fluid_static, lin_fluid_static);

    auto lin_heat_static= form1(_test= m_heat.Th);
    auto bil_heat_static= form2(_test= m_heat.Th, _trial= m_heat.Th);
    m_heat.init_matrix(bil_fluid_static, lin_fluid_static);

    auto tmp= m_ns.Vph->element();
    auto beta= expr<FEELPP_DIM, 1>(m_poiseuille);

    // les derniers modeles avec les equations de m_fluide
    auto flowToConv= opInterpolation( 
            _domainSpace= m_ns.m_fluidt.element<0>().functionSpace(), 
            _imageSpace= m_heat.beta.functionSpace()
            );
    toc("init form");


    auto expCase= exporter(
            _mesh= m_mesh, 
            _name= soption("Exporter.save")
            );

    for(double t = 0; t<Tfinal; t+= dt)
    {
        tic();
        Q.setParameterValues({{"t", t}});
        beta.setParameterValues({{"t", t}});

        tic();
        //tmp.element<0>().on(_range= elements(m_ns.m_mesh), _expr= zero<FEELPP_DIM, 1>());
        tmp= m_ns.m_fluidPrec;
        for(int i= 0;i<2;i++)
        {
            bil_fluid= bil_fluid_static;
            lin_fluid= lin_fluid_static;

            //auto expr_temp= m_ns.m_rho * inner(idv(m_ns.m_fluidPrec.element<0>()), idt(m_ns.m_fluidt.element<0>()))/dt;
            //lin_fluid+= integrate(
            //        _range= elements(m_ns.m_mesh),
            //        _expr= expr_temp
            //        );
            m_ns.build_time_linear(lin_fluid,dt);

            m_ns.run_picard( tmp, beta, dt);
        }

        double error= 1;
        bool test_newton= true;
        for(int i= 0;i<13 && test_newton;i++)
        {
            bil_fluid= bil_fluid_static;
            lin_fluid= lin_fluid_static;

            //auto expr_temp= m_ns.m_rho * inner(idv(m_ns.m_fluidPrec.element<0>()), idt(m_ns.m_fluidt.element<0>()))/dt;
            //lin_fluid+= integrate(
            //        _range= elements(m_ns.m_mesh),
            //        _expr= expr_temp
            //        );
            m_ns.build_time_linear(lin_fluid,dt);

            error= m_ns.run_newton( tmp, beta, dt);// NE MARCHE PAS
            if(error < 1e-8)
                test_newton= false;
            else
                Feel::cout << "--> error= " << error << "\n";
        }
        m_ns.m_fluidPrec= m_ns.m_fluidt;

        expCase->step(t)->add("fluid_velocity", m_ns.m_fluidPrec.element<0>());
        expCase->step(t)->add("fluid_pressure", m_ns.m_fluidPrec.element<1>());

        m_ns.m_fluidPrec= m_ns.m_fluidt;
        toc(" FLUID ");

        tic();
        flowToConv->apply(
                m_ns.m_fluidt.element<0>(), 
                m_heat.beta
                );
        bil_heat= bil_heat_static;
        lin_heat= lin_heat_static;

        m_heat.run( Q, dt);
        m_heat.uPrec= m_heat.ut;

        expCase->step(t)->add("heat", m_heat.uPrec);
        toc(" HEAT  ");
        expCase->save();

        add_sortie(t);

        std::ostringstream ostr_time;
        ostr_time << "time : " << std::setprecision(3) << t << " s";
        toc(ostr_time.str());
        Feel::cout<< "\n";
    }



    affiche_sortie();

}
#endif






void Projet::add_sortie(double t)
{
    m_sortie << std::scientific
        << std::setw(15) << t
        << std::setw(15) << m_heat.get_heat_IC1()
        << std::setw(15) << m_heat.get_heat_IC2()
        << std::setw(15) << m_heat.get_heat_out()
        << "\n";
}


void Projet::affiche_sortie()
{
    if(Environment::isMasterRank())
    {
        std::cout << "contenue du fichier ...\n" << m_sortie.str() << "\n";


        std::fstream file("heat_evolution.dat", std::ios::out|std::ios::trunc);
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
            << "\n\t|rho    : " << doption("Air.rho") << " kg/m^3"
            << "\n\t|mu     : " << doption("Air.mu") << " kg/m/s"
            << "\n\t|Re/D   : " << doption("Geo.hPCB")/doption("Geo.eAIR")*doption("Air.rho")/doption("Air.mu");
            ;
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
            << "\n\t|dt     : " << doption("Time.dt") << " sec";
    Feel::cout
        << "\n\t|===================================== "
        << "\n\t|stab   : " << std::boolalpha << boption("Modele.GaLS")
        << "\n\t|time   : " << boption("Time.time")
        << "\n\t|save   : " << soption("Exporter.save")<< "\n\n";

}
#endif
#endif
