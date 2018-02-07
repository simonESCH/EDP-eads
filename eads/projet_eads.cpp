/**
 * \file md_projet.cpp
 */


//#include "md_projet.hpp"
//#include "md_ns.hpp"
//#include "md_heat.hpp"
//#include "projet_eads.hpp"
#include "projet_eads2.hpp"





/// \fn main
int main(int argc, char* argv[])
{
    Environment env( 
            _argc= argc, 
            _argv= argv, 
            _desc= makeOptions(), 
            _about= about( _name= "projetEDP-M1-CSMI", 
                _author= "Simon ESCHACH", 
                _email= ""));


    //affichage();

    tic();
#if TEST_PROJET_HPP
    Projet eads;
    eads.run();
#else
    auto mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>, 
            _desc= createGMSH(), 
            _update= MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES
            );
    toc("load mesh  ");

    auto modele= init_modele();
    NavierStokes projet_NS(mesh, modele);
    Heat projet_Heat(mesh, modele);


    std::string beta_str= init_edge_in();
    Feel::cout<<"\nbeta est alors : \n\tu|in= "<<beta_str<<"\n";

    auto beta_expr= expr<FEELPP_DIM, 1>(beta_str);
    auto Q= expr(soption("Proc.Q"));

    auto NS_to_Heat= opInterpolation(
            _domainSpace= projet_NS.m_fluidt.element<0>().functionSpace(), 
            _imageSpace= projet_Heat.beta.functionSpace()//, 
            //_range= markedelements(projet_Heat.m_mesh, "AIR")
            );

    double t= DBL_MAX;
    Q.setParameterValues({{"t", t}});
    beta_expr.setParameterValues({{"t", t}});


    projet_NS.run(beta_expr);

    NS_to_Heat->apply(
            projet_NS.m_fluidt.element<0>(), 
            projet_Heat.beta);

    projet_Heat.run(Q);
    
    auto exp= exporter(_mesh= mesh, _name= soption("Exporter.save"));
    exp->add("heat", projet_Heat.ut);
    exp->add("param_capacity", projet_Heat.m_rc);
    exp->add("param_conductivity", projet_Heat.m_k);
    exp->add("fluid_velocity", projet_NS.m_fluidt.element<0>());
    exp->add("fluid_pressure", projet_NS.m_fluidt.element<1>());
    exp->save();

    Feel::cout
        <<"la temperature de IC1 est      "<<projet_Heat.get_heat_IC1()<<" K\n"
        <<"la temperature de IC2 est      "<<projet_Heat.get_heat_IC2()<<" K\n"
        <<"la temperature a la sortie est "<<projet_Heat.get_heat_out()<<" K\n";

#endif
    return 0;
}

