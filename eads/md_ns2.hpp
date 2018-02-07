#ifndef __MD_NS_HPP__
#define __MD_NS_HPP__ 1


/**
 *\file md_ns.hpp
 */

//#include <feel/feel.hpp>
//#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <fstream>
//#include <sstream>
//#include <string>
#include "md_commun.hpp"

#define MAX_LOOP_PICARD 10

using namespace Feel;
using namespace vf;

//typedef Feel::Expr


//! class NavierStokes
//! \brief approxime l'ecoulemnt du fluide

class NavierStokes
{
    // typedef m_mesh
    typedef Mesh<MyMesh_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // typedef bases
    typedef Lagrange<1, Scalar> bases_pressure;
    typedef Lagrange<2, Vectorial, Continuous> bases_vect;
    typedef bases<bases_vect, bases_pressure> bases_fluid;

    // typedef functionSpace
    typedef boost::parameter::void_ param_void;
    typedef FunctionSpace< mesh_type, bases<bases_vect> > space_conv_type;
    typedef FunctionSpace< mesh_type, bases_fluid> space_fluid_type;

    // typedef pointer of functionSpace
    typedef boost::shared_ptr<space_conv_type> space_conv_ptrtype;
    typedef boost::shared_ptr<space_fluid_type> space_fluid_ptrtype;

    // typedef of functionSpace's element
    typedef typename space_conv_type::element_type element_conv_type;
    typedef typename space_fluid_type::element_type element_fluid_type;

    // typedef pointer of element
    typedef boost::shared_ptr<element_conv_type> element_conv_ptrtype;
    typedef boost::shared_ptr<element_fluid_type> element_fluid_ptrtype;

    // typedef matrix
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename m_backend::sparse_matrix_type matrix_type;
    typedef typename m_backend::vector_type vector_type;

    // typedef ptr matrix
    typedef typename m_backend::sparse_matrix_ptrtype matrix_ptrtype;
    typedef typename m_backend::vector_ptrtype vector_ptrtype;

    public:



    //parametre
    double m_mu;
    double m_rho;

    // m_modele du m_modele
    Modele_type m_modele= modele0;


    // m_mesh
    mesh_ptrtype m_mesh;



    // function of heat convection
    space_conv_ptrtype Th_vect;

    // function of navier stokes
    space_fluid_ptrtype Vph;
    element_fluid_type m_fluid;
    element_fluid_type m_fluidt;
    element_fluid_type m_fluidPrec;

    //m_backend
    backend_ptrtype m_backend;

    //matrix
    matrix_ptrtype matrix;
    vector_ptrtype vector;




    //function

#if TEST_PROJET_HPP
    NavierStokes(){}
    void init(mesh_ptrtype mesh, Modele_type mod);
    void init_fluid();



    void init_matrix();

    template<typename L>
        void build_time_linear(L & linear, double dt);

    template<typename B, typename V>
        void build_fluid_dynamic_bilinear(B & bilinear, V velocity_test, V velocity_trial, V velocity_loop);

    template<typename myexpr_type>
        void run(myexpr_type D);

    template<typename myexpr_type>
void run_picard(element_fluid_type & tmp, myexpr_type flow, double dt, double maxError= 1e-3);
    
    template<typename myexpr_type>
double run_newton(element_fluid_type & tmp, myexpr_type flow, double dt, double maxError= 1e-3);
    
    template<typename F>
void run(F flow, double dt);

template<typename E>
        void init_beta(E expr);

    double get_error_Lu();

};




#if TEST_PROJET_HPP

void NavierStokes::init(mesh_ptrtype mesh, Modele_type mod)
{
    m_mesh= createSubmesh(mesh, markedelements(mesh, "AIR"));
    m_modele= mod;

    //    init_conv();
    Vph= space_fluid_type::New(_mesh= m_mesh);

    m_fluid= Vph->element();
    m_fluidt= Vph->element();
    m_fluidPrec= Vph->element();


    m_mu= doption("Air.mu");
    m_rho= doption("Air.rho");

    m_fluidPrec.element<0>().on(
            _range= elements(m_mesh), 
            _expr= zero<FEELPP_DIM, 1>()
            );

    m_backend= m_backend(_name= "m_backend");
    matrix= m_backend->newMatrix(Vph, Vph);
    vector= m_backend->newVector(Vph);
    
    matrix->zero();
    vector->zero();
}


// function of the class NavierStokes============================================ 




/// \fn init_matrix
/// \brief add to the bilinear form the static part
/// \f$ += \int_\Omega \bigl(\frac{\nabla u+\nabla u^T}{2}\bigl):\nabla v\f$
/// \f$ += \int_\Omega q\nabla u+p\nabla v\f$
template<typename bilinear_type, typename linear_type>
void NavierStokes::init_matrix(bilinear_type & bilinear, linear_type & linear)
{
    tic();
    bilinear.zero();
    linear.zero();

    auto v= m_fluid.template element<0>();
    auto u= m_fluidt.template element<0>();
    auto q= m_fluid.template element<1>();
    auto p= m_fluidt.template element<1>();

    auto expr_diffusion= inner(grad(v), gradt(u));//inner(sym(gradt(u)), grad(v));
    auto expr_pressure= id(q)*divt(u)-idt(p)*div(v);

    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_mu*expr_diffusion
                +expr_pressure
            );

    if(boption("Time.time"))
        bilinear+= integrate(
                _range= elements(m_mesh), 
                _expr= m_rho/doption("Time.dt")*inner(idt(u), id(v))
                );
    toc("NS static  ");
}


    template<typename linear_type>
void NavierStokes::build_time_linear(linear_type & linear, double dt)
{
    auto v= m_fluid.element<0>();
    auto uPrec= m_fluidPrec.element<0>();
    linear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_rho * inner(id(v), idv(uPrec))/dt
            );
}





    template<typename myexpr_type>
void NavierStokes::run_picard(
        element_fluid_type & tmp, 
        myexpr_type flow, 
        double dt, 
        double maxError
        )
{
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto u_tmp= tmp.element<0>();

    auto linear= form1(_test= Vph, _vector= vector);
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix);


    auto autoconv= trans(idv(u_tmp))*gradv(u);//picard
    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_rho * autoconv * id(v) // CONDITION AU BORD !!!!!!!!!!!!!!!!!!!!!!
            );

    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    bilinear+= on(
            _range= markedfaces(m_mesh, {"borderFluid", "wall"}), 
            _rhs= linear, 
            _element= u, 
            _expr= zero<FEELPP_DIM, 1>()
            );

    m_backend->solve(
            _solution= m_fluidt, 
            _matrix= matrix, 
            _rhs= vector
            );
    u_tmp= u;
}

    template<typename myexpr_type>
double NavierStokes::run_newton(
        element_fluid_type & tmp, 
        myexpr_type flow, 
        double dt, 
        double maxError)
{
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto u_tmp= tmp.element<0>();

    auto linear= form1(_test= Vph, _vector= vector);
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix);

    auto autoconv= trans(idv(u_tmp))*gradt(u);//picard
    //auto autoconv= trans(idv(u_tmp))*gradt(u) + 
    //                trans(idt(u))*gradv(u_tmp);//newton
    //auto autoconv_lin= trans(idv(u_tmp))*gradv(u_tmp);

    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_rho * autoconv * id(v)
            );

    //linear+= integrate(
    //        _range= elements(m_mesh), 
    //        _expr= m_rho * autoconv_lin * id(v)
    //        );


    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    bilinear+= on(
            _range= markedfaces(m_mesh, {"borderFluid", "wall"}), 
            _rhs= linear, 
            _element= u, 
            _expr= zero<FEELPP_DIM, 1>()
            );

    m_backend->solve(
            _solution= m_fluidt, 
            _matrix= matrix, 
            _rhs= vector
            );

    double error= normL2(
            _range= elements(m_mesh), 
            _expr= idv(u_tmp)-idv(u)
            );
    u_tmp= u;
    return error;
}




    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void NavierStokes::run(myexpr_type flow, double dt)
{

    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto uPrec= m_fluidPrec.element<0>();

    auto linear= form1(_test= Vph, _vector= vector);
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix);

    linear+= integrate(// terme de memoire
            _range= elements(m_mesh), 
            _expr= m_rho/dt*inner(idv(uPrec), id(v))
            );

    bilinear+= on(// flow d'entree
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    bilinear+= on(// condition au bord
            _range= markedfaces(m_mesh, {"borderFluid", "wall"}), 
            _rhs= linear, 
            _element= u, 
            _expr= zero<FEELPP_DIM, 1>()
            );

    m_backend->solve(
            _solution= m_fluidt, 
            _matrix= matrix, 
            _rhs= vector
            );
}



    template<typename myexpr_type>
void NavierStokes::run(myexpr_type flow)
{
    auto bilinear= form2(
            _test= Vph, 
            _trial= Vph, 
            _matrix= matrix
            );
    auto linear= form1(
            _test= Vph, 
            _vector= vector
            );

    init_matrix(bilinear, linear);
    run(bilinear, linear, flow);
}








#endif
