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
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;

    // typedef ptr matrix
    typedef typename backend_type::sparse_matrix_ptrtype matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;


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

    //backend
    backend_ptrtype backend_fluid;

    //matrix
    matrix_ptrtype m_matrix;
    vector_ptrtype m_vector;




    //function

#if TEST_PROJET_HPP
    NavierStokes(){}
    void init(mesh_ptrtype mesh, Modele_type mod);
#else
    NavierStokes(mesh_ptrtype mesh, Modele_type mod);
#endif
    void init_fluid();



    void init_matrix();

    template<typename L, typename V, typename P>
        void build_fluid_dynamic_linear(L & linear, V velocity, V velocity_previous, double dt);

    template<typename B, typename V>
        void build_fluid_dynamic_bilinear(B & bilinear, V velocity_test, V velocity_trial, V velocity_loop);

    template<typename myexpr_type>
        void run(myexpr_type D);

    //template<typename myexpr_type>
    //    void run_step(myexpr_type D, double dt);

    //template<typename B, typename L, typename E>
    //    void resolve_fluid(B & bilinear, L & linear, E & flow, double error= 1e-3);

    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void run_picard(bilinear_type & bilinear, linear_type & linear,element_fluid_type & tmp, myexpr_type flow, double dt,double maxError=1e-3);
    
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
double run_newton(bilinear_type & bilinear, linear_type & linear,element_fluid_type & tmp, myexpr_type flow, double dt,double maxError=1e-3);
    
    template<typename B, typename L, typename F>
void run(B & bilinear, L & linear, F flow, double dt=-1,double maxError=1e-3);

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

    backend_fluid= backend(_name= "backend_fluid");
    m_matrix= backend_fluid->newMatrix(Vph, Vph);
    m_vector= backend_fluid->newVector(Vph);
    init_matrix();
}


#else


NavierStokes::NavierStokes(mesh_ptrtype i_mesh, Modele_type mod)
{
    m_mesh= createSubmesh(i_mesh, markedelements(i_mesh, "AIR"));
    m_modele= mod;

    //    init_conv();
    Vph= space_fluid_type::New(_mesh= m_mesh);

    m_fluid= Vph->element();
    m_fluidt= Vph->element();
    m_fluidPrec= Vph->element();


    m_fluidt.element<0>().on(
            _range= elements(m_mesh), 
            _expr= zero<FEELPP_DIM, 1>()
            );

    backend_fluid= backend(_name= "backend_fluid");
    m_matrix= backend_fluid->newMatrix(Vph, Vph);
    m_vector= backend_fluid->newVector(Vph);


    m_mu= doption("Air.mu");
    m_rho= doption("Air.rho");
}
#endif


// function of the class NavierStokes============================================ 




/// \fn init_matrix
/// \brief add to the bilinear form the static part
/// \f$ += \int_\Omega \bigl(\frac{\nabla u+\nabla u^T}{2}\bigl):\nabla v\f$
/// \f$ += \int_\Omega q\nabla u+p\nabla v\f$
void NavierStokes::init_matrix()
{
    tic();
    m_matrix->zero();
    m_vector->zero();
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= m_matrix);

    auto v= m_fluid.template element<0>();
    auto u= m_fluidt.template element<0>();
    auto q= m_fluid.template element<1>();
    auto p= m_fluidt.template element<1>();

    auto expr_diffusion=  inner(grad(v), gradt(u));//sym(gradt(u))
    auto expr_pressure= id(q)*divt(u)-idt(p)*div(v);

    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_mu*expr_diffusion
            +expr_pressure
            );

    if(boption("Time.time"))
        bilinear+= integrate(
                _range= elements(m_mesh),
                _expr= m_rho/doption("Time.dt")*inner(idt(u),id(v))
                );
    toc("NS static  ");
}










    template<typename linear_type, typename velocity_type, typename pressure_type>
void NavierStokes::build_fluid_dynamic_linear(linear_type & linear, velocity_type v, velocity_type v_prec, double dt)
{

}

    template<typename bilinear_type, typename velocity_type>
void NavierStokes::build_fluid_dynamic_bilinear(bilinear_type & bilinear, velocity_type v, velocity_type vt, velocity_type vl)
{

}


/// \brief resolve the navier-stokes with de method of Picard 
//    template<typename bilinear_type, typename linear_type, typename expr_type>
//void NavierStokes::resolve_fluid(bilinear_type & bilinear, linear_type & linear, expr_type & flow, double error)
//{
//    auto bilinear_static= form2(_test= Vph, _trial= Vph, _matrix= m_matrix);
//    auto linear_static= form1(_test= Vph, _vector= m_vector);
//
//    auto v= m_fluid.element<0>();
//    auto u= m_fluidt.element<0>();
//    auto u_tmp= m_fluidPrec.element<0>();
//
//
//    //auto expr_navier= ( trans(idt(u))*gradv(u_tmp) )*id(v);
//    
//    double loop_err= 0;
//    double cpt= 0;
//    auto exp_debug=exporter(_mesh=m_mesh,_name="debug");
//
//    do
//    {
//        tic();
//
//        tic();
//        bilinear= bilinear_static;
//        linear= linear_static;
//        u= m_fluidt.element<0>();
//
//        // partie d'auto-convection u\nabla u
//        auto expr_navier_lin= trans(idv(u_tmp))*gradv(u_tmp);
//        auto expr_navier_bil= trans(idt(u))*gradv(u_tmp) + trans(idv(u_tmp))*gradt(u);
//        //auto expr_navier= ( trans(idv(u_tmp))*gradt(u) )*id(v);
//        //auto expr_navier= ( trans(idv(u_tmp))*gradt(u) + trans(idt(u))*gradv(u_tmp) )/4.;
//        linear+= integrate(
//                _range= elements(m_mesh), 
//                _expr= m_rho * expr_navier_lin * id(v)
//                );
//        bilinear+= integrate(
//                _range= elements(m_mesh), 
//                _expr= m_rho * expr_navier_bil * id(v)
//                );
//
//
//        //condition au bord
//        bilinear+= on(//entre d'air
//                _range= markedfaces(m_mesh, {"in1", "in2"}), 
//                _rhs= linear, 
//                _element= u, 
//                _expr= flow
//                );
//        bilinear+= on(//condition au bord
//                _range= markedfaces(m_mesh, {"borderFluid", "wall"}), 
//                _rhs= linear, 
//                _element= u, 
//                _expr= zero<FEELPP_DIM, 1>()
//                );
//
//        toc("build");
//
//        tic();
//        bilinear.solve(
//                _solution= m_fluidt, 
//                _rhs= linear
//                );
//        toc("solve");
//
//
//        if(m_modele != modele2)
//        {
//            loop_err= normL2(
//                    _range= elements(m_mesh), 
//                    _expr= (idv(u)-idv(u_tmp))
//                    );
//            Feel::cout<<"erreur= "<<loop_err<<"\n";
//            cpt++;
//            m_fluidPrec= m_fluidt;
//        }
//
//        exp_debug->step((double) cpt)->add("velocity",u);
//        std::ostringstream ostr;
//        ostr<<"LOOP : "<<cpt;
//        toc(ostr.str());
//
//    } while((error<loop_err)&&(cpt<MAX_LOOP_PICARD));
//    exp_debug->save();
//}




    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void NavierStokes::run_picard(bilinear_type & bilinear, linear_type & linear,element_fluid_type & tmp, myexpr_type flow, double dt,double maxError)
{
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto u_tmp= tmp.element<0>();

    auto autoconv= trans(idv(u_tmp))*gradt(u);//picard
    bilinear+= integrate(
            _range= elements(m_mesh),
            _expr= m_rho * autoconv * id(v)
            );

    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    bilinear+= on(
            _range= markedfaces(m_mesh, {"borderFluid","wall"}),
            _rhs= linear,
            _element= u,
            _expr= zero<FEELPP_DIM,1>()
            );

    bilinear.solve(
            _solution= m_fluidt, 
            _rhs= linear
            );
    u_tmp=u;
}

    template<typename bilinear_type, typename linear_type, typename myexpr_type>
double NavierStokes::run_newton(bilinear_type & bilinear, linear_type & linear,element_fluid_type & tmp, myexpr_type flow, double dt,double maxError)
{
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto u_tmp= tmp.element<0>();

    auto autoconv= trans(idv(u_tmp))*gradt(u)+trans(idt(u))*gradv(u_tmp);//newton
    auto autoconv_lin= trans(idv(u_tmp))*gradt(u_tmp);
    Feel::cout << "BALISE NEWTON 1\n";
    bilinear+= integrate(
            _range= elements(m_mesh),
            _expr= m_rho * autoconv * id(v)
            );

    Feel::cout << "BALISE NEWTON 2\n";
    linear+= integrate(
            _range= elements(m_mesh),
            _expr=- m_rho * autoconv_lin * id(v)
            );
    
    Feel::cout << "BALISE NEWTON 3\n";
    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    Feel::cout << "BALISE NEWTON 4\n";
    bilinear+= on(
            _range= markedfaces(m_mesh, {"borderFluid","wall"}),
            _rhs= linear,
            _element= u,
            _expr= zero<FEELPP_DIM,1>()
            );

    Feel::cout << "BALISE NEWTON 5\n";
    bilinear.solve(
            _solution= m_fluidt, 
            _rhs= linear
            );

    Feel::cout << "BALISE NEWTON 6\n";
    double error=normL2(
            _range= elements(m_mesh),
            _expr= idv(u_tmp)-idv(u)
            );
    
    Feel::cout << "BALISE NEWTON 7\n";
    u_tmp=u;
    return error;
}




    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void NavierStokes::run(bilinear_type & bilinear, linear_type & linear, myexpr_type flow, double dt,double maxError)
{

    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto uPrec= m_fluidPrec.element<0>();
    //auto u_tmp= m_fluid.element<0>();


    double error= 0;
    int cpt= 0;
    //u_tmp.on(_range= elements(m_mesh), _expr= zero<FEELPP_DIM,1>());

    linear+=integrate(
            _range= elements(m_mesh),
            _expr= m_rho/dt*inner(idv(uPrec),id(v))
            );

    if(m_modele==modele3)
    {
        //auto autoconv= trans(idv(uPrec))*gradt(u);//Picard
        //auto autoconv_lin= zero<FEELPP_DIM,1>();
        //auto autoconv= (trans(idv(uPrec))*gradt(u)+trans(idt(u))*gradv(uPrec));//newton
        //auto autoconv= zero<FEELPP_DIM,1>()*idt(u);
        //auto autoconv_lin= -trans(idv(uPrec))*gradv(uPrec);


        //bilinear+= integrate(
        //        _range= elements(m_mesh),
        //        _expr= m_rho * autoconv * id(v)
        //        );
        //linear+= integrate(
        //        _range=elements(m_mesh),
        //        _expr= m_rho *autoconv_lin * id(v)
        //        );
    }

    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );

    bilinear+= on(
            _range= markedfaces(m_mesh, {"borderFluid","wall"}),
            _rhs= linear,
            _element= u,
            _expr= zero<FEELPP_DIM,1>()
            );

    bilinear.solve(
            _solution= m_fluidt, 
            _rhs= linear
            );

    //if(m_modele == modele3)
    //{
    //    error= normL2(
    //            _range= elements(m_mesh),
    //            _expr= idv(u)-idv(u_tmp)
    //            );
    //    u_tmp= u;
    //    Feel::cout << "iteration " << cpt << ",\terreur= " << error <<"\n";
    //}
    //cpt++;
    //}while(( maxError<error*0 ) && (cpt<MAX_LOOP_PICARD));


    }

///rename run to run
//    template<typename bilinear_type, typename linear_type, typename myexpr_type>
//void NavierStokes::run(bilinear_type & bilinear, linear_type & linear, myexpr_type flow, double dt)
//{
//
//    //resolve_fluid(bilinear, linear, flow);
//
//    auto v= m_fluid.element<0>();
//    auto u= m_fluidt.element<0>();
//    auto u_tmp= m_fluidPrec.element<0>();
//    auto q= m_fluid.element<1>();
//    auto p= m_fluidt.element<1>();
//
//    if( dt>=0 )
//    {// composante temporelle
//        linear+= integrate(
//                _range= elements(m_mesh), 
//                _expr= m_rho * inner(idv(u_tmp),id(v)) /dt
//                );
//        bilinear+= integrate(
//                _range= elements(m_mesh), 
//                _expr= m_rho * inner(idt(u), id(v))/dt
//                )
//
//            if(m_modele == modele3)
//            {// terme d'auto-convection
//                auto expr_navier_lin= cst(0.);//trans(idv(u_tmp))*gradv(u_tmp)*(dt-1);
//                auto expr_navier_bil= trans(idt(u))*gradv(u_tmp);//(trans(idt(u))*gradv(u_tmp) + trans(idv(u_tmp))*gradt(u))*dt;
//
//                //linear+= integrate(
//                //        _range= elements(m_mesh), 
//                //        _expr= m_rho * expr_navier_lin * id(v)
//                //        );
//                bilinear+= integrate(
//                        _range= elements(m_mesh), 
//                        _expr= m_rho * expr_navier_bil * id(v)
//                        );
//            }
//    }
//
//
//    tic();
//    bilinear+= on(//condition of Dirichlet in the edge 
//            _range= markedfaces(m_mesh, {"borderFluid", "wall"}), 
//            _rhs= linear, 
//            _element= u, 
//            _expr= zero<FEELPP_DIM, 1>()
//            );
//
//    bilinear+= on(//2
//            _range= markedfaces(m_mesh, {"in1", "in2"}), 
//            _rhs= linear, 
//            _element= u, 
//            _expr= flow
//            );
//    toc("  edges  ");
//
//
//    tic();
//    bilinear.solve(
//            _solution= m_fluidt, 
//            _rhs= linear
//            );
//    toc("  solve  ");
//}




    template<typename myexpr_type>
void NavierStokes::run(myexpr_type flow)
{
    auto bilinear= form2(
            _test= Vph,
            _trial= Vph,
            _matrix= m_matrix
            );
    auto linear= form1(
            _test= Vph,
            _vector= m_vector
            );

    run(bilinear, linear, flow);
}








#endif
