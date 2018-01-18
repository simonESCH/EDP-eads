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

#define MAX_LOOP_PICARD 10000

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
    element_fluid_type m_fluid_tmp;

    //backend
    backend_ptrtype backend_fluid;

    //matrix
    matrix_ptrtype matrix_fluid;
    vector_ptrtype vector_fluid;




    //function

#if TEST_PROJET_HPP
    NavierStokes(){}
    void init(mesh_ptrtype mesh, Modele_type mod);
#else
    NavierStokes(mesh_ptrtype mesh, Modele_type mod);
#endif
    void init_fluid();



    void build_fluid_static();

    template<typename L, typename V, typename P>
        void build_fluid_dynamic_linear(L & linear, V velocity, V velocity_previous, double dt);

    template<typename B, typename V>
        void build_fluid_dynamic_bilinear(B & bilinear, V velocity_test, V velocity_trial, V velocity_loop);

    template<typename myexpr_type>
        void run(myexpr_type D);

    template<typename myexpr_type>
        void run_step(myexpr_type D,double dt);

    template<typename B, typename L, typename E>
        void resolve_fluid(B & bilinear, L & linear, E & flow, double error= 1e-3);

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
    Vph = space_fluid_type::New(_mesh= m_mesh);

    m_fluid = Vph->element();
    m_fluidt = Vph->element();
    m_fluid_tmp = Vph->element();


    m_fluidt.element<0>().on(
            _range= elements(m_mesh),
            _expr= zero<FEELPP_DIM,1>()
            );

    backend_fluid= backend(_name= "backend_fluid");
    matrix_fluid= backend_fluid->newMatrix(Vph,Vph);
    vector_fluid= backend_fluid->newVector(Vph);

}


#else


NavierStokes::NavierStokes(mesh_ptrtype i_mesh,Modele_type mod)
{
    m_mesh= createSubmesh(i_mesh, markedelements(i_mesh, "AIR"));
    m_modele= mod;

    //    init_conv();
    Vph = space_fluid_type::New(_mesh= m_mesh);

    m_fluid = Vph->element();
    m_fluidt = Vph->element();
    m_fluid_tmp = Vph->element();


    m_fluidt.element<0>().on(
            _range= elements(m_mesh),
            _expr= zero<FEELPP_DIM,1>()
            );

    backend_fluid= backend(_name= "backend_fluid");
    matrix_fluid= backend_fluid->newMatrix(Vph,Vph);
    vector_fluid= backend_fluid->newVector(Vph);

}
#endif


// function of the class NavierStokes============================================ 




/// \fn build_fluid_static
/// \brief add to the bilinear form the static part
/// \f$ += \int_\Omega \bigl(\frac{\nabla u+\nabla u^T}{2}\bigl):\nabla v\f$
/// \f$ += \int_\Omega q\nabla u+p\nabla v\f$
void NavierStokes::build_fluid_static()
{
    tic();
    matrix_fluid->zero();
    vector_fluid->zero();
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix_fluid);

    auto v= m_fluid.template element<0>();
    auto u= m_fluidt.template element<0>();
    auto q= m_fluid.template element<1>();
    auto p= m_fluidt.template element<1>();

    auto mat_ID= eye<FEELPP_DIM, FEELPP_DIM>();
    double m_mu= doption("Air.m_mu");
    /// \f$ \bigl(\nabla u -pI\!d\bigl)v\cdot n\f$ is the edge's term
    //auto deft= m_mu*sym(gradt(u))-idt(p)*mat_ID;
    //auto deft= m_mu*gradt(vt)-idt(pt)*mat_ID;

    //auto stokes_expr= inner( 2* m_mu* deft - idt(p)*mat_ID, grad(v))+
    //                    id(q)* divt(u);

    //Feel::cout<<"test "<<stokes_expr<<"\n";

    bilinear+= integrate(
            _range= elements(m_mesh),
            _expr= m_mu*inner(grad(v), sym(gradt(u)))
            +id(q)*divt(u)-idt(p)*div(v)
            );

    toc("NS static  ");
    //bilinear+= integrate(
    //        _range= elements(m_mesh),
    //        _expr= id(p)*divt(vt)-idt(pt)*div(v)
    //        );
}










    template<typename linear_type, typename velocity_type, typename pressure_type>
void NavierStokes::build_fluid_dynamic_linear(linear_type & linear, velocity_type v,velocity_type v_prec, double dt)
{

}

    template<typename bilinear_type, typename velocity_type>
void NavierStokes::build_fluid_dynamic_bilinear(bilinear_type & bilinear, velocity_type v, velocity_type vt, velocity_type vl)
{

}


/// \brief resolve the navier-stokes with de method of Picard 
    template<typename bilinear_type, typename linear_type, typename expr_type>
void NavierStokes::resolve_fluid(bilinear_type & bilinear, linear_type & linear, expr_type & flow,double error)
{
    auto bilinear_static= form2(_test= Vph,_trial= Vph,_matrix= matrix_fluid);
    auto linear_static= form1(_test= Vph,_vector= vector_fluid);

    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto u_tmp= m_fluidt.element<0>();

    double loop_err;
    double cpt= 0;

    do
    {
        m_fluid_tmp= m_fluidt;
        bilinear= bilinear;
        linear= linear;

        // partie d'auto-convection u\nabla u
        auto expr_navier= ( trans(idv(u_tmp))*gradt(u) )*id(v);
        bilinear+= integrate(
                _range= elements(m_mesh),
                _expr= m_rho*expr_navier
                );

        bilinear+= on(//2
                _range= markedfaces(m_mesh, {"in1", "in2"}),
                _rhs= linear,
                _element= u,
                _expr= flow
                );

        bilinear.solve(
                _solution= m_fluidt,
                _rhs= linear
                );

        loop_err= normL2(
                _range= elements(m_mesh),
                _expr= (idv(u)-idv(u_tmp))
                );

        cpt++;
        // a completer
    } while((error<loop_err)&&(cpt<MAX_LOOP_PICARD));
}

















///rename run to run
    template<typename myexpr_type>
void NavierStokes::run(myexpr_type flow)
{
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();
    auto q= m_fluid.element<1>();
    auto p= m_fluidt.element<1>();


    build_fluid_static();

    //matrix_fluid->zero();
    //vector_fluid->zero();

    auto bilinear= form2( _test= Vph, _trial= Vph, _matrix= matrix_fluid);
    auto linear= form1( _test= Vph, _vector= vector_fluid);

    tic();
    bilinear+= on(//condition of Dirichlet in the edge 
            _range= markedfaces(m_mesh, {"borderFluid","wall"}),
            _rhs= linear,
            _element= u,
            _expr= zero<FEELPP_DIM,1>()
            );

    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}),
            _rhs= linear,
            _element= u,
            _expr= flow
            );
    toc("NS edge     ");

    tic();
    bilinear.solveb(
            _solution= m_fluidt, 
            _rhs= linear,
            _backend= backend_fluid);
    toc("NS solve   ");



}

    template<typename myexpr_type>
void run_step(myexpr_type D,double dt)
{
    // vide pour l'instant
}




#endif
