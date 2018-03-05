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
//#include "commun.hpp"

#define MAX_LOOP_NAVIER 10
#define MIN_ERROR_NAVIER 1e-4

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

    //m_backend
    backend_ptrtype m_backend;

    //matrix
    matrix_ptrtype matrix;
    matrix_ptrtype matrix_static;
    vector_ptrtype vector;
    vector_ptrtype vector_static;




    //function


    NavierStokes(){}
    void init(mesh_ptrtype mesh, Modele_type mod);
    void init_fluid();


    void init_matrix();

    void build_time_linear();

    void reset_dynamic();


    template<typename myexpr_type>
        void run(myexpr_type D);

    template<typename myexpr_type>
        void run_stokes(myexpr_type flow);

    template<typename myexpr_type>
        double run_navier(myexpr_type flow);



};





void NavierStokes::init(mesh_ptrtype mesh, Modele_type mod)
{
    tic();
    m_mesh= createSubmesh(mesh, markedelements(mesh, "AIR"));
    m_modele= mod;

    //    init_conv();
    Vph= space_fluid_type::New(_mesh= m_mesh);

    m_fluid= Vph->element(); //fonction servant a stocker la solution du solver 
    m_fluidt= Vph->element();// fonction servant temporairement a stocker la solution
    m_fluidPrec= Vph->element();// fonction servant a stocker l'etat precedent

    // parametre du fluide
    m_mu= doption("Air.mu");
    m_rho= doption("Air.rho");

    m_fluidPrec.element<0>().on(
            _range= elements(m_mesh), 
            _expr= zero<FEELPP_DIM, 1>()
            );

    m_backend= backend( _name= "backend_fluid");
    
    matrix= m_backend->newMatrix(Vph, Vph);
    vector= m_backend->newVector(Vph);

    matrix_static= m_backend->newMatrix(Vph, Vph);
    vector_static= m_backend->newVector(Vph);
    
    // mise en place de la matrice matrix_static
    init_matrix();
    toc("init NAVIERSTOKES");
}


// function of the class NavierStokes============================================ 




/// \fn init_matrix
/// \brief add to the bilinear form the static part
/// \f$ += \int_\Omega \bigl(\frac{\nabla u+\nabla u^T}{2}\bigl):\nabla v\f$
/// \f$ += \int_\Omega q\nabla u+p\nabla v\f$
void NavierStokes::init_matrix()
{
    auto bilinear= form2(_test= Vph,_trial= Vph, _matrix= matrix_static);
    auto linear= form1(_test= Vph, _vector= vector_static);
    bilinear.zero();
    linear.zero();

    auto v= m_fluid.template element<0>();
    auto u= m_fluid.template element<0>();
    auto q= m_fluid.template element<1>();
    auto p= m_fluid.template element<1>();
    auto elmts= elements(m_mesh);
    
    auto deft= sym(gradt(u));
    auto Id= eye<FEELPP_DIM,FEELPP_DIM>();
    auto sigmat=-idt(p)*Id+m_mu*gradt(u);
    //auto sigmat=-idt(p)*Id + 2*m_mu*deft;

    // terme principal
    bilinear+= integrate(
            _range= elmts,
            _expr= inner(sigmat,grad(v))
            );

    // condition d'incompressibilité
    bilinear+= integrate(
            _range= elmts,
            _expr= m_rho * divt(u)*id(q)
            );

    // terme de reaction temporelle
    if(boption("Time.time"))
    {
        double dt=doption("Time.dt");
        bilinear+= integrate(
                _range= elements(m_mesh), 
                _expr= m_rho * inner(idt(u), id(v)) /dt
                );
    }
}


void NavierStokes::reset_dynamic()
{
    matrix->zero();
    matrix->addMatrix(1,matrix_static);
    
    vector->zero();
    vector= vector_static;
}





void NavierStokes::build_time_linear()
{
    auto linear= form1(_test= Vph, _vector= vector);
    auto v= m_fluid.element<0>();
    auto uPrec= m_fluidPrec.element<0>();
    auto dt= doption("Time.dt");

    // terme d'inertie du systeme
    linear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_rho * inner(id(v), idv(uPrec))/dt
            );
}





    template<typename myexpr_type>
void NavierStokes::run_stokes(myexpr_type flow)
{
    auto v= m_fluid.element<0>();
    auto u= m_fluid.element<0>();
    auto u_tmp= m_fluidt.element<0>();

    // re... de la matrice matrix et du vecteur vector
    reset_dynamic();
    auto linear= form1(_test= Vph, _vector= vector);
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix);
    if(boption("Time.time"))
        build_time_linear();


    // condition au bord
    bilinear+= on(// condition a l'entree d'air :dirichlet
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );
    bilinear+= on(// condition sur les parois
            _range= markedfaces(m_mesh, "borderfluid"), 
            _rhs= linear, 
            _element= u, 
            _expr= zero<FEELPP_DIM, 1>()
            );

    // resolution
    m_backend->solve(
            _solution= m_fluid, 
            _matrix= matrix, 
            _rhs= vector,
            _rebuild= true
            );
}




    template<typename myexpr_type>
double NavierStokes::run_navier(myexpr_type flow)
{
    tic();
    auto v= m_fluid.element<0>();
    auto u= m_fluidt.element<0>();// terme temporaire pour faire la difference
    auto u_tmp= m_fluid.element<0>();
    //u.on(_range=elements(m_mesh), _expr= zero<FEELPP_DIM,1>());

    // re... du systeme
    reset_dynamic();
    auto linear= form1(_test= Vph, _vector= vector);
    auto bilinear= form2(_test= Vph, _trial= Vph, _matrix= matrix);
    if(boption("Time.time"))
        build_time_linear();


    // terme d'auto-convection
    //auto autoconv= inner(gradt(u)*idv(u_tmp)+gradv(u_tmp)*idt(u),id(v)); 
    auto autoconv= inner(gradt(u)*idv(u_tmp),id(v)); 
    //auto autoconv_lin= inner(gradv(u_tmp)*idv(u_tmp),id(v)); 

    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= m_rho * autoconv//m_rho * autoconv * id(v)
            );
    //linear+= integrate(
    //        _range= elements(m_mesh), 
    //        _expr= m_rho * autoconv_lin//m_rho * autoconv * id(v)
    //        );


    //condition au bord
    bilinear+= on(//2
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= u, 
            _expr= flow
            );
    bilinear+= on(
            _range= markedfaces(m_mesh, "borderfluid"), 
            _rhs= linear, 
            _element= u, 
            _expr= zero<FEELPP_DIM, 1>()
            );

    // resolution
    m_backend->solve(
            _solution= m_fluidt, 
            _matrix= matrix, 
            _rhs= vector
            );

    //evaluation de l'erreur
    double error= normL2(
            _range= elements(m_mesh), 
            _expr= idv(u_tmp)-idv(u)
            );
    // remplacement de la variabla temporaire pa la solution
    u_tmp= u;
    toc("run_navier");
    return error;
}




    template<typename myexpr_type>
void NavierStokes::run(myexpr_type flow)
{
    tic();


    if(m_modele== modele2)
        run_stokes(flow);
    else
    {
        init_matrix();
        bool est_non_fini= true;
        double error;
        m_fluidt=m_fluidPrec;
        
        for(int i= 0;(i<MAX_LOOP_NAVIER) && est_non_fini;i++)
        {
            error=run_navier(flow);
            if(error<MIN_ERROR_NAVIER)
                est_non_fini=false;
            Feel::cout << i << ": erreur= " << error << "\n";
        }
    }

    auto v_tmp= m_fluid.element<0>();
    auto p_tmp= m_fluid.element<1>();
    //auto v_norm= (
    //        v_tmp.comp<ComponentType::X>().pow(2)+
    //        v_tmp.comp<ComponentType::Y>().pow(2)
    //        ).sqrt();

    //Feel::cout << "marqueur 1\n";
    //double v_mean=  normL2(
    //        _range= elements(m_mesh),
    //        _expr= idv(v_norm)
    //        ); 

    //Feel::cout << "marqueur 2\n";
    double p_mean= mean(
            _range=elements(m_mesh),
            _expr= idv(p_tmp)
            )(0,0);

    //Feel::cout << "marqueur 3\n";
    p_tmp.on(
            _range= elements(m_mesh),
            _expr= idv(p_tmp)-p_mean);

    //Feel::cout << "la vitesse moyenne est \n"<< v_mean << "\n";

    m_fluidPrec.zero();
    m_fluidPrec= m_fluid;
    toc("run Fluid");
}










#endif
