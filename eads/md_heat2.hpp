#ifndef __MD_HEAT_HPP__
#define __MD_HEAT_HPP__ 1


/**
 *\file md_projet.hpp
 */

//#include <feel/feel.hpp>
//#include <fstream>
//#include <sstream>
//#include <string>
#include "md_commun.hpp"


using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> M_type;
//typedef Feel::Expr




//// metode utlisation
// Heat heat;
// auto bilinear= form2(_test=heat.Th, _trial= heat.Th,_matrix=heat.matrix);
// auto linear= form1(test=heat.Th, _vector= heat.vector);
// heat.init(bilinear,linear);
// heat.run(...);
//
// ou
//
// Heat heat;
// auto bilinear= form2(_test=heat.Th, _trial= heat.Th,_matrix=heat.matrix);
// auto linear= form1(test=heat.Th, _vector= heat.vector);
// auto bilinear_static= form2(_test=heat.Th, _trial= heat.Th);
// auto linear_static= form1(test=heat.Th);
// heat.init(bilinear_static,linear_static);
// while()
// {
//      bilinear.zero()
//      linear.zero()
//      bilinear+=bilinear_static
//      linear+=linear_static
//      heat.run(...);
// }
// 


class Heat
{
    // typedef m_mesh
    typedef Mesh<M_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // typedef bases
    typedef Lagrange<0, Scalar, Discontinuous> bases_param;
    typedef Lagrange<2, Scalar, Continuous> bases_P2;
    typedef Lagrange<2, Vectorial, Continuous> bases_vect;

    // typedef functionSpace
    typedef boost::parameter::void_ param_void;
    typedef FunctionSpace< mesh_type, bases<bases_param> > space_param_type;
    typedef FunctionSpace< mesh_type, bases<bases_P2>> space_heat_type; 
    typedef FunctionSpace< mesh_type, bases<bases_vect> > space_conv_type;

    // typedef pointer of functionSpace
    typedef boost::shared_ptr<space_param_type> space_param_ptrtype;
    typedef boost::shared_ptr<space_heat_type> space_heat_ptrtype;
    typedef boost::shared_ptr<space_conv_type> space_conv_ptrtype;

    // typedef of functionSpace's element
    typedef typename space_param_type::element_type element_param_type;
    typedef typename space_heat_type::element_type element_heat_type;
    typedef typename space_conv_type::element_type element_conv_type;

    // typedef pointer of element
    typedef boost::shared_ptr<element_param_type> element_param_ptrtype;
    typedef boost::shared_ptr<element_heat_type> element_heat_ptrtype;
    typedef boost::shared_ptr<element_conv_type> element_conv_ptrtype;


    // typedef matrix
    //template<namespace value_type>
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;

    // typedef ptr matrix
    typedef typename backend_type::sparse_matrix_ptrtype matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    private:

    void init_param();
    void init_heat();
    void init_conv();

    public:


    //parametre

    // m_modele du m_modele
    Modele_type m_modele=modele0;


    // m_mesh
    mesh_ptrtype m_mesh;

    // parameter of the differential equation
    space_param_ptrtype P0;
    element_param_type m_rc;
    element_param_type m_k;

    // function of heat
    space_heat_ptrtype Th;
    element_heat_type u;
    element_heat_type ut;
    element_heat_type uPrec;

    // function of heat convection
    space_conv_ptrtype Th_vect;
    element_conv_type beta;

    //backend
    backend_ptrtype m_backend;

    //matrix
    matrix_ptrtype matrix;
    vector_ptrtype vector;




    // private function


    template<typename bilinear_type, typename linear_type, typename myexpr_type>
        void build_heat_stab(bilinear_type & bilinear, linear_type & linear, myexpr_type Q);

    template<typename bilinear_type, typename linear_type>
        void init_matrix(bilinear_type & bilinear,linear_type & linear);

    template<typename expr_type>
        void run(expr_type Q, double dt);

    template<typename myexpr_type>
        void run(myexpr_type Q);


    /// \fn Heat
    /// \brief constructor
    ///
    Heat(){}
    void init(mesh_ptrtype, Modele_type);

    template<typename E>
        void init_beta(E expr);

    double get_heat_IC1();
    double get_heat_IC2();
    double get_heat_out();
    double get_error_Lu();

};






void Heat::init(mesh_ptrtype theMesh, Modele_type mod)
{
    m_modele=mod;
    this->m_mesh= createSubmesh( theMesh, elements(theMesh));
    init_param();
    init_heat();
    
    matrix->zero();
    vector->zero();
}





// function of the class Heat============================================



//! \fn init_param
//! \brief init the function's Space of the parameter 
//! create the parameter m_rc and m_k like space function P0's element
void Heat::init_param()
{
    tic();
    P0 = space_param_type::New(_mesh=m_mesh);

    // coefficient of thermic condictivity
    m_k = P0->element("thermic conductivity");
    m_k.on(_range= markedelements(m_mesh, "PCB"), _expr= cst(doption("PCB.k")));
    m_k.on(_range= markedelements(m_mesh, "IC1"), _expr= cst(doption("Proc.k")));
    m_k.on(_range= markedelements(m_mesh, "IC2"), _expr= cst(doption("Proc.k")));
    m_k.on(_range= markedelements(m_mesh, "AIR"), _expr= cst(doption("Air.k")));

    // coefficient of capacity 
    m_rc = P0->element("capacity");
    m_rc.on(_range=markedelements(m_mesh, "PCB"), _expr= cst(doption("PCB.rc")));
    m_rc.on(_range=markedelements(m_mesh, "IC1"), _expr= cst(doption("Proc.rc")));
    m_rc.on(_range=markedelements(m_mesh, "IC2"), _expr= cst(doption("Proc.rc")));
    m_rc.on(_range=markedelements(m_mesh, "AIR"), _expr= cst(doption("Air.rc"))); 
    toc("init P0   ");
}


//! \fn init_heat
//! \brief init the function's Space for resolution of the Heat equation
//! create the space function Th=\f$ P_2[X]\f$ for the heat equation
//! init :
//! * u     : the test function
//! * ut    : the trial function
//! * uPrec : the memory function( heat of the time t-dt)(init : ambiant heat)
//! * beta  : air flow (init at zero)
void Heat::init_heat()
{
    tic();
    // Space function for the equation of heat
    Th = space_heat_type::New(_mesh= m_mesh);
    ut = Th->element("trial Heat");
    uPrec = Th->element("temporate Heat");
    u  = Th->element("test Heat ");

    Th_vect = space_conv_type::New(_mesh= m_mesh);
    beta = Th_vect->element("Heat convection");
    beta.on(_range=elements(m_mesh), _expr=zero<FEELPP_DIM, 1>());
    uPrec.on(_range=elements(m_mesh), _expr=cst(doption("Modele.Tamb")));

    m_backend= backend(_name="m_backend");
    matrix= m_backend->newMatrix(Th, Th);
    vector= m_backend->newVector(Th);

    toc("Th/Th_vect");
}












//! \fn get_heat_IC1
//! \brief give the average of the heat for the first processor
double Heat::get_heat_IC1()
{
    return mean(
            _range= markedelements(m_mesh, "IC1"), 
            _expr= idv(ut)).value();
}


//! \fn get_heat_IC2
//! \brief give the average of the heat for the second processor
double Heat::get_heat_IC2()
{
    return mean(
            _range= markedelements(m_mesh, "IC2"), 
            _expr= idv(ut)).value();
}


//! \fn get_heat_out
//! \brief give the average of the heat for the out of the air conduct 
double Heat::get_heat_out()
{
    return mean(
            _range= markedfaces(m_mesh, "out"), 
            _expr= idv(ut)).value();
}

//! \fn get_error_Lu()
//! \brief give the result of Lu in norm L2
double Heat::get_error_Lu()
{
    auto L=-idv(m_k)*laplacianv(ut) + idv(m_rc)*gradv(ut)*idv(beta);
    return normL2(
            _range= elements(m_mesh), 
            _expr= L
            );
}














//! \fn build_heat_bilinear_stab
//! \brief use method of Galerkin Least Square for stabilise the system
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Heat::build_heat_stab(bilinear_type & bilinear, linear_type & linear, myexpr_type Q)
{
    tic();
    // forme bilineaire de notre probleme

    if(boption("Time.time"))
    {
        auto L= -idv(m_k)*laplacian(u) + idv(m_rc)*( grad(u)*idv(beta) ) + idv(m_rc)*id(u)/doption("Time.dt");
        auto Lt= -idv(m_k)*laplaciant(ut) +  idv(m_rc) * ( gradt(ut)*idv(beta) ) + idv(m_rc)*idt(ut)/doption("Time.dt"); 
        auto f = Q + idv(m_rc)*idv(uPrec)/doption("Time.dt");

        auto delta= doption("Modele.epsilon") * cst(1.)/( 1./h() + idv(m_k)/(h()*h()) );

        bilinear+=integrate(
                _range= elements(m_mesh), 
                _expr= delta*L*Lt
                );    

        linear+=integrate(
                _range= markedelements(m_mesh, {"IC1", "IC2"}), 
                _expr= delta*L*f
                );
    }
    else
    {

        auto L= -idv(m_k)*laplacian(u) + idv(m_rc)*( grad(u)*idv(beta) );
        auto Lt= -idv(m_k)*laplaciant(ut) +  idv(m_rc) * ( gradt(ut)*idv(beta) );
        auto f=Q;


        auto delta= doption("Modele.epsilon") * cst(1.)/( 1./h() + idv(m_k)/(h()*h()) );

        bilinear+=integrate(
                _range= elements(m_mesh), 
                _expr= delta*L*Lt
                );    

        linear+=integrate(
                _range= markedelements(m_mesh, {"IC1", "IC2"}), 
                _expr= delta*L*f
                );
    }
    toc("stabilisation");
}




//! \fn init_matrix
//! \brief build the matrix matrix for the term with no-time dependance
    template<typename bilinear_type, typename linear_type>
void Heat::init_matrix(bilinear_type & bilinear, linear_type & linear)
{

    linear.zero();
    bilinear.zero();
    bilinear+= integrate(
            _range= elements(m_mesh), 
            _expr= idv(m_k) * inner(grad(u), gradt(ut))
            );// diffusion of heat
}

//template<typename bilinear_type, typename linear_type>
//void update(bilinear_type & bilinear,linear_type & linear)
//{
//    matrix.
//}



// RUN //
    template<typename myexpr_type>
void Heat::run(myexpr_type Q)
{
    tic();
    double dt= doption("Time.time");

    tic();
    auto linear= form1(_test= Th, _matrix= matrix);
    auto bilinear= form2(_test= Th, _trial= Th, _vector=vector);
    toc("init bil/lin");

    toc();

    tic();
    linear+= integrate(
            _range= markedelements(m_mesh, {"IC1", "IC2"}), 
            _expr= Q*id(u)
            );//heating of the processors


    if( m_modele != modele0 )
    {
        bilinear+= integrate(
                _range= markedelements(m_mesh, "AIR"), 
                _expr= idv(m_rc) * (gradt(ut) * idv(beta)) *id(u)
                );//convection of the heat

        if(boption("Modele.GaLS"))
            build_heat_stab(bilinear, linear, Q);
    }

    if(boption("Time.time"))
    {
        linear+= integrate( // terme de memoire/inertie
                _range= elements(m_mesh), 
                _expr= idv(m_rc) * id(u) * idv(uPrec)/dt
                );

        bilinear+= integrate( // terme de reaction
                _range= elements(m_mesh),
                _expr= idv(m_rc) * id(u) * idt(ut)/dt
                );
    }


    bilinear+=on( // temperature du flot d'entree
            _range= markedfaces(m_mesh, {"in1", "in2"}), 
            _rhs= linear, 
            _element= ut, 
            _expr= cst(doption("Modele.Tamb"))
            );
    toc("  build  "); 

    //build_heat_dynamic(bilinear, linear, Q);
    //build_heat_diric_edge(bilinear, linear);

    tic();
    m_backend->solve(
            _solution= ut,
            _matrix=matrix,
            _rhs=vector
            );
    toc("  solve  ");   
    toc("HEAT");
}













//! \fn run_step
    template<typename myexpr_type>
void Heat::run( myexpr_type Q)
{
    auto bilinear= form2( _test= Th, _trial= Th, _matrix= matrix);
    auto linear= form1( _test= Th,_vector=vector);

    run(bilinear, linear, Q);
}




















#endif
