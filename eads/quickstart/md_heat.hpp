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


    /// \enum modele_enum
    /// \brief modele entre 4 modele different
    /// * modele0 : pas de convection du a l'air
    /// * modele1 : on suppose que le flux d'air suit un profile de poiseuille 
    /// * modele2 : on suppose que l'air suit les equations de modele2
    /// * navier modele2 : on suppose que l'air suit les equation de navier-modele2

    public:


    //parametre

    // modele du modele
    Modele_type modele=modele0;


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
    element_heat_type u_tmp;

    // function of heat convection
    space_conv_ptrtype Th_vect;
    element_conv_type beta;

    //backend
    backend_ptrtype backend_heat;

    //matrix
    matrix_ptrtype matrix_heat;
    vector_ptrtype vector_heat;




    // private function

    void init_param();
    void init_heat();
    void init_conv();




    template<typename bilinear_type, typename linear_type>
        void build_heat_diric_edge(bilinear_type & bilinear, linear_type & linear);

    template<typename bilinear_type, typename linear_type, typename myexpr_type>
        void build_heat_stab(bilinear_type & bilinear, linear_type & linear, myexpr_type Q);

    void build_heat_static(double dt=-1);

    template<typename bilinear_type, typename linear_type,typename myexpr_type>
        void build_heat_dynamic(bilinear_type & bilinear, linear_type & linear, myexpr_type Q, double dt=-1);

    template<typename myexpr_type>
        void run(myexpr_type Q);

    template<typename myexpr_type>
        void run_step( myexpr_type Q, double dt=-1);


    /// \fn Heat
    /// \brief constructor
    ///
    Heat(mesh_ptrtype mesh, Modele_type mod);

    template<typename E>
        void init_beta(E expr);

    double get_heat_IC1();
    double get_heat_IC2();
    double get_heat_out();
    double get_error_Lu();

};







Heat::Heat(mesh_ptrtype theMesh,Modele_type mod):modele(mod)
{
    this->m_mesh= createSubmesh( theMesh, elements(theMesh));
    init_param();
    init_heat();
}


// function of the class Heat============================================



//! \fn init_param
//! \brief init the function's Space of the parameter 
//! create the parameter m_rc and m_k like space function P0's element
void Heat::init_param()
{
    tic();
    P0 = space_param_type::New(_mesh=m_mesh,_extended_doftable=true);

    // coefficient of thermic condictivity
    m_k = P0->element("thermic conductivity");
    m_k.on(_range= markedelements(m_mesh, "PCB"),_expr= cst(doption("PCB.k")));
    m_k.on(_range= markedelements(m_mesh, "IC1"),_expr= cst(doption("Proc.k")));
    m_k.on(_range= markedelements(m_mesh, "IC2"),_expr= cst(doption("Proc.k")));
    m_k.on(_range= markedelements(m_mesh, "AIR"),_expr= cst(doption("Air.k")));

    // coefficient of capacity 
    m_rc = P0->element("capacity");
    m_rc.on(_range=markedelements(m_mesh, "PCB"),_expr= cst(doption("PCB.rc")));
    m_rc.on(_range=markedelements(m_mesh, "IC1"),_expr= cst(doption("Proc.rc")));
    m_rc.on(_range=markedelements(m_mesh, "IC2"),_expr= cst(doption("Proc.rc")));
    m_rc.on(_range=markedelements(m_mesh, "AIR"),_expr= cst(doption("Air.rc"))); 
    toc("init P0   ");
}


//! \fn init_heat
//! \brief init the function's Space for resolution of the Heat equation
//! create the space function Th=\f$ P_2[X]\f$ for the heat equation
//! init :
//! * u : the test function
//! * ut: the trial function
//! * u_tmp :the memory function( heat of the time t-dt)
void Heat::init_heat()
{
    tic();
    // Space function for the equation of heat
    Th = space_heat_type::New(_mesh= m_mesh);
    ut = Th->element("trial Heat");
    u  = Th->element("test Heat ");

    Th_vect = space_conv_type::New(_mesh= m_mesh);
    beta = Th_vect->element("Heat convection");
    beta.on(_range=elements(m_mesh),_expr=zero<FEELPP_DIM,1>());

    backend_heat= backend(_name="backend_heat");
    matrix_heat= backend_heat->newMatrix(Th,Th);
    vector_heat= backend_heat->newVector(Th);

    if(boption("Time.time"))// A CHANGER)
        build_heat_static(doption("Time.dt"));
    else
        build_heat_static();

    toc("Th/Th_vect");
}












//! \fn get_heat_IC1
//! \brief give the average of the heat for the first processor
double Heat::get_heat_IC1()
{
    return mean(
            _range= markedelements(m_mesh,"IC1"),
            _expr= idv(ut)).value();
}


//! \fn get_heat_IC2
//! \brief give the average of the heat for the second processor
double Heat::get_heat_IC2()
{
    return mean(
            _range= markedelements(m_mesh,"IC2"),
            _expr= idv(ut)).value();
}


//! \fn get_heat_out
//! \brief give the average of the heat for the out of the air conduct 
double Heat::get_heat_out()
{
    return mean(
            _range= markedfaces(m_mesh,"out"),
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















//! \param bilinear
//! \param linear
    template<typename bilinear_type, typename linear_type>
void Heat::build_heat_diric_edge(bilinear_type & bilinear, linear_type & linear)
{
    bilinear+=on(
            _range= markedfaces(m_mesh, {"in1", "in2"}),
            _rhs= linear,
            _element= ut,
            _expr= cst(doption("Modele.Tamb"))
            );
}




//! \fn build_heat_bilinear_stab
//! \brief use method of Galerkin Least Square for stabilise the system
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Heat::build_heat_stab(bilinear_type & bilinear, linear_type & linear, myexpr_type Q)
{
    auto L= -idv(m_k)*laplacian(u) + idv(m_rc)*( grad(u)*idv(beta) );
    auto Lt= -idv(m_k)*laplaciant(ut) +  idv(m_rc) * ( gradt(ut)*idv(beta) );
    auto delta= doption("Modele.epsilon") * cst(1.)/( 1./h() + idv(m_k)/(h()*h()) );
    
    bilinear+=integrate(
            _range= elements(m_mesh),
            _expr= delta*L*Lt
            );    

    linear+=integrate(
            _range= markedelements(m_mesh,{"IC1", "IC2"}),
            _expr= L*delta*Q
            );
}

//! \fn build_heat_static
//! \brief build the matrix matrix_heat for the term with no-time dependance
void Heat::build_heat_static(double dt)
{
    matrix_heat->zero();
    vector_heat->zero();

    auto bilinear= form2(_test= Th, _trial= Th, _matrix= matrix_heat);
    bilinear+= integrate(
            _range= elements(m_mesh),
            _expr= idv(m_k) * inner(grad(u), gradt(ut))
            );// diffusion of heat
    if( dt>0 )
    {
        bilinear+= integrate(
                _range= elements(m_mesh),
                _expr= idv(m_rc) / dt * id(u) * idt(ut)
                ); // term on time
    }
}



//! \fn build_heat_dynamic
//! \brief add the part dependant ot time
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Heat::build_heat_dynamic(bilinear_type & bilinear, linear_type & linear, myexpr_type Q, double dt)
{
    linear+= integrate(
            _range= markedelements(m_mesh,{"IC1", "IC2"}),
            _expr= Q*id(u)
            );//heating of the processors

    if( modele != modele0 )
    {
        bilinear+= integrate(
                _range= markedelements(m_mesh,"AIR"),
                _expr= idv(m_rc) * (gradt(ut) * idv(beta)) *id(u)
                );//convection of the heat
        if(boption("Modele.GaLS"))
            build_heat_stab(bilinear, linear, Q);
    }
    if(dt>0)
    linear+= integrate(
            _range= elements(m_mesh),
            _expr= idv(m_rc) /dt * id(u) * idv(u_tmp)
            );// term of time's memory 

}









// RUN //


    template<typename myexpr_type>
void Heat::run(myexpr_type Q)
{
    Feel::cout<<"run\n";

    tic();

    auto bilinear= form2(_test= Th, _trial= Th, _matrix= matrix_heat);
    auto linear= form1(_test= Th,  _vector= vector_heat);

    build_heat_dynamic(bilinear, linear, Q);
    build_heat_diric_edge(bilinear, linear);
    toc("build heat");

    tic();
    bilinear.solveb(
            _solution= ut,
            _rhs=linear,
            _backend= backend_heat
            );
    toc("solve heat");
}


//! \fn run_step
    template<typename myexpr_type>
void Heat::run_step( myexpr_type Q, double dt)
{
    Feel::cout<<"step\n";

    tic();
    auto bilinear_static= form2( _test= Th, _trial= Th, _matrix= matrix_heat);
    auto bilinear= form2( _test= Th, _trial= Th);
    auto linear= form1( _test= Th);

    bilinear= bilinear_static;
    build_heat_dynamic(bilinear, linear, dt, Q);
    build_heat_diric_edge(bilinear, linear);

    toc("build heat");


    tic();
    bilinear.solve(
            _solution= ut,
            _rhs=linear
            );
    toc("solve heat");

}

#if 0
//! \fn run
    template<typename heat_type,typename conv_type>
void Heat::run(heat_type Q, conv_type beta_expr, double T, std::string export_name)
{
    tic();
    auto envoi= exporter(_mesh= m_mesh, _name= export_name);
    double dt= doption("Time.dt");


    tic();
    element_heat_type u_tmp= Th->element();
    //auto f= form2( _test= Th, _trial= Th);
    auto f= form2( _test= Th, _trial= Th, _matrix= matrix_heat);
    build_heat_static( f, dt);
    u_tmp.on( 
            _range= elements(m_mesh),
            _expr = cst(doption("Modele.Tamb"))
            );//la temperature initiale est la temperature ambiante
    toc("init dynamic");



    Feel::cout<<"les resultats seront dans le dossier '"<<export_name<<"'\n";

    for(double t=0;t<T;t+=dt)
    {
        tic();

        // init for the time t the value for the edge and the heat quantity Q
        beta_expr.setParameterValues({{"t",t}});
        Q.setParameterValues({{"t",t}});
        init_beta(beta_expr);

        run_step(f, u_tmp, Q, dt);

        envoi->step(t)->add("heat",u_tmp);
        u_tmp= ut;
        std::ostringstream ostr;
        ostr << "t= " << t << " s";

        toc(ostr.str());

        envoi->save();
        Feel::cout<<"la temperature du processeur IC2 est : "<<get_heat_IC2()<<" K\n"
            <<"la norme de u a qui on a applique L est :"<<get_error_Lu()<<"\n"
            <<"\n";
    }
}
#endif






























#endif

