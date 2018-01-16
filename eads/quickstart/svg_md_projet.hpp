#ifndef __MD_PROJET_HPP__
#define __MD_PROJET_HPP__


/**
 *\file md_projet.hpp
 */

#include <feel/feel.hpp>
#include <fstream>
#include <sstream>
#include <string>

using namespace Feel;
using namespace vf;

typedef Simplex<FEELPP_DIM> M_type;
//typedef Feel::Expr





/// \fn makeOptions
/// \brief declaration of the option 
inline
    po::options_description
makeOptions()
{
    po::options_description myappOptions( "My app options" );
    myappOptions.add_options()
        ("Air.rc", po::value< double>()->default_value( 0 ),
         "Rho*Capa of air")
        ("Proc.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of processor")
        ("PCB.rc", po::value< double>()->default_value( 0 ), 
         "Rho*Capa of motherboard")

        ("Air.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of air")
        ("Proc.k", po::value< double>()->default_value( 1), 
         "conductivity thermic of processor")
        ("PCB.k", po::value< double>()->default_value( 1), 
         "conductivity thermic  of motherboard")

        ("Air.mu",po::value<double>()->default_value( 1),
        "viscosity of Air")
        ("Air.rho",po::value<double>()->default_value( 0),
        "density of Air")

        ("Modele.Tamb", po::value< double>()->default_value(293.15), 
         "temperature of reference")
        ("Proc.Q", po::value< std::string>()->default_value("0."), 
         "Quantity of heat of the processor")
        ("Modele.flux",po::value< std::string>()->default_value("0."), 
         "profile du flux d'air a l'entree")
        ("Modele.modele",po::value< std::string>()->default_value("modele0"),
         "modele du Modele")
        ("Air.D",po::value<std::string>()->default_value("1e3"), 
         "power of the air flow")
        ("Modele.epsilon",po::value<double>()->default_value(1), 
         "epsilon de GaLS")

        ("Modele.Tmax",po::value< double>()->default_value(0), 
         "temps maximal")
        ("Modele.dt",po::value< double>()->default_value(1), 
         "pas de temps")
        ("Modele.GaLS",po::value<bool>()->default_value( false ), 
         "activation de la stabilisation GaLS")

        //("Exporter.save",po::value< std::string>()->default_value(""),
        // "chemin de la sauvegarde")
        //("Exporter.load",po::value< std::string>()->default_value(""),
        // "chemin de l'approximation fine")
        ;
    myappOptions.add(backend_options("backend_heat"));
    myappOptions.add(backend_options("backend_fluid"));
    return myappOptions.add( feel_options() );
}





/// \fn affichage
/// \brief affichage of the option in the file of config
void affichage()
{

    Feel::cout<< "\n\t========================"
        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/mm^3"
        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/mm^3"
        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/mm^3"
        << "\n\t|======================="
        << "\n\t|Air.k  : " << doption("Air.k") << " W/K/mm"
        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/mm"
        << "\n\t|PCB.k  : " << doption("PCB.k") << " W/K/mm"
        << "\n\t|======================="
        << "\n\t|Tamb   : " << doption("Modele.Tamb") << " K"
        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/mm^3"
        << "\n\t|ventil : " << soption("Modele.flux")
        << "\n\t|modele  : " << soption("Modele.modele");

    //if(soption("Exporter.load").compare(""))
    //    Feel::cout<<"\n\t|load   : "<<soption("Exporter.load");
    //if(soption("Exporter.save").compare(""))
    //    Feel::cout<<"\n\t|save   : "<<soption("Exporter.save");

    Feel::cout<< "\n\t|======================="
        << "\n\t|hsize  : " << doption("gmsh.hsize")<<" mm"
        << "\n\t|Tmax   : " << doption("Modele.Tmax")<<" sec"
        << "\n\t|dt     : " << doption("Modele.dt")<<" sec"        
        << "\n\t========================\n\n";
    Feel::cout<<"\n\t|chemin .geo : "<<soption("gmsh.filename")<<"\n";
}






#if 1
class Projet
{
    // typedef mesh
    typedef Mesh<M_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    // typedef bases
    typedef Lagrange<0, Scalar, Discontinuous> bases_param;
    typedef Lagrange<1, Scalar> bases_pressure;
    typedef Lagrange<2, Scalar, Continuous> bases_P2;
    typedef Lagrange<2, Vectorial, Continuous> bases_vect;
    typedef bases<bases_vect, bases_pressure> bases_fluid;

    // typedef functionSpace
    typedef boost::parameter::void_ param_void;
    typedef FunctionSpace< mesh_type, bases<bases_param> > space_param_type;
    typedef FunctionSpace< mesh_type, bases<bases_P2>> space_heat_type; 
    //auto th = space_heat_type::New(_mesh=mesh);
    //typedef Pch_type<mesh_type,2> space_heat_type;
    typedef FunctionSpace< mesh_type, bases<bases_vect> > space_conv_type;
    typedef FunctionSpace< mesh_type, bases_fluid> space_fluid_type;

    // typedef pointer of functionSpace
    typedef boost::shared_ptr<space_param_type> space_param_ptrtype;
    typedef boost::shared_ptr<space_heat_type> space_heat_ptrtype;
    typedef boost::shared_ptr<space_conv_type> space_conv_ptrtype;
    typedef boost::shared_ptr<space_fluid_type> space_fluid_ptrtype;

    // typedef of functionSpace's element
    typedef typename space_param_type::element_type element_param_type;
    typedef typename space_heat_type::element_type element_heat_type;
    typedef typename space_conv_type::element_type element_conv_type;
    typedef typename space_fluid_type::element_type element_fluid_type;

    // typedef pointer of element
    typedef boost::shared_ptr<element_param_type> element_param_ptrtype;
    typedef boost::shared_ptr<element_heat_type> element_heat_ptrtype;
    typedef boost::shared_ptr<element_conv_type> element_conv_ptrtype;
    typedef boost::shared_ptr<element_fluid_type> element_fluid_ptrtype;


    // typedef exporter
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

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
    /// * no_areation : pas de convection du a l'air
    /// * no_stokes : on suppose que le flux d'air suit un profile de poiseuille 
    /// * stokes : on suppose que l'air suit les equations de stokes
    /// * navier stokes : on suppose que l'air suit les equation de navier-stokes
    enum modele_enum{no_areation=0, no_stokes, stokes, navier_stokes};




    //parametre

    // modele du modele
    modele_enum modele=no_areation;


    // mesh
    mesh_ptrtype mesh;
    mesh_ptrtype mesh_TTL;
    mesh_ptrtype mesh_AIR;

    // parameter of the differential equation
    space_param_ptrtype P0;
    element_param_type rc;
    element_param_type k;

    // function of heat
    space_heat_ptrtype Th;
    element_heat_type u;
    element_heat_type ut;
    //element_heat_type u_tmp;

    // function of heat convection
    space_conv_ptrtype Th_vect;
    element_conv_type beta;

    // function of navier stokes
    space_fluid_ptrtype Vph;
    element_fluid_type fluid;
    element_fluid_type fluidt;
    element_fluid_type fluid_tmp;

    //backend
    backend_ptrtype backend_heat;
    backend_ptrtype backend_fluid;

    //matrix
    matrix_ptrtype matrix_heat;
    matrix_ptrtype matrix_fluid;
    vector_ptrtype vector_heat;
    vector_ptrtype vector_fluid;

    // export
    export_ptrtype exp;



    // private function

    void init_mesh();
    void init_param();
    void init_heat();
    void init_conv();
    void init_fluid();
    void init_choice();


    public:

    template<typename L,typename E>
        void build_heat_heat_proc( L & linear, E heat);

    template<typename B> // 4 ligne donc a fusionner avec un autre
        void build_heat_conv(B & bileaire, element_conv_type const & convection);

    template<typename B,typename L>
        void build_heat_diric_edge(B & bilinear, L & linear);

    template<typename L>
        void build_heat_time_linear(L & linear,element_heat_type & u_prec,double dt);

    template<typename B,typename L, typename E>
        void build_heat_stab(B & bilinear, L & linear, element_conv_type conv, E heat);

    template<typename B>
        void build_fluid_static(B & bilinear, double dt=-1);
    
    template<typename L, typename V, typename P>
        void build_fluid_dynamic_linear(L & linear, V velocity, V velocity_previous, double dt);
    
    template<typename B, typename V>
        void build_fluid_dynamic_bilinear(B & bilinear, V velocity_test, V velocity_trial, V velocity_loop);


    template<typename B, typename L, typename E>
        void build_heat_dynamic(B & bilinear, L & linear,element_conv_type beta, E heat);

    template<typename B, typename L, typename E>
        void build_heat_dynamic(B & bilinear, L & linear, element_heat_type u_prec, double dt, element_conv_type conv, E heat);

    template<typename B>
        void build_heat_static(B & bilinear, double dt=-1);

    template<typename bilinear_type, typename myexpr_type>
        void run_heat_step(bilinear_type & bilinear_static, element_heat_type u_prec, myexpr_type Q, double dt);
    
    //element_fluid_type init_in();

    template<typename myexpr_type>
    void run_fluid(myexpr_type D);

    
    std::string init_edge_in(std::string vitesse_flux);

    /// \fn Projet
    /// \brief constructor
    ///
    Projet()
    {
        init_mesh();
        init_choice();
        init_param();
        init_heat();
        init_conv();
        init_fluid();
    }

    template<typename E>
        void init_beta(E expr);

    double get_heat_IC1();
    double get_heat_IC2();
    double get_heat_out();
    double get_error_Lu();

    template<typename myexpr_type>
        void run_heat(myexpr_type Q,std::string="staticSolve");

    template<typename H, typename B>
        void run_heat(H heat, B beta, double T, std::string="dynamicSolve");
};










// function of the class Projet============================================



/// \fn init_mesh()
/// \brief init the mesh of the system
/// create a main mesh( mesh)
/// and 2 submesh : mesh_TTL for the heat equation and mesh_AIR for fluid equation
void Projet::init_mesh()
{
    tic();// construction des maillages
    mesh = loadMesh(_mesh = new Mesh<M_type>);
    mesh_AIR = createSubmesh(mesh, markedelements(mesh,"AIR"));
    mesh_TTL = createSubmesh(mesh, elements(mesh));
    toc("init mesh");
}


/// \fn init_choice
/// \brief set the using model
void Projet::init_choice()
{
    std::string choice=soption("Modele.modele");
    if(!choice.compare("modele1"))
        modele= no_stokes;
    else if(!choice.compare("modele2"))
        modele= stokes;
    else if(!choice.compare("modele3"))
        modele= navier_stokes;
    else if(!choice.compare("modele0"))
        modele= no_areation;
    else
    {
        Feel::cout<<"attention : le modele choisi n'est pas repertorie\n\tle modele choisi sera par defaut\n";
        modele= no_areation;
    }
}


/// \fn init_param
/// \brief init the function's Space of the parameter 
/// create the parameter rc and k like space function P0's element
void Projet::init_param()
{
    tic();
    P0 = space_param_type::New(_mesh=mesh_TTL,_extended_doftable=true);

    // coefficient of thermic condictivity
    k = P0->element("thermic conductivity");
    k.on(_range= markedelements(mesh_TTL, "PCB"),_expr= cst(doption("PCB.k")));
    k.on(_range= markedelements(mesh_TTL, "IC1"),_expr= cst(doption("Proc.k")));
    k.on(_range= markedelements(mesh_TTL, "IC2"),_expr= cst(doption("Proc.k")));
    k.on(_range= markedelements(mesh_TTL, "AIR"),_expr= cst(doption("Air.k")));

    // coefficient of capacity 
    rc = P0->element("capacity");
    rc.on(_range=markedelements(mesh_TTL, "PCB"),_expr= cst(doption("PCB.rc")));
    rc.on(_range=markedelements(mesh_TTL, "IC1"),_expr= cst(doption("Proc.rc")));
    rc.on(_range=markedelements(mesh_TTL, "IC2"),_expr= cst(doption("Proc.rc")));
    rc.on(_range=markedelements(mesh_TTL, "AIR"),_expr= cst(doption("Air.rc"))); 
    toc("init P0  ");
}


/// \fn init_heat
/// \brief init the function's Space for resolution of the Heat equation
/// create the space function Th=\f$ P_2[X]\f$ for the heat equation
/// init :
/// * u : the test function
/// * ut: the trial function
/// * u_tmp :the memory function( heat of the time t-dt)
void Projet::init_heat()
{
    tic();
    // Space function for the equation of heat
    Th = space_heat_type::New(_mesh=mesh_TTL);
    ut = Th->element("trial Heat");
    u  = Th->element("test Heat ");

    backend_heat=backend(_name="backend_heat");
    matrix_heat=backend_heat->newMatrix(Th,Th);
    vector_heat=backend_heat->newVector(Th);
    toc("Th       ");
}


/// \fn init_conv
/// \brief init the vector of heat's convection
/// init the version vectorial of the heat
void Projet::init_conv()
{
    tic();
    // Space function for the convection of the Heat 
    Th_vect = space_conv_type::New(_mesh= mesh_TTL);
    beta = Th_vect->element("Heat convection");
    toc("Th_vect  ");
}


/// \fn init_fluid
/// \brief init of function's Space for the resolution of flow equation
/// create the space function 
void Projet::init_fluid()
{
    tic();
    // Space funtion for Navier-Stockes
    Vph = space_fluid_type::New(_mesh= mesh_AIR);
    fluid = Vph->element();

    // backend
    backend_fluid=backend(_name="backend_fluid");
    matrix_fluid=backend_fluid->newMatrix(Vph,Vph);
    vector_fluid=backend_fluid->newVector(Vph);
    toc("Vph      ");
}










/// \fn get_heat_IC1
/// \brief give the average of the heat for the first processor
double Projet::get_heat_IC1()
{
    return mean(
            _range= markedelements(mesh_TTL,"IC1"),
            _expr= idv(ut)).value();
}


/// \fn get_heat_IC2
/// \brief give the average of the heat for the second processor
double Projet::get_heat_IC2()
{
    return mean(
            _range= markedelements(mesh_TTL,"IC2"),
            _expr= idv(ut)).value();
}


/// \fn get_heat_out
/// \brief give the average of the heat for the out of the air conduct 
double Projet::get_heat_out()
{
    return mean(
            _range= markedfaces(mesh_TTL,"out"),
            _expr= idv(ut)).value();
}

/// \fn get_error_Lu()
/// \brief give the result of Lu in norm L2
double Projet::get_error_Lu()
{
    auto L=-idv(k)*laplacianv(ut) + idv(rc)*gradv(ut)*idv(beta);
    return normL2(
            _range= elements(mesh_TTL),
            _expr= L
            );
}













/// \fn build_heat_heat_proc
/// \brief add to the linear form the heat of processor : \f$Q\f$
/// \param linear
/// \param Q
//    template<typename linear_type, typename myexpr_type>
//void Projet::build_heat_heat_proc(linear_type & linear, myexpr_type Q)
//{
//    linear= integrate(
//            _range= markedelements(mesh_TTL,{"IC1", "IC2"}),
//            _expr= Q*id(u)
//            );
//}

/// \fn build_heat_conv
/// \brief add to the bilinear form the convection : \f$\rho C_{air}(\beta\cdot \nabla u)\f$
/// \param bilinear
/// \param conv 
    template<typename bilinear_type>
void Projet::build_heat_conv(bilinear_type & bilinear, element_conv_type const & conv)
{
    bilinear+= integrate(
            _range= markedelements(mesh_TTL,"AIR"),
            _expr= idv(rc) * (gradt(ut) * idv(conv)) *id(u)
            );
}

/// \fn build_heat_diric_edge
/// \brief add to the linear form the edge condition of dirichlet in the edge "in"
/// \param bilinear
/// \param linear
    template<typename bilinear_type, typename linear_type>
void Projet::build_heat_diric_edge(bilinear_type & bilinear, linear_type & linear)
{
    bilinear+=on(
            _range= markedfaces(mesh_TTL, {"in1", "in2"}),
            _rhs= linear,
            _element= ut,
            _expr= cst(doption("Modele.Tamb"))
            );
}

/// \fn build_heat_time_linear
/// \brief add the temporal term
/// \param linear : reference of the linear form of the low solution
/// \param u_prec : the repartition of the heat for previous step
/// \param dt : time step
    template<typename linear_type>
void Projet::build_heat_time_linear(linear_type & linear, element_heat_type & u_prec, double dt)
{
    linear+= integrate(
            _range= elements(mesh_TTL),
            _expr= idv(rc) /dt * id(u) * idv(u_prec)
            );
}



/// \fn init_beta
/// \brief init the value of the convection
    template<typename myexpr_type>
void Projet::init_beta(myexpr_type beta_expr)
{
    beta.on(
            _range= markedelements(mesh_TTL,"AIR"),
            _expr= beta_expr
           );
}

/// \fn build_heat_bilinear_stab
/// \brief use method of Galerkin Least Square for stabilise the system
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Projet::build_heat_stab(bilinear_type & bilinear, linear_type & linear, element_conv_type conv, myexpr_type Q)
{
    auto L= -idv(k)*laplacian(u) + idv(rc)*( grad(u)*idv(conv) );
    auto Lt= -idv(k)*laplaciant(ut) +  idv(rc) * ( gradt(ut)*idv(conv) );
    auto delta= doption("Modele.epsilon") * cst(1.)/( 1./h() + idv(k)/(h()*h()) );
    bilinear+=integrate(
            _range= elements(mesh_TTL),
            _expr= delta*L*Lt
            );    
    linear+=integrate(
            _range= markedelements(mesh_TTL,{"IC1", "IC2"}),
            _expr= L*delta*Q
            );
}

/// \fn build_heat_static
/// \brief build the matrix matrix_heat for the term with no-time dependance
    template<typename bilinear_type>
void Projet::build_heat_static(bilinear_type & bilinear, double dt)
{
    bilinear+= integrate(
            _range= elements(mesh_TTL),
            _expr= idv(k) * inner(grad(u), gradt(ut))
            );// diffusion of heat
    if( dt>0 )
    {
        bilinear+= integrate(
                _range= elements(mesh_TTL),
                _expr= idv(rc) / dt * id(u) * idt(ut)
                ); // term on time
    }
}


/// \fn build_heat_dynamic
/// \brief add the part dependant ot time
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Projet::build_heat_dynamic(bilinear_type & bilinear, linear_type & linear, element_conv_type conv, myexpr_type Q)
{
    //build_heat_proc(linear,Q);
    linear+= integrate(
            _range= markedelements(mesh_TTL,{"IC1", "IC2"}),
            _expr= Q*id(u)
            );

    if( modele!=no_areation )
    {
        build_heat_conv(bilinear, beta);
        if(boption("Modele.GaLS"))
            build_heat_stab(bilinear, linear, beta, Q);
    }
}

/// \fn build_heat_dynamic
/// \brief add the part dependant ot time
    template<typename bilinear_type, typename linear_type, typename myexpr_type>
void Projet::build_heat_dynamic(bilinear_type & bilinear, linear_type & linear, element_heat_type u_prec,double dt, element_conv_type conv, myexpr_type Q)
{
    linear+= integrate(
            _range= markedelements(mesh_TTL,{"IC1", "IC2"}),
            _expr= Q*id(u)
            );

    if( modele != no_areation )
    {
        build_heat_conv(bilinear, beta);
        if(boption("Modele.GaLS"))
            build_heat_stab(bilinear, linear, beta, Q);
    }

    linear+= integrate(
            _range= elements(mesh_TTL),
            _expr= idv(rc) /dt * id(u) * idv(u_prec)
            );// term of time's memory 

}


/// \fn build_fluid_static
/// \brief add to the bilinear form the static part
/// \f$ += \int_\Omega \bigl(\frac{\nabla u+\nabla u^T}{2}\bigl):\nabla v\f$
/// \f$ += \int_\Omega q\nabla u+p\nabla v\f$
/// if dt is not the default value : 
/// \f$ += \int_\Omega \rho \frac{uv}{dt} \f$
    template<typename bilinear_type>
void Projet::build_fluid_static(
        bilinear_type & bilinear, 
        double dt)
{
    auto v= fluid. template element<0>();
    auto vt= fluidt.element<0>();
    auto p= fluid.element<1>();
    auto pt= fluidt.element<1>();

    //auto mat_id= eye<FEELPP_DIM, FEELPP_DIM>();

    Feel::cout<<"test "<<doption("Air.mu")<<"\n";
     
    Feel::cout<<"static point 1\n";
    //bilinear = integrate(
    //        _range= elements(mesh_AIR),
    //        _expr= doption("Air.mu") * inner( sym(gradt(vt)), grad(v) )
    //        );

    Feel::cout<<"static point 2\n";
    
    bilinear+= integrate(
            _range= elements(mesh_AIR),
            _expr= id(p)*divt(vt)-idt(pt)*div(v)
            );
    
    Feel::cout<<"static point 3\n";
    
    if(dt>0)
        bilinear+= integrate(
                _range= elements(mesh_AIR),
                _expr= doption("Air.rho")*inner(id(v),idt(vt))/dt
                );
    
    Feel::cout<<"static point 4\n";
}

    template<typename linear_type, typename velocity_type, typename pressure_type>
void Projet::build_fluid_dynamic_linear(linear_type & linear, velocity_type v,velocity_type v_prec, double dt)
{

}

    template<typename bilinear_type, typename velocity_type>
void Projet::build_fluid_dynamic_bilinear(bilinear_type & bilinear, velocity_type v, velocity_type vt, velocity_type vl)
{

}




















// RUN //


    template<typename myexpr_type>
void Projet::run_heat(myexpr_type Q,std::string export_name)
{
    tic();
    matrix_heat->zero();
    vector_heat->zero();

    auto bilinear= form2(_test= Th, _trial= Th, _matrix= matrix_heat);
    auto linear= form1(_test= Th, _trial= Th, _vector= vector_heat);
    build_heat_static(bilinear);
    build_heat_dynamic(bilinear, linear, beta, Q);
    build_heat_diric_edge(bilinear, linear);
    toc("build heat");

    tic();
    bilinear.solve(
            _solution= ut,
            _rhs=linear
            );
    toc("solve heat");

    Feel::cout
        <<"les resultats seront affiche dans le dossier '"
        <<export_name<<"'\n";
    auto envoi= exporter(_mesh=mesh, _name=export_name);
    envoi->add("parameter_conductivity",k);
    envoi->add("parameter_capacity",rc);
    envoi->add("heat",ut);
    envoi->save();
}


/// \fn run_heat_step
    template<typename bilinear_type, typename myexpr_type>
void Projet::run_heat_step(bilinear_type & bilinear_static, element_heat_type u_prec, myexpr_type Q, double dt)
{
    tic();
    //matrix_heat->zero();
    //vector_heat->zero();
    auto bilinear= form2( _test= Th, _trial= Th);//, _matrix= matrix_heat);
    auto linear= form1( _test= Th);//, _vector= vector_heat);
    
    bilinear= bilinear_static;
    build_heat_dynamic(bilinear, linear, u_prec, dt, beta, Q);
    build_heat_diric_edge(bilinear, linear);

    toc("build heat");


    tic();
    bilinear.solve(
            _solution= ut,
            _rhs=linear
            );
    toc("solve heat");

}


/// \fn run_heat
    template<typename heat_type,typename conv_type>
void Projet::run_heat(heat_type Q, conv_type beta_expr, double T, std::string export_name)
{
    tic();
    auto envoi= exporter(_mesh= mesh, _name= export_name);
    double dt= doption("Modele.dt");


    tic();
    element_heat_type u_tmp= Th->element();
    //auto f= form2( _test= Th, _trial= Th);
    auto f= form2( _test= Th, _trial= Th, _matrix= matrix_heat);
    build_heat_static( f, dt);
    u_tmp.on( 
            _range= elements(mesh_TTL),
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

        run_heat_step(f, u_tmp, Q, dt);

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



    template<typename myexpr_type>
void Projet::run_fluid(myexpr_type D)
{
    auto file=exporter(_mesh=mesh,_name="stokesSolve");

    auto v= fluid.element<0>();
    auto vt= fluidt.element<0>();
    auto p= fluid.element<1>();
    auto pt= fluidt.element<1>();

    matrix_fluid->zero();
    vector_fluid->zero();
    auto bilinear= form2( _test= Vph, _trial= Vph);
    auto linear= form1( _test= Vph);
    
    Feel::cout<<"run point 1\n";
//#if 0
    build_fluid_static(bilinear);
    Feel::cout<<"run point 2\n";
#if 0
    bilinear+= on(//1
            _range= markedfaces(mesh_AIR, "wall"),
            _rhs= linear,
            _element= vt,
            _expr= zero<FEELPP_DIM,1>()
            );
    Feel::cout<<"run point 3\n";

    bilinear+= on(//2
            _range= markedfaces(mesh_AIR, {"in1", "in2"}),
            _rhs= linear,
            _element= vt,
            _expr= D
            );
    Feel::cout<<"run point 4\n";
    bilinear.solve(_solution= fluid, _rhs= linear);
    file->add("velocity", fluid.element<0>());
    file->add("pressure", fluid.element<1>());
#endif
    auto v_expr= expr<FEELPP_DIM,1>("{0,x*x*x+x*y+y*y*y}:x:y");
    v.on(_range=elements(mesh_AIR), _expr= v_expr);
    file->add("test1",v);
    file->save();

}






// \fn init_edge_in
// \brief init the boundary condition on the edge "in"
// \return a string with the expression of the flow on the edge "in"
std::string Projet::init_edge_in(std::string D)
{
    int n_d= D.find(":t");
    std::ostringstream ostr;
    ostr<<"{0, -"<<D.substr(0,n_d)<<"*";

    if( (modele==no_areation) || (modele==no_stokes) )
    {
        auto mesh_in= createSubmesh(mesh_AIR,markedfaces(mesh,"in2"));
        auto bord_in= boundaryfaces(mesh_in);
        double min=DBL_MAX;
        for(auto const & e:bord_in)
        {
            auto const& elt = unwrap_ref( e );
            double val=elt.point(0)[0];
            ostr<<"(x-"<<val<<")*";
            min= (min>val)? val : min;
        }
        ostr<<"(x>"<<min<<")";
    }
    else
    {
        auto mesh_in= createSubmesh(mesh_AIR, markedfaces(mesh,{"in1","in2"}));
        auto bord_in= boundaryfaces(mesh_in);
        for(auto const & e : bord_in)
        {
            auto const & elt= unwrap_ref(e);
            ostr<<"(x-"<<elt.point(0)[0]<<")*";
        }
        ostr<<1;
    }
    ostr<<"}:x";
    if(n_d!=D.size())
    {   
        Feel::cout<<n_d<<", "<<D.size()<<"\n";
        ostr<<":t";
    }
    return ostr.str();
}

























#endif

#endif
