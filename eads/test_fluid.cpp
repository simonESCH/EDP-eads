#include "fluide.hpp"




inline
    po::options_description
makeOpt()
{
    po::options_description myapplOpt("mes options");
    myapplOpt.add_options()
        ("Geo.longsize", po::value< double>()->default_value( 1 ), 
         "largeur")
        ("Geo.thicksize", po::value< double>()->default_value( 1 ), 
         "longueur")
        //("gmsh.rebuild", po::value<bool>()->default_value(false),"")

        ("Fluid.rho", po::value< double>()->default_value( 1 ))
        ("Fluid.mu", po::value< double>()->default_value( 1.8e-5 ))
        ("Fluid.debit", po::value< double>()->default_value( 1 ))

        ("Time.time", po::value< bool>()->default_value(false))
        ("Time.Tfinal", po::value< double>()->default_value( 0 ))
        ("Time.dt", po::value< double>()->default_value( 1 ))
        ("Time.save", po::value< double>()->default_value( 0 ))
#if TEST_CONVERGENCE
        ("Test.solution_velocity", po::value< std::string>()->default_value("{0,0}:x:y"))
        ("Test.solution_pressure", po::value< std::string>()->default_value("0:x:y"))
        ("Test.f", po::value< std::string>()->default_value("0:x:y"))
        ("Test.neumann", po::value< std::string>()->default_value("{0,0}:x:y"))

#endif
        ("Exporter.save", po::value<std::string>()->default_value("testFluid"), 
         "");

    myapplOpt.add(backend_options("backend_fluid"));
    return myapplOpt.add(feel_options());
}


gmsh_ptrtype
createKovasznay()
{
    double dt=doption("gmsh.hsize");

    gmsh_ptrtype geomet(new Gmsh);

    CHECK(dt<1.5) << "le pas de temps est trop grand\n";
    geomet->setCharacteristicLength(dt);

    std::ostringstream desc;
    desc
        << "//exemple donne dans le manuel de Feelpp(https://book.feelpp.org/math/fem/#sec:mixed-finite-element)\n"
        << "//  u(x,y)={1-exp(lambda*x)*cos(2*pi*y) , lambda/(2*pi)*exp(lambda*x)*sin(2*pi*y)}\n"
        << "//  p(x,y)=-exp(2*lambda*x)/2\n"
        << "//  f(x,y)={exp(lambda*x)*((lambda^2-4*pi^2)*nu*cos(2*pi*y)-lambda*exp(lambda*x)),\n"
        << "//          exp(lambda*x)*nu*sin(2*pi*y)*(-lambda^2+4*pi^2)                      }\n"
        << "//  lambda=1/(2*mu)-sqrt(1/(4*mu**2)+4*pi**2)\n"
        << "\n//-----------------Preambule------------------\n"
        << geomet->preamble() << "\n"
        << "\n//----------------Description-----------------\n"
        << "Point(1)    = {-0.5,  -0.5,  0, h};\n"
        << "Point(2)    = {1,     -0.5,  0, h};\n"
        << "Point(3)    = {1,     1.5,   0, h};\n"
        << "Point(4)    = {-0.5,  1.5,   0, h};\n"
        << "\n"

        << "Line(1) = {1, 2};\n"
        << "Line(2) = {2, 3};\n"
        << "Line(3) = {3, 4};\n"
        << "Line(4) = {4, 1};\n"

        << "Line Loop(10) = {1, 2, 3, 4};\n"
        << "Plane Surface(15) = {10};\n"

        << "Physical Surface(\"AIR\") = {15};\n\n"
        << "Physical Line(\"in\") = {1,2,3,4};\n";

        std::ostringstream name;
        name << "Kovasznay.geo";
        geomet->setDescription(desc.str());
        geomet->setPrefix(name.str());
        return geomet;
}

// mise en place du maillage
    gmsh_ptrtype
createCarre()
{
    gmsh_ptrtype geomet(new Gmsh);
    geomet->setCharacteristicLength(doption("gmsh.hsize"));

    std::ostringstream ostr;
    ostr
        << "// test avec l'exemple du carre :\n"
        << "// On a une boite carre et on soufle sur l'un des bords.\n"
        << "// On doit alors avoir un tourbillon centrale et des petits\n"
        << "// tourbillons dans les coins.\n\n"


        << geomet->preamble() << "\n\n"

        //<< "h = "<<doption("gmsh.hsize")<<"//m;\n"
        << "largeur = " << doption("Geo.longsize") << ";\n"
        << "longueur = " << doption("Geo.thicksize") << ";\n"

        << "Point(1)    = {0,         0,        0, h};\n"
        << "Point(2)    = {largeur,   0,        0, h};\n"
        << "Point(3)    = {largeur,   longueur, 0, h};\n"
        << "Point(4)    = {0,         longueur, 0, h};\n"
        << "\n"

        << "Line(1) = {1, 2};\n"
        << "Line(2) = {2, 3};\n"
        << "Line(3) = {3, 4};\n"
        << "Line(4) = {4, 1};\n"

        << "Line Loop(10) = {1, 2, 3, 4};\n"
        << "Plane Surface(15) = {10};\n"
        << "Physical Surface(\"AIR\") = {15};\n\n"

        << "Physical Line(\"in\") = {3};\n"
        << "Physical Line(\"borderfluid\") = {1,2,4};\n";
    std::ostringstream nameStr;
    nameStr << "geo_test_fluid_carre.geo";

    //Feel::cout << "description :\n" << ostr.str();


    geomet->setPrefix(nameStr.str());
    geomet->setDescription(ostr.str());
    Feel::cout << "prefix : " << geomet->prefix() << "\n";

    return geomet;


}



int main(int argc,char* argv[])
{

    Environment env( 
            _argc= argc, 
            _argv= argv, 
            _desc= makeOpt(), 
            _about= about( _name= "test_fluid", 
                _author= "Simon ESCHBACH", 
                _email= ""));

    tic();
#if TEST_CONVERGENCE
    auto gmsh= createKovasznay();
    std::string modele_str("modele2");
#else
    auto gmsh= createCarre();
    std::string modele_str("modele3");
#endif

    auto mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>,
            _desc= gmsh,
            _rebuild=boption("gmsh.rebuild")
            );
    toc("init mesh");

    std::ostringstream ostr_exp;
    ostr_exp << soption("Exporter.save");
    if(TEST_CONVERGENCE)
        ostr_exp << "_poiseuille";
    else
        ostr_exp << "_carre";
    if(boption("Time.time"))
        ostr_exp << "_dynamic";
    else
        ostr_exp << "_fixed";

    auto exp= exporter(_mesh=mesh, _prefix="test_fluid", _name=ostr_exp.str());

    auto modele=init_modele(modele_str);

    std::ostringstream ostr_souffle;
#if TEST_CONVERGENCE
    Feel::cout
        << "force    : " << soption("Test.f")
        << "\nneumann  : " << soption("Test.neumann")
        << "\nvitesse  : " << soption("Test.solution_velocity")
        << "\npression : " << soption("Test.solution_pressure");
    ostr_souffle << "{" << doption("Fluid.debit") <<"*y*(" <<doption("Geo.thicksize") << "-y) ,0}:x:y";
#else
    ostr_souffle << "{" << doption("Fluid.debit") << ",0}";
#endif
    auto souffle= expr<FEELPP_DIM,1>(ostr_souffle.str());



    NavierStokes fluid;
    fluid.init(mesh, modele);
    Feel::cout << "Nb Reynold : "
        << doption("Fluid.debit") * doption("Geo.thicksize") * doption("Fluid.rho")/doption("Fluid.mu") << "\n";


    if(! boption("Time.time"))
    {
        fluid.run(souffle);
        exp->add("velocity", fluid.m_fluid.element<0>());
        exp->add("pressure", fluid.m_fluid.element<1>());
        exp->save();

#if TEST_CONVERGENCE
        auto diff_vitesse= idv(fluid.m_fluid.element<0>())-expr<FEELPP_DIM,1>(soption("Test.solution_velocity"));
        auto diff_pressure= idv(fluid.m_fluid.element<1>())-expr(soption("Test.solution_pressure"));
        double errorL2_vitesse=
            integrate(
                    _range= elements(mesh),
                    _expr=inner(diff_vitesse,diff_vitesse)
                    ).evaluate()(0,0);

        double errorL2_pressure=normL2(
                _range= elements(mesh),
                _expr= idv(fluid.m_fluid.element<1>())
                );
        Feel::cout << "vitesse  : " << sqrt(errorL2_vitesse)
            << "\npressure : " << errorL2_pressure << "\n";
#endif

    }
    else
    {
        double t;
        double dt=doption("Time.dt");
        double time_save= doption("Time.save");
        Feel::cout << "number of iteration per save : " << time_save << "\n";
        int cpt_save=1;

        for(t=0; t<doption("Time.Tfinal");t+=dt)
        {
            tic();
            souffle.setParameterValues({{"t",t}});
            //Feel::cout << "marqueur 1\n";
            fluid.run(souffle);
            int momentsave=(int)(t/dt);
            if( (t <= time_save*cpt_save)&&(t+dt > time_save*cpt_save) )
            {
                std::ostringstream ostr;
                ostr << "save: " << cpt_save/time_save; 
                tic();
                exp->step(t)->add("velocity", fluid.m_fluid.element<0>());
                exp->step(t)->add("pressure", fluid.m_fluid.element<1>());
                exp->save();
                toc(ostr.str());
                cpt_save++;
                //fluid.init_matrix();
            }
            std::ostringstream time;
            time << "time: " << std::setw(6) << t << " s";
            toc(time.str());
            Feel::cout << "\n";
        }
    }

    return 0;
}
