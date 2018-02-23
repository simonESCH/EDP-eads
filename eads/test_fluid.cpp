#include "fluide.hpp"




inline
    po::options_description
makeOpt()
{
    po::options_description myapplOpt("mes options");
    myapplOpt.add_options()
        ("Air.mu", po::value<double>()->default_value(1.8e-5),"viscosity")
        ("Air.rho", po::value<double>()->default_value(1.184),"density")
        ("Modele.flux",po::value<std::string>()->default_value("{0,1}"),"flow")
        ("Modele.modele", po::value<std::string>()->default_value("modele0"))
        ("Time.time",po::value<bool>()->default_value(false))
        ("Time.dt",po::value<double>()->default_value(.01),"dt")
        ("Time.Tfinal",po::value<double>()->default_value(10),"temps final");
    myapplOpt.add(backend_options("backend_fluid"));
    return myapplOpt.add(feel_options());
}


// mise en place du maillage
    gmsh_ptrtype
createCarre()
{
    gmsh_ptrtype desc(new Gmsh);

    std::ostringstream ostr;
    ostr
        << "// test avec l'exemple du carre :\n"
        << "// On a une boite carre et on soufle sur l'un des bords.\n"
        << "// On doit alors avoir un tourbillon centrale et des petits\n"
        << "// tourbillons dans les coins.\n\n"


        << desc->preamble() << "\n\n"

        << "h = "<<doption("gmsh.hsize")<<"//m;\n"
        << "Point(1)    = {0, 0, 0, h};\n"
        << "Point(2)    = {2, 0, 0, h};\n"
        << "Point(3)    = {2, 2, 0, h};\n"
        << "Point(4)    = {1, 2, 0, h};\n"
        << "Point(5)    = {0, 2, 0, h};\n"
        << "\n"

        << "Line(1) = {1, 2};\n"
        << "Line(2) = {2, 3};\n"
        << "Line(3) = {3, 4};\n"
        << "Line(4) = {4, 5};\n"
        << "Line(5) = {5, 1};\n"

        << "Line Loop(10) = {1, 2, 3, 4, 5};\n"
        << "Plane Surface(15) = {10};\n"

        << "Physical Line(\"in1\") = {3};\n"
        << "Physical Line(\"in2\") = {4};\n"
        << "Physical Line(\"borderfluid\") = {1,2,5};\n"
        << "Physical Surface(\"AIR\") = {15};\n";

    std::ostringstream nameStr;
    nameStr << "geo_test_fluid";
    //Feel::cout << ostr.str();

    desc->setCharacteristicLength(doption("gmsh.hsize"));
    desc->setDescription(ostr.str());

    return desc;


}



int main(int argc,char* argv[])
{

    Environment env( 
            _argc= argc, 
            _argv= argv, 
            _desc= makeOpt(), 
            _about= about( _name= "test_fluid_M1-CSMI", 
                _author= "Simon ESCHBACH", 
                _email= ""));

    auto mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>,
            _desc= createCarre()
            );
    auto modele=init_modele();
    auto souffle= expr<FEELPP_DIM,1>(soption("Modele.flux"));

    NavierStokes fluid;
    fluid.init(mesh, modele);

    if(! boption("Time.time"))
        fluid.run(souffle);
    else
    {
        Feel::cout << "on se concentre deja sur le cas statique\n";
    }
    //auto exp= exporter(_mesh= mesh, _name= "lala");
    //exp->add("velocity", fluid.m_fluid.element<0>());
    //exp->save();

    return 0;
}
