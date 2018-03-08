#include "fluide.hpp"





po::options_description
makeFluidOption()
{
    po::options_description myappOptions( "My app options" );
    myappOptions.add_options()
        ("Geo.largsize", po::value< double>()->default_value( 2 ), 
         "largeur")
        ("Geo.longsize", po::value< double>()->default_value( 2 ), 
         "largeur")

        ("Fluid.rho", po::value< double>()->default_value( .1 ), 
         "")
        ("Fluid.mu", po::value< double>()->default_value( 1.8e-5 ), 
         "")
        ("Fluid.debit", po::value< double>()->default_value( 1 ), 
         "")
        ("Modele.modele", po::value<std::string>()->default_value("modele2"), 
         "")

        ("Time.time", po::value< bool>()->default_value(false), 
         "")
        ("Time.Tfinal", po::value< double>()->default_value( 0 ), 
         "")
        ("Time.dt", po::value< double>()->default_value( 1 ), 
         "")
        ("Time.save", po::value< double>()->default_value( 0 ), 
         "")
        
        ("Exporter.save", po::value<std::string>()->default_value("test_fluid"), 
         "");


}

// mise en place du maillage
gmsh_ptrtype
createGMSH()
{
    gmsh_ptrtype desc(new Gmsh);

    std::ostringstream ostr;
    ostr
        << "// test avec l'exemple du carre :\n"
        << "// On a une boite carre et on soufle sur l'un des bords.\n"
        << "// On doit alors avoir un tourbillon centrale et des petits\n"
        << "// tourbillons dans les coins.\n\n"

        << "h = "<<doption("gmsh.hsize")<<";//m\n"

        << desc->preamble() << "\n\n"

        << "larg = " << doption("Geo.largsize") << ";\n"
        << "long = " << doption("Geo.longsize") << ";\n"

        << "Point(1)    = {0,      0,    0, h};\n"
        << "Point(2)    = {larg,   0,    0, h};\n"
        << "Point(3)    = {larg,   long, 0, h};\n"
        << "Point(4)    = {larg/2, long, 0, h};\n"
        << "Point(5)    = {0,      long, 0, h};\n"
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
        << "Physical Line(\"borderfluid\") = {1,2,5};\n";

    std::ostringstream nameStr;
    nameStr << "geo_test_fluid";

    desc->setCharacteristicLength(doption("gmsh.hsize"));
    desc->setDescription(ostr_desc.str());

    return desc;


}



int main(int argc,char* argv[])
{

    Environment env( 
            _argc= argc, 
            _argv= argv, 
            _desc= makeOptions(), 
            _about= about( _name= "test_fluid_M1-CSMI", 
                _author= "Simon ESCHACH", 
                _email= ""));

    auto mesh= createGMSHMesh(
            _mesh= new Mesh<MyMesh_type>,
            _desc= createGMSH()
            );
    auto modele=init_modele();

    std::ostringstream ostr_souffle;
    ostr_souffle << "{" << doption("Fluid.debit") << ",0}";
    auto souffle= expr<FEELPP_DIM,1>(ostr_souffle.str());

    NavierStokes fluid;
    fluid.init(mesh,modele);

    if(! boption("Time.time"))
        run(souffle);
    else
    {
        feel::cout << "on se concentre deja sur le cas statique\n";
    }


    return 0;
}
