#include "fluide.hpp"


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
    auto souffle= expr<FEELPP_DIM,1>("{1,0}")

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
