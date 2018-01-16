#include <feel/feel.hpp>

//#include <feel/feelmesh/mesh2d.hpp>
//#include <feel/feelmesh/simplex.hpp>
//#include <feel/feelmesh/complement.hpp>
//#include <feel/feelmesh/marker.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
//#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <string>
#include <fstream>

//#define ORDRE_MESH 1
#define ORDRE_TEMPERATURE 1
/// \warning marche pour (MESH, TEMPERATURE) =(1, 1), (2, 2)      /(1,2) et (2,1) ne marche pas
//(ordre mesH=2 ET THch=1 marche pas, marche avec mesh=1)



using namespace Feel;
//simplex peut aussi s'ecrire :
// Simplex< d,o> avec d la dimension et o l'ordre de l'espace d'approximation
using typeMesh = Simplex< FEELPP_DIM>;//, ORDRE_MESH >;//<FEELPP_DIM, ORDRE_MESH, FEELPP_DIM> ;


//double variation(double t, double t_param, double val_max)
//{
//    return val_max *(1- math::exp(-t/t_param));
//}




//template<typename SpacePtrType>

//}




//_____________________________________________________________________________

int main(int argc, char* argv[])
{

    po::options_description my_options( "project options");

    my_options.add(backend_options("backend air"));
    my_options.add(backend_options("backend temperature"));

    // declaration of the option in the file .cfg 
    my_options.add_options()
        ("Air.rc", po::value< double>()->default_value( 0 ), "Rho*Capa of air")
        ("Proc.rc", po::value< double>()->default_value( 0 ), "Rho*Capa of processor")
        ("PCB.rc", po::value< double>()->default_value( 0 ), "Rho*Capa of motherboard")

        ("Air.k", po::value< double>()->default_value( 1), "conductivity thermic of air")
        ("Proc.k", po::value< double>()->default_value( 1), "conductivity thermic of processor")
        ("PCB.k", po::value< double>()->default_value( 1), "conductivity thermic  of motherboard")

        ("Modele.Tamb", po::value< double>()->default_value(293.15), "temperature of reference")
        ("Proc.Q", po::value< std::string>()->default_value("0."), "Quantity of heat of the processor")
        ("Modele.flux",po::value< std::string>()->default_value(""), "profile de poiseuille")
        ("Modele.choix",po::value< std::string>()->default_value("NO_AREATION"), "choix du Modele")

        ("Modele.Tmax",po::value< double>()->default_value(1), "temps maximal")
        ("Modele.dt",po::value< double>()->default_value(1), "pas de temps")
        ("Modele.GaLS",po::value<bool>()->default_value( false ), "activation dee la stabilisation GaLS")

        ("Exporter.save",po::value< std::string>()->default_value("")," chemin de la sauvegarde")
        ("Exporter.load",po::value< std::string>()->default_value(""),"chemin de l'approximation fine")
        //("Exporter.name",po::value< std::string>()->default_value(""), "nom du fichier ou sera le .case")
        ;

    //declaration of the environment
    Environment env( _argc=argc, _argv=argv,
            _desc=my_options,
            _about=about(  _name="qs_projet",
                _author="",
                _email="feelpp-devel@feelpp.org"));


    double T_ambiante = doption("Modele.Tamb");
    auto v_expr = expr<FEELPP_DIM,1>(soption(_name="Modele.flux"));
    auto Q_proc = expr(soption( _name = "Proc.Q"));

    //print of the using parameter 
    Feel::cout<< "\n\t========================"
        << "\n\t|Air.rc : " << doption("Air.rc") << " J/K/mm^3"
        << "\n\t|Proc.rc: " << doption("Proc.rc") << " J/K/mm^3"
        << "\n\t|PCB.rc : " << doption("PCB.rc") << " J/K/mm^3"
        << "\n\t|======================="
        << "\n\t|Air.k  : " << doption("Air.k") << " W/K/mm"
        << "\n\t|Proc.k : " << doption("Proc.k") << " W/K/mm"
        << "\n\t|PCB.k  : " << doption("PCB.k") << " W/K/mm"
        << "\n\t|======================="
        << "\n\t|Tamb   : " << T_ambiante << " K"
        << "\n\t|chauffe: " << soption("Proc.Q") << " en W/mm^3"
        << "\n\t|ventil : " << soption("Modele.flux")
        << "\n\t|choix  : " << soption("Modele.choix");
    
    if(soption("Exporter.load").compare(""))
        Feel::cout<<"\n\t|load   :"<<soption("Exporter.load");
    if(soption("Exporter.save").compare(""))
        Feel::cout<<"\n\t|save   :"<<soption("Exporter.save");

    Feel::cout<< "\n\t|======================="
        << "\n\t|hsize  : " << doption("gmsh.hsize")
        << "\n\t|Tmax   : " << doption("Modele.Tmax")
        << "\n\t|dt     : " << doption("Modele.dt")
        << "\n\t========================\n\n";


    // declaration of the backend (test)
    auto backend_air = backend( _name = "backend air");
    auto backend_ttl = backend( _name = "backend temperature");

    
    //Feel::cout<<"mesh_type : "<<Feel::tag::typeMesh::convex_type::name()<<"\n";


    tic();
    //creation of the mesh
    auto mesh = loadMesh(_mesh=new Mesh< typeMesh >);
    auto mesh_Air = createSubmesh(mesh, markedelements(mesh,"AIR"));
    //auto mesh_ttl = createSubmesh(mesh, elements(mesh));
    toc("loadMesh");


    tic();
    //initialisation of the space for equation of navier-stockes
    //auto Vph = new FunctionSpace< Mesh<typeMesh>, 
    //     bases< Lagrange<2, Vectorial>, 
    //     Lagrange<1, Scalar>> >
    //         (mesh_Air);
    //auto F = Vph->element();
    auto Vph = THch<1>(mesh_Air);
    auto F = Vph->element();

    //Fonction trial
    auto vt = F.element<0>("vitesse vent");
    auto qt = F.element<1>("pression");

    //Condition of the systeme in the time t-dt
    auto Fk = Vph->element();
    auto v_tmp = Fk.element<0>("vitesse temporaire");

    //Fonction test
    auto v = F.element<0>("vitesse test");
    auto q = F.element<1>("pression test");

    qt.zero();
    vt.on(_range = markedelements(mesh_Air,"POI"), _expr = v_expr); 
    toc("Vph");


    tic();
    //initialisation of Temperature's space and function test
    auto Th = Pch< ORDRE_TEMPERATURE >(mesh);//P continue
    auto ut = Th->element("temperature");
    auto u = Th->element("temperature test");
    ut.on(_range = elements(mesh), _expr = cst(T_ambiante));
    auto app_sol = Th->element("load");
    toc("Th");


    tic();
    //initialision of the vt's interpolate
    auto Th_vect = Pchv< ORDRE_TEMPERATURE >(mesh);
    auto beta = Th_vect->element();
    toc("Th_vect");


    tic();
    // interpolation  permet de faire passer le courant de mesh_Air dans le maillage mesh 
    /// \bug je n'arrive pas a le faire fonctionner quand test avec Pchv, Pdhv
    auto AirToBeta = opInterpolation( _domainSpace = vt.functionSpace(), _imageSpace = beta.functionSpace(), _range = markedelements(mesh,"AIR") );
    //auto AirToBeta = operatorinterpolation(vt.functionSpace(), beta.functionSpace, 
    //beta.on(_range = markedelements(mesh, "POI"), _expr = idv(vt) );
    //AirToBeta->apply(vt, beta);
    //beta.on( _range = markedelements(mesh, "POI"), _expr = v_expr);
    toc("interpolation");


    if(soption("Exporter.load").compare(""))
    {
        tic();
        auto mesh_tmp = loadMesh(_mesh = new Mesh<typeMesh>,_filename= soption("Exporter.load")+"mesh.json");
        auto TMP = Pch<1>(mesh_tmp);
        auto tmp = TMP->element();
        tmp.load(_path=soption("Exporter.load"), _type="hdf5");
        app_sol.on( _range=elements(mesh), _expr=idv(tmp));
        toc("load affine");
    }//end if LOAD

    tic();
    auto l = form1( _test=Th, _backend = backend_ttl);
    auto rapport_qk = Q_proc / doption("Proc.k");
    l = integrate(_range = markedelements(mesh, "IC1"), _expr = rapport_qk * id(u));
    l+= integrate(_range = markedelements(mesh, "IC2"), _expr = rapport_qk * id(u));
    toc("l");


    tic();
    
    auto a = form2( _trial=Th, _test=Th, _backend = backend_ttl);
    a = integrate( _range = elements(mesh), _expr =  gradt(ut)* trans(grad(u)));//diffusion de chaleur

    if(soption("Modele.choix").compare("NO_AREATION"))
    {
        auto N_out = oneY(); 
        auto rapport_rck = cst(doption("Air.rc")/doption("Air.k"));

        a+= integrate( _range = markedelements(mesh, "POI"), _expr = rapport_rck * (gradt(ut) * idv(beta)) * id(u));//convection de la temperature
        a+= integrate( _range = markedfaces(mesh, "out"), _expr = rapport_rck * inner(idv(beta), N_out) * idt(ut) * id(u));//condition de sortie
#if 0
        if (boption("Modele.stabilisation"))
        {
            tic();
            /// \warning grad(u) est un vecteur ligne donc pas produit scalaire
            auto L = rapport_rck * grad(u) * idv(beta) - laplacian(u);
            auto Lt = rapport_rck * gradt(ut) * idv(beta) - laplaciant(ut);
            double epsilon = 1;
            auto h_mesh = h();
            auto delta = cst(1.)/(1./h_mesh + epsilon/(h_mesh*h_mesh));
            a+= integrate( _range = markedelements(mesh,"AIR"), 
                    _expr = delta * L * Lt );
            l+= integrate( _range = markedelements(mesh, "AIR"), 
                    _expr = delta * L * cst(rapport_qk));
            toc("stabilisation");
        }
#endif
    }//end if AREATION
    /// \remarque div(bT)=div(b)T+b*grad(T)
    /// \warning N() est definit par rapport au .geo donc eviter d'utiliser N()

    a+= on( _range = markedfaces(mesh, "in"), _rhs = l, _element = ut, _expr = cst(T_ambiante));//condition d'entre
    toc("a");


    tic();
    a.solveb( _rhs = l, _solution = ut, _backend=backend_ttl);
    toc("a is solved");


    tic();
    auto e =exporter(_mesh=mesh);
    e->add("temperature", ut);
    if(soption("Exporter.load").compare(""))
        e->add("temperatureAff",app_sol);
    e->add("flot", vt);
    e->save();
    toc("export");


    tic();
    auto Tproc1 = mean( _range = markedelements(mesh,"IC1"), _expr = idv(ut));
    Feel::cout<<"Tproc1\n";
    auto Tproc2 = mean( _range = markedelements(mesh,"IC2"), _expr = idv(ut));
    Feel::cout<<"Tproc2\n";
    //auto Tout = mean( _range = markedfaces(mesh,"out"), _expr = idv(ut));
    toc("mean proc");


    //affiche la temperature des processeurs
    if(Environment::isMasterRank())
        std::cout<< "________________________.._________" << std::endl
            << "| temperature de IC1    ||" << Tproc1.value() << " K" << std::endl
            << "| temperature de IC2    ||" << Tproc2.value() << " K" << std::endl
            //<< "| temperature de sortie ||" << Tout.value() << " K" << std::endl
            << "|_______________________||_________" << std::endl;


    // recuperer l'erreur par rapport a l'approximation fine de load dans un fivhier .csv
    if(soption("Exporter.load").compare(""))
    {

        auto hMax = mesh->hMax();
        auto diffL2 = normL2(_range=elements(mesh), _expr=idv(ut)-idv(app_sol));
        Feel::cout<< "la difference avec l'approximation fine est :  "
            << diffL2 << " K\n";
    //if(Environment::isMasterRank())
    //{
    //    std::ofstream fichier("/data/atlas_home/atlas_eschbach/feel/qs_projet/resultat/tableau.csv",ios::out|ios::app);
    //    if(fichier)
    //    {
    //        fichier<<hMax<< "," << diffL2<< std::endl;
    //        fichier.close();
    //    }
    //    else
    //        std::cerr<< "ERREUR : le fichier ~/feel/qs_projet/resultat/tableau.csv ne c'est pas ouvert" << std::endl;
    //}
    }//end if TABLEAU

    if(soption("Exporter.save").compare(""))
    {
        tic();//ATTENTION ON NE SAUVE QUE SI ON VERIFIE LE BON NOMBRE DE COEUR
        ut.save(_path=soption("Exporter.save"), _type="hdf5");
        mesh->saveHDF5(soption("Exporter.save")+"mesh.json");
        toc("sauvegarde type = hdf5");
    }//end if SAVE

    return 0;
}
