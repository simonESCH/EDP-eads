/**
 * \file md_projet.cpp
 */


#include "md_projet.hpp"







/// \fn main
int main(int argc,char* argv[])
{
    Environment env( 
            _argc=argc,
            _argv=argv,
            _desc=makeOptions(),
            _about=about( _name="projetEDP-M1-CSMI",
                _author="Simon ESCHACH",
                _email=""));


    affichage();


    Projet projet;

    projet.export_init();
    projet.export_param();
    projet.build_heat();
    projet.solve_heat();
    projet.export_heat();
    projet.export_save();
    Feel::cout<<"la temperature de IC1 est "<<projet.get_heat_IC1()<<" K\n";
    
    return 0;
    projet.export_init("cinematics");
    projet.export_evolution();
    projet.export_save();
#if 0
    auto v_expr= expr<FEELPP_DIM,1>( soption(_name = "Modele.flux") );
    auto Q_proc = expr( soption(_name = "Proc.Q") );
    double dt=doption("Modele.dt");


    tic();// construction des maillages
    auto mesh = loadMesh(_mesh = new Mesh<Mesh_type>);
    auto mesh_AIR = createSubmesh(mesh, markedelements(mesh,"AIR"));
    auto mesh_TTL = createSubmesh(mesh, elements(mesh));
    toc("init mesh");












    tic();
    //Feel::cout<<Feel::tag::convex_type::name()<<"\n";

    auto P0 = Pdh<0>(mesh_TTL);

    // coefficient of thermic condictivity
    auto k = P0->element("thermic conductivity");
    k.on(_range= markedelements(mesh_TTL, "PCB"),_expr= cst(doption("PCB.k")));
    k.on(_range= markedelements(mesh_TTL, "IC1"),_expr= cst(doption("Proc.k")));
    k.on(_range= markedelements(mesh_TTL, "IC2"),_expr= cst(doption("Proc.k")));
    k.on(_range= markedelements(mesh_TTL, "AIR"),_expr= cst(doption("Air.k")));

    // coefficient of capacity 
    auto rc = P0->element("capacity");
    rc.on(_range=markedelements(mesh_TTL, "PCB"),_expr= cst(doption("PCB.rc")));
    rc.on(_range=markedelements(mesh_TTL, "IC1"),_expr= cst(doption("Proc.rc")));
    rc.on(_range=markedelements(mesh_TTL, "IC2"),_expr= cst(doption("Proc.rc")));
    rc.on(_range=markedelements(mesh_TTL, "AIR"),_expr= cst(doption("Air.rc"))); 

    toc("init P0  ");


    tic();
    // approximation's Space of Heat
    auto Th = Pch<2>(mesh_TTL);
    auto ut = Th->element("trial Heat");
    auto u  = Th->element("test Heat ");
    auto u_tmp = Th->element("Heat of the time (t-dt)");
    u_tmp.on( 
            _range= elements(mesh_TTL),
            _expr = cst(doption("Modele.Tamb"))
            );//la temperature initiale est la temperature ambiante
    toc("Th       ");


    tic();
    // approximation's Space of the Heat's convection 
    auto Th_vect = Pchv<2>(mesh_TTL);
    auto beta = Th_vect->element("Heat convection");
    toc("Th_vect  ");


    tic();
    // approximation's Space of Navier-Stockes
    auto Vph = THch<1>(mesh_AIR);
    auto F = Vph->element();

    //trial function
    auto vt = F.element<0>("trial flow velocity");
    auto qt = F.element<1>("trial pressure");
    //vt.zero();
    //vt.on(_range=elements(mesh_AIR),_expr=v_expr);


    //test function
    auto v = F.element<0>("test flow velocity");
    auto q = F.element<1>("test pressure");
    toc("Vph      ");











    tic();
    // build of a interpolation's oper√ror within NS and Heat_Convection
    auto flowToConv = opInterpolation( 
            _domainSpace = vt.functionSpace(),
            _imageSpace = beta.functionSpace(),
            _range= markedelements(mesh_TTL, "AIR")
            );
    toc("interpol ");


    tic();
    auto l_Tp = form1(_test = Th);
    //auto l_vp = form1(_test = Vph);
    toc("init Lin ");

    tic();
    auto a_Tp = form2(_trial = Th,_test  = Th);
    toc("init Bil ");





    //partie calcul
#define TEST_CALCUL 1
#if TEST_CALCUL
    for(double t=0; t<=doption("Modele.Tmax"); t+=dt)
    {
        Feel::cout<<"a l'instant t= "<< t << "\n";
        Q_proc.setParameterValues({"t", t});
        v_expr.setParameterValues({"t", t});
        l_Tp.zero();
        a_Tp.zero();

#endif

        vt.on(_range=elements(mesh_AIR), _expr=v_expr);
        tic();
        flowToConv->apply(vt, beta);
        toc("flw->conv");


        tic();
        l_Tp = integrate(
                _range= markedelements(mesh_TTL, "IC1"),
                _expr=  Q_proc *id(u));//chauffe du processeur IC1

        l_Tp+= integrate(
                _range= markedelements(mesh_TTL, "IC2"),
                _expr=  Q_proc *id(u));//chauffe du processeur IC2
#if TEST_CALCUL
        if(doption("Modele.Tmax")<doption("Modele.dt"))
            l_Tp+= integrate(
                    _range= elements(mesh_TTL),
                    _expr= idv(rc) * idv(u_tmp)/dt * id(u)
                    );//ajout du terme de temps
#endif
        toc("Lin temp ");


        tic();
        a_Tp = integrate( 
                _range= elements(mesh_TTL), 
                _expr= idv(k) * gradt(ut) * trans(grad(u))
                );//diffusion de la chaleur
        if (soption("Modele.choix").compare("NO_AREATION"))
        {
            a_Tp+= integrate( 
                    _range= markedelements(mesh_TTL, "AIR"),
                    _expr= idv(rc) *gradt(ut) *idv(beta) *id(u)
                    );//convection de la temperature
            a_Tp+= integrate( 
                    _range= markedelements(mesh_TTL,"out"),
                    _expr= idv(rc) * inner(idv(beta),oneY()) * idt(ut) *id(u)
                    );//condition de sortie


            if(boption("Modele.GaLS"))
            {
                tic();
                auto L = idv(rc) *( grad(u)*idv(beta) ) - idv(k)*laplacian(u);
                auto Lt = idv(rc) *( gradt(ut)*idv(beta) ) -idv(k)*laplaciant(ut);
                double eps = doption("Modele.epsilon");
                auto delta = cst(1.) / (1./h() + eps/(h()*h()));
                a_Tp+= integrate(
                        _range= markedelements(mesh_TTL,"AIR"),
                        _expr= delta * L *Lt );
                toc("stab GaLS");
            }
        }
#if TEST_CALCUL
        a_Tp+=integrate(
                _range= elements(mesh_TTL),
                _expr= idv(rc) * idt(ut)/dt *id(u)
                );
#endif
        a_Tp+= on (
                _range= markedfaces(mesh_TTL,{"in1","in2"}), 
                _rhs = l_Tp,
                _element = ut,
                _expr=cst(doption("Modele.Tamb"))
                );//condition de dirichlet a l'entre d'air
        toc("Bil temp ");









        tic();
        a_Tp.solve(_rhs = l_Tp, _solution = ut);
        toc("heat equation is solved");


        //affichage des temperatures
        tic();
        auto Tproc1 = mean(
                _range= markedelements(mesh_TTL,"IC1"),
                _expr= idv(ut));//temperature du processeur 1
        auto Tproc2 = mean(
                _range= markedelements(mesh_TTL,"IC2"),
                _expr= idv(ut));//temperature du processeur 2
        auto T_out = mean(
                _range= markedfaces(mesh_TTL,"out"),
                _expr= idv(ut));//temperature de la sortie d'air
        toc("param out");

        if(Environment::isMasterRank())
            std::cout<<"\t_______________________________________"
                <<"\n\t|temperature de IC1    ::"<< Tproc1.value()<<" K"
                <<"\n\t|temperature de IC2    ::"<< Tproc2.value()<<" K"
                <<"\n\t|temperature de sortie ::"<< T_out.value()<<" K"
                <<"\n\t_______________________________________"<<std::endl;

        tic();
        auto test_expr = idv(rc) * (gradv(ut) * idv(beta)) - idv(k) *laplacianv(ut);
        auto testbis = P0->element("test proc");
        auto test =P0->element("test");
        testbis.zero();
        testbis.on(_range= markedelements(mesh_TTL,{"IC1","IC2"}), _expr=-Q_proc);
        test.on(_range=elements(mesh_TTL), _expr=test_expr);

        toc("test sol");

        tic();
        auto e =exporter(mesh,_name="visualSolve");
        e->add("temperature",ut);
        e->add("flow",vt);
        e->add("pressure",qt);
        //e->add("test",test-testbis);
        e->save();
        toc("exportSol");

#if TEST_CALCUL
        u_tmp=ut;
    }
#endif

    auto lu_L2 = normL2(
            _range=elements(mesh_TTL),
            _expr= idv(test)-idv(testbis));
    Feel::cout<<"norme L2 de Lu : "<<lu_L2<<"\n";


    tic();
    auto exp = exporter(mesh_TTL,_name="visualParameter");
    exp->add("p_capacity",rc);
    exp->add("p_conductivity",k);
    exp->save();
    toc("exported parameter");
#endif
    return 0;
}







