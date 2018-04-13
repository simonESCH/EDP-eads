/**
 * \file md_projet.cpp
 */


//#include "md_projet.hpp"
//#include "md_ns.hpp"
//#include "md_heat.hpp"
#include "projet_eads.hpp"
//#include "projet_eads2.hpp"




/// \fn main
int main(int argc, char* argv[])
{
    Environment env( 
            _argc= argc, 
            _argv= argv, 
            _desc= makeOptions(), 
            _about= about( _name= "projetEDP-M1-CSMI", 
                _author= "Simon ESCHACH", 
                _email= ""));

    Projet eads;
    eads.run();
    return 0;
}

