#! /bin/bash


# FONCTION __________________________________________________________
my_exit()
{
    echo "une erreur c'est produite donc arret du programme"
        echo "les resultats sont : "
        cat $chemin_tab

        cp $chemin_tab $document
        rm $chemin_tab
        echo "le tableau a ete mis a : "$document
        cd $chemin_origine
        exit 1
}

dossier()
{   
    actuel=$PWD
    cd $1
    if [ ! -d "$2" ];
        then
        mkdir $2
        echo "le dossier "$2" a ete cree"
    fi
    cd $actuel
}



# CHOIX _____________________________________________________________
choix="laplacian"
#choix="projet"


chemin_origine=$PWD


nom_fichier="convergence"
nom_tableau=tableau.csv

if [ $choix = "laplacian" ] 
then    
document=~/feel/rapport/laplacian/$nom_fichier
executable=feelpp_qs_laplacian_2d
config=/ssd/atlas_eschbach/feelpp/quickstart/laplacian/circle/circle-all.cfg 
chemin_tab=/data/atlas_home/atlas_eschbach/feel/qs_laplacian/resultat/tableau.csv
type_valeur="Hmax,norme L2,norme H1,temps maillage,temps resolution"

dossier ~/feel rapport
dossier ~/feel/rapport laplacian
dossier ~/feel/rapport/laplacian $nom_fichier
dossier /data/atlas_home/atlas_eschbach/feel/qs_laplacian resultat

elif [ $choix = "projet" ]
then

document=~/feel/rapport/projet/$nom_fichier
executable=feelpp_qs_projet_2d
config=/ssd/atlas_eschbach/feelpp/quickstart/projet/projet2d/projet.cfg
chemin_tab=/data/atlas_home/atlas_eschbach/feel/qs_projet/resultat/tableau.csv
type_valeur="Hmax,norme L2"

dossier ~/feel rapport
dossier ~/feel/rapport laplacian
dossier ~/feel/rapport/laplacian $nom_fichier
dossier /data/atlas_home/atlas_eschbach/feel/qs_projet resultat

else

echo "vous n'avez pas choisi entre laplacian et projet" 
exit 1

fi




#EXECUTION __________________________________________________________


cd /ssd/atlas_eschbach/build/quickstart

touch $chemin_tab
echo $type_valeur>$chemin_tab
#echo "Hmax,norme L2">$chemin_tab




#for h in .5 .3 .2 .1 .06 .04 .02 
for h in .5 .3 .2 .1 .06 .04 .02 .01 .006 .004
do

echo -e "\033[34;1mapplication avec h ="$h"\033[0m"
mpirun -np 4 $executable --config-file $config --gmsh.hsize $h
trap my_exit ERR

done 






echo "les simulations ont ete execute"

cp $chemin_tab  $document/$mon_tableau
rm $chemin_tab
echo "le tableau a ete mis a : "$document
cd $document









echo -e "\033[34;1mmise en graphique( format PNG)\033[0m"

if [ $choix = "laplacian" ] 
then
gnuplot<<EOF 
set term png
set output "graphique.png"
set datafile separator ','
set logscale x
set logscale y
set xlabel "hMax(mm)"
set ylabel "erreur"
set key on inside left top
f(x)=a_L2*x**ordre_L2
fit f(x) "$nom_tableau" using 1:2 via a_L2,ordre_L2
h(x)=a_H1*x**ordre_H1
fit h(x) "$nom_tableau" using 1:3 via a_H1,ordre_H1
set title "erreur selon le pas"
plot "$nom_tableau" using 1:2 with points title " norme L2", x with lines, "$nom_tableau" using 1:3 with points title "norme H1" , x**2 with lines
set grid
quit
EOF

gnuplot<<EOF 
set term png
set output "vitesse.png"
set datafile separator ','
set logscale x
set logscale y
set xlabel "hMax(mm)"
set ylabel "temps(s)"
set key on inside left top
set title "erreur selon le pas"
plot "$nom_tableau" using 1:4 with linespoints title "temps creation mesh", "$nom_tableau" using 1:5 with linespoints title "temps resolution" 
set grid
quit
EOF

else
gnuplot<<EOF 
set term png
set output "graphique.png"
set datafile separator ','
set logscale x
set logscale y
set xlabel "hMax(mm)"
set ylabel "erreur"
set key on inside left top
f(x)=a*x**b
fit f(x) "$nom_tableau" using 1:2 via a,b
set title "erreur selon le pas"
plot "$nom_tableau" using 1:2 with linespoints title " norme L2", x with lines, x**2 with lines
set grid
quit
EOF
fi

cd $chemin_origine
