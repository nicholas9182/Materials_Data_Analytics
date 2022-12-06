cd FES_CM1
rm FES*
plumed sum_hills --hills ../HILLS --outfile FES --kt 2.47 --stride 1000 --idw CM1 --bin 200,200

cd ../FES_CM1_D1
rm FES*
plumed sum_hills --hills ../HILLS --outfile FES --kt 2.47 --mintozero --bin 200,200

cd ../FES_D1
rm FES*
plumed sum_hills --hills ../HILLS --outfile FES --kt 2.47 --stride 1000 --idw D1 --bin 200,200
