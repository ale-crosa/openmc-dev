name="res_0.125"

# Diamond
cd Diamond
cd ./2;   openmc > "$name.out"
cd ../5;  openmc > "$name.out"
cd ../10; openmc > "$name.out"
cd ../20; openmc > "$name.out"
cd ../..

# Gyroid
cd Gyroid
cd ./2;   openmc > "$name.out"
cd ../5;  openmc > "$name.out"
cd ../10; openmc > "$name.out"
cd ../20; openmc > "$name.out"
cd ../..

# Schwarz_p
cd Schwarz_p
cd ./2;   openmc > "$name.out"
cd ../5;  openmc > "$name.out"
cd ../10; openmc > "$name.out"
cd ../20; openmc > "$name.out"
cd ../..