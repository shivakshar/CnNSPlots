echo "Hello"
echo "What is the detector threshold in eVnr?"
read th
echo "What is the overall avg Signal to Background ratio?"
read s_b
echo "What is the systematic error (b/w 0 to 1)?"
read sys
echo "How many energy bins are considered? (10 or greater recommended)"
read bins
echo "How many kgs of Germanium detector are you using?"
read m
echo "How far is the detector from the reactor in meters?"
read L
echo "Great! Plotting..."
python chi_plotter.py $th $s_b $sys $bins $L $m
echo "	Chi-square plot done"
python cs_plotter.py $th $s_b $sys $bins $L $m
echo "	Cross Section Resolution plot done"
python SterileOsscillation.py $th $s_b $sys $bins $L $m
echo "	Sterile Sensitivity plot done"
echo "You will find all the plots in the 'Plots' folder"
echo "You should see an interactive window that shows # of CNS\n with and without Sterile Oscillations soon..."
python SterileDraw.py $th $s_b $sys $bins $L $m  
