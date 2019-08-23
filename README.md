License CC BY-NC-SA

# Closure-phases

First generate a bubble field with ipython (for a periodic 2D field with 512x512 pixels, with 20 bubbles of radius 20 pixels, allowed to overlap):

from Binary_bubbles_nb import *

cube=RandomBubbles(DIM=512, nb=20, radius=20, NDIM = 2, nooverlap = False, periodic=True) 

RandomBubbles.write_ionisation_field(cube) #saves real field in a text file

RandomBubbles.write_ionisation_field_k(cube) #saves FT in two text files, one for the real part, one for the imaginary part

The text files will be:
   Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.100_N512_2D.txt
   Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.100_N512_2D_imagpart.txt
   Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.100_N512_2D_realpart.txt
  
Then compile the .cpp file with the correct information: box side length to 400Mpc, computation parallelised on 20 cores:

g++ -std=c++11 -fopenmp -DSIZE=512 -DLENGTH=400 -DNTHREADS=20 -Dfilename=\"Field_20binary_bubbles_nooverlap=False_radius=20_xhII0.100_N512_2D\" -o ./SC.o ./Spherical_correlations2d.cpp
    
Run the executable:  ./SC.o

The output will be a textfile with the r considered (hard coded, in Mpc), the real part and the imaginary part of the TCF for each r.
