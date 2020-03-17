License CC BY-NC-SA

# Closure-phases

First generate a bubble field with ipython (for a periodic 2D field with 512x512 pixels, with bubbles of radius 20 pixels, to reach a given filling fraction, allowed to overlap):

from Random_bubbles import *

cube=RandomBubbles(DIM=512, fillingfraction=0.05, radius=20, NDIM = 2, nooverlap = False, periodic=True) 

RandomBubbles.write_ionisation_field(cube) #saves real field in a text file

The text files will be:
   Field_11binary_bubbles_nooverlap=False_radius=20_xhII0.05_N512_2D.txt
  
Then modify the header file SC.h with the correct information: box side length, computation parallelised on 20 cores, range of correlation scales considered (in Mpc)...
Compile
	make 2d for 2D TCF
	make 3d for 3D TCF    
Run the executable:  ./SC.o

The output will be a textfile with the r considered (hard coded, in Mpc), the real part and the imaginary part of the TCF for each r, and the number of modes used to compute the TCF for each r.
