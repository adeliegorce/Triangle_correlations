###################################################################
# Class for producing box containing randomly sized disks/bubbles #
###################################################################
#
# Currently doesn't handle periodic boundary conditions
#
# Written for 2D by Jonathan Pritchard (2017) - upgraded to 3D by Adelie Gorce (November 2017)
#

import numpy as np
import numpy.random as npr
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import sqrt, pi
import os
import sys, getopt
import pylab as pl
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import ascii
from astropy.table import Table, Column

class RandomBubbles:
    """
    Produce a random 2D box filled with disks in 2D
    """

    def __init__(self, DIM=256, fillingfraction = 0.3, params = [10.0, 2.0], NDIM = 2, nooverlap = True):
        """
        Initialise a DIM x DIM box with random bubbles until
        fillingfraction of pixels have been labelled

        if NDIM =3 make a 3D box
        """

        self.NDIM = NDIM  #number of dimensions e.g. 2D or 3D
        self.DIM = DIM   #Number of cells per dimension
        self.nooverlap = nooverlap
        self.fillingfraction=fillingfraction
        
        print("initialising a %iD box with %i cells" %(self.NDIM, self.DIM))
        print("filling faction is %f and nooverlap is %i" %(fillingfraction, nooverlap))

        if self.NDIM == 2:
            self.box = np.ones([DIM, DIM])
        elif self.NDIM ==3:
            self.box = np.ones([DIM, DIM, DIM])
        else:
            raise Exception("NDIM must be 2 or 3" % (NDIM))

        self.Rmean = params[0] #mean radius
        self.sigR = params[1] #deviation from mean

        self.bubbles = []
        self.sizes = []

        #Add bubbles to get to target filling fraction
        self.increase_filling_fraction(fillingfraction)


    def summary(self):
        """
        Update summary statistics
        """
        print("Number of bubbles in=", len(self.bubbles))
        print("box mean = ", self.box.mean())
        #self.sizes = np.array(self.sizes)
        overlap = np.sum(np.pi * np.power(np.array(self.sizes), 2.0))
        overlap /= float(np.power(self.DIM, 2.0))
        print("volume in, mean volume, overlap = %f, %f, %f" % (overlap, 1.0 - self.box.mean(), overlap/(1.0 - self.box.mean())))

        #Show the slice
        plt.ion()
        colours=['lightyellow','midnightblue']
        bounds=[0,1]
        cmap = mpl.colors.ListedColormap(colours)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        fig = plt.figure(1, (15., 15.))
        if (self.NDIM==2):
            plt.imshow(self.box,cmap=cmap)
        elif (self.NDIM==3):
            grid = ImageGrid(fig, 111, # similar to subplot(111)
                        nrows_ncols = (2, 2), # creates 2x2 grid of axes
                        axes_pad=0.12, # pad between axes in inch.
                        label_mode = "L",
                        share_all = "False",
                        )
            img0=grid[0].imshow(self.box[:,int(self.DIM/2),:],cmap=cmap)
            img1=grid[1].imshow(self.box[:,int(self.DIM/2)+10,:],cmap=cmap)
            img2=grid[2].imshow(self.box[int(self.DIM/2),:,:],cmap=cmap)
            img3=grid[3].imshow(self.box[:,:,int(self.DIM/2)],cmap=cmap)

        plt.savefig('Bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')
        
        plt.figure()
        radii=np.array(self.sizes)
        values, bins , _ = plt.hist(radii,bins=max(radii)-min(radii))
        plt.title('Radius distribution')
        plt.text(2,values.max(),'Mean radius : '+str(np.mean(radii)))
        plt.savefig('Radius_distribution_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')

    def increase_filling_fraction(self, target):
        """
        Add randomly located and sized bubbles until reach
        desired filling fraction
        """

        if (1.0 - target > self.box.mean()):
            print("Box already more ionised than target")
            return
        
        count=0

        while (self.box.mean() > 1.0 - target):

            #Gaussian distribution of bubble sizes
            #forced to be >+1 to avoid zero radius bubbles
            R = max(1, int(npr.normal(self.Rmean, self.sigR)))
            x = npr.randint(self.DIM, size = self.NDIM)

            # #uniform distribution for bubble locations
            # if self.nooverlap:
            #     y,accept = self.uniform_distribution_nooverlap(R)
            #     if accept==False:
            #         x=[]
            #         continue
            #     else:
            #     	x=y
            # else:
            #     x = npr.randint(self.DIM, size = self.NDIM)
            #     #outputs an array with three values, corresponding to the x- y- and z- coordinates@ of the bubble

            # #Use mask to ID pixels in this bubble (dim condition in bubble_mask)
            # mask = self.bubble_mask(x, R)
            # self.box[mask] = 0

            #Use mask to ID pixels in this ellipse
            mask = self.bubble_mask(x, R)
            test=0
            if self.nooverlap:
                if np.all(self.box[mask]):
                    test=0
                else:
                    test+=1
                #Check overlap with edges of the cube
            if self.noedges:
                for i in range(0,self.NDIM):
                    if ((x[i]<=R) or (self.DIM-x[i]<=R)):
                        test+=1
                        
            if self.noedges or self.nooverlap:            
                if test==0:
                    self.box[mask]=0
                else:
                    pass
            else:
                self.box[mask] = 0

            #Store bubbles so can assemble true size PDF
            #Important if number of bubbles is small
            self.bubbles.append(x)
            self.sizes.append(R)

            count=count+1
            if count % 100 == 0:
                #print(x, R)
                print(self.box.mean())
                print(str(count)+' bubbles')
                print(' ')

        self.summary()


    def uniform_distribution_nooverlap(self, R):
        """
        Uniform distribution for bubbles without overlap
        This tests for overlap with other bubbles and also with the edge of the box
        """
        test=[]
        accept=False
        #generate coordinates for the centre of the bubble, the radius is an input
        x = npr.randint(self.DIM, size = self.NDIM)

		#if first bubble, automatically accepted
        if len(self.bubbles)==0:
            test.append(True)
        else:
            for i in range(0,len(self.bubbles)) :
                if (self.NDIM==3):
                    #test if distance between centres is larger than sum of radii
                    if (sqrt((self.bubbles[i][0]-x[0])**2 + (self.bubbles[i][1]-x[1])**2 +(self.bubbles[i][2]-x[2])**2)>= self.sizes[i]+R):
                        test.append(True)
                    else:
                        test.append(False)
                        pass
                elif (self.NDIM==2):
                    #test if distance between centres is larger than sum of radii
                    if (sqrt((self.bubbles[i][0]-x[0])**2 + (self.bubbles[i][1]-x[1])**2)>= self.sizes[i]+R):
                        test.append(True)
                    else:
                        test.append(False)
                        pass
                else:
                    raise Exception ("NDIM is not 2 or 3")
            
        #Check overlap with edges of the cube
        for i in range(0,self.NDIM):
            if ((x[i]>R) and (self.DIM-x[i]>R)):
                test.append(True)
            else:
                test.append(False)
                pass
			
        if np.all(test):
            accept=True
            return x,accept
        else:
            accept=False
            return x,accept


    def bubble_mask(self, x, R):
        #wrapper to handle different dimensionality
        #NOT FULLY WRITTEN FOR 3D yet
        if (self.NDIM == 2):
            return self.disk_mask(x,R)
        elif(self.NDIM == 3):
            return self.sphere_mask(x,R)
        else:
            raise Exception ("NDIM is not 2 or 3")

    def disk_mask(self, pos, R):
        #This gives pixels within R of the centre (x0, y0)
        #Here just a 2D circle
        #
        #Memory inefficient at present
        #Disks are truncated at boundaries

        #pos is coordinates of the centre of the bubble
        #R is its radius

        full_struct = np.zeros([self.DIM, self.DIM]).astype(np.bool)

        #Creates a disk at centre of array to allow for large R
        #without hitting boundaries
        structsize = int(2*R+5) #avoids to generate another DIM**3 box : just enough to contain the disk
        struct = np.zeros((structsize, structsize))
        x, y = np.indices((structsize, structsize))
        mask = (x - structsize/2)**2 + (y - structsize/2)**2 <= R**2 #puts the disk in the middle of new box
        struct[mask] = 1
        
        #Now work out coordinate shift to move centre to pos	
        xmov = [pos[0] - int(structsize/2),pos[0]+int(structsize/2)+1]
        ymov = [pos[1] - int(structsize/2),pos[1]+int(structsize/2)+1]

        #if struct goes out of the box
        xmin=max(xmov[0],0)
        xmax=min(xmov[1],self.DIM-1)
        #print(xmin,xmax)

        ymin=max(ymov[0],0)
        ymax=min(ymov[1],self.DIM-1)
        #print(ymin,ymax)

        #truncated struct if some part is outside the full struct
        small_struct=struct[abs(xmov[0]-xmin):structsize-abs(xmov[1]-xmax), abs(ymov[0]-ymin):structsize-abs(ymov[1]-ymax)]
        #add to previous box in case some intermediate structures overlap
        full_struct[xmin:xmax,ymin:ymax] = np.add(full_struct[xmin:xmax,ymin:ymax],small_struct.astype(np.bool)) 

        plt.ion()
        plt.imshow(full_struct)
        
        return full_struct.astype(np.bool)


    def disk_mask_periodic(self, pos, n):
        #This gives pixels within n of the centre (x0, y0)
        #Here just a 2D circle, but in such a way that
        #periodic boundary conditions are respected

        #Create a disk at centre of array to allow for large n
        #without hitting boundaries
        struct = np.zeros((self.DIM, self.DIM))
        x, y = np.indices((self.DIM, self.DIM))
        mask = (x - self.DIM/2)**2 + (y - self.DIM/2)**2 <= n**2
        struct[mask] = 1

        #Now work out coordinate shift to move centre to pos
        xmov = self.DIM/2 - pos[0]
        ymov = self.DIM/2 - pos[1]

        #use np.roll to move center to pos respecting perodicity
        struct=np.roll(np.roll(struct,shift=-xmov,axis=0),shift=-ymov,axis=1)
        return struct.astype(np.bool)


    def sphere_mask(self, pos, R):
        #This gives pixels within n of the centre (x0, y0, z0)
        #Here just a 3D sphere
        #This is very slow as wastes lots of time on not needed positions

        full_struct = np.zeros([self.DIM, self.DIM, self.DIM]).astype(np.bool)

        #Create a sphere at centre of array to allow for large R
        #without hitting boundaries
        structsize = int(2*R+5)
        struct = np.zeros((structsize, structsize, structsize))
        x, y, z = np.indices((structsize, structsize, structsize))
        mask = (x - structsize/2)**2 + (y - structsize/2)**2 + (z - structsize/2)**2<= R**2
        struct[mask] = 1

        #Now work out coordinate shift to move centre to pos
        xmov = [pos[0] - int(structsize/2),pos[0]+int(structsize/2)+1]
        ymov = [pos[1] - int(structsize/2),pos[1]+int(structsize/2)+1]
        zmov = [pos[2] - int(structsize/2),pos[2]+int(structsize/2)+1]

        #if struct goes out of the box
        xmin=max(xmov[0],0)
        xmax=min(xmov[1],self.DIM-1)
        #print(xmin,xmax)

        ymin=max(ymov[0],0)
        ymax=min(ymov[1],self.DIM-1)
        #print(ymin,ymax)

        zmin=max(zmov[0],0)
        zmax=min(zmov[1],self.DIM-1)
        #print(zmin,zmax)

        small_struct=struct[abs(xmov[0]-xmin):structsize-abs(xmov[1]-xmax), abs(ymov[0]-ymin):structsize-abs(ymov[1]-ymax),abs(zmov[0]-zmin):structsize-abs(zmov[1]-zmax)] #truncated struct if some part is outside the full struct
        #print(abs(xmov[0]-xmin),structsize-abs(xmov[1]-xmax), abs(ymov[0]-ymin),structsize-abs(ymov[1]-ymax),abs(zmov[0]-zmin),structsize-abs(zmov[1]-zmax))
        full_struct[xmin:xmax,ymin:ymax,zmin:zmax] = np.add(full_struct[xmin:xmax,ymin:ymax,zmin:zmax],small_struct.astype(np.bool)) #add to previous box in case some intermediate structures overlap

        return full_struct

    def phase_analysis(self,box):

        import matplotlib.pyplot as plt
        import numpy as np
        plt.close('all')
        plt.ion()

        meand=box.mean()
        overdensity=(box-meand)/meand

        #Fourier transform
        if (self.NDIM==3) :
            overdensity_k=np.fft.fftn(overdensity,axes=(0,1,2))
        elif self.NDIM==2:
            overdensity_k=np.fft.fftn(overdensity,axes=(0,1))
        #hermitian condition is already verified: density[-1,-1,-1]=density[1,1,1]*

        #power spectrum
        PS=np.power(np.abs(overdensity_k),2)
        #phase spectrum
        phase_k=np.angle(overdensity_k)
        #phase_k=np.fft.fftshift(phase_k)

        if self.NDIM==2:
            Dx=np.diff(phase_k,axis=0)[:,:-1] #removes last columns since they won't be used in the overall increment
            Dy=np.diff(phase_k,axis=1)[:-1,:]
            #overall phase increment
            Dk=np.power(np.power(Dx,2)+np.power(Dy,2),1/2)/sqrt(3)
        elif self.NDIM==3:
            Dx=np.diff(phase_k,axis=0)[:,:-1,:-1] #removes last columns since they won't be used in the overall increment
            Dy=np.diff(phase_k,axis=1)[:-1,:,:-1]
            Dz=np.diff(phase_k,axis=2)[:-1,:-1,:]
            #overall phase increment
            Dk=np.power(np.power(Dx,2)+np.power(Dy,2)+np.power(Dz,2),1/2)/sqrt(3)


        Dk=Dk.reshape((self.DIM-1)**self.NDIM)

        #plot the distribution with von Mises model as a comparison
        fig = plt.figure(facecolor='white', figsize=(12,8))
        plt.xlim(-2,6)
        plt.hist(Dk,1000,histtype='step',normed='True',label=r'Phase increment $D_\mathrm{k}$ for simulated box with bubbles')
        mu=np.mean(Dk)
        x=np.linspace(-2,6,1000)
        y=von_mises(x,mu,1.35)
        plt.plot(x,y,label=r'Von Mises distribution for $\mu= \langle D_\mathrm{k} \rangle$ and $\kappa=1.35$')

        plt.text(4.,0.4,'Statistical properties')
        plt.text(4.,0.38,'Mean: '+str(mu))
        plt.text(4.,0.36,'Median: '+str(np.median(Dk)))
        plt.text(4.,0.34,'Std: '+str(np.std(Dk)))

        plt.title(r'Phase increment distribution for ionised bubbles box with $Q_\mathrm{H_{II}}$ = '+str(self.fillingfraction))
        plt.legend(loc='best',frameon=False)
        plt.savefig('Distribution_phase_increment_bubblebox_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')

        ## ENTROPY
        ## integral of Dk*log(Dk) as a function of phase increment over range 0, 2pi 
        from scipy.special import xlogy
        DlogD=np.multiply(Dk,xlogy(np.sign(Dk),Dk))
        S=-np.sum(DlogD)
        Q=np.log(2*pi)-S #phase structure quantity

        print('Entropy: '+str(S))
        print('Mean increment: '+str(np.mean(Dk)))

        self.phase_k = phase_k 
        self.power_spectrum = PS

    def real_phase_plot(self):

        import matplotlib.pyplot as plt
        import numpy as np

        from mpl_toolkits.axes_grid1 import ImageGrid

        PS, phase_k, Dk = self.phase_analysis(self.box)
        if self.NDIM==3:
            phase=np.fft.irfftn(phase_k,axes=(0,1,2),s=[self.DIM,self.DIM,self.DIM])
        elif self.NDIM==2:
            phase=np.fft.irfftn(phase_k,axes=(0,1),s=[self.DIM,self.DIM])

        fig = plt.figure()
        if (self.NDIM==2):
            plt.imshow(phase)
        elif (self.NDIM==3):
            grid = ImageGrid(fig, 111, # similar to subplot(111)
                        nrows_ncols = (2, 2), # creates 2x2 grid of axes
                        axes_pad=0.12, # pad between axes in inch.
                        label_mode = "L",
                        share_all = "False",
                        cbar_location = "right",
                        cbar_mode="single",
                        cbar_pad = 0.05,
                        cbar_size = "2%"
                        )
            img0=grid[0].imshow(phase[:,100,:])
            img1=grid[1].imshow(phase[:,150,:])
            img2=grid[2].imshow(phase[100,:,:])
            img3=grid[3].imshow(phase[:,:,100])
            grid.cbar_axes[0].colorbar(img0)

            plt.savefig('Phase_realspace_bubblebox_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')

    def power_spectrum_plot(self):

        import numpy as np

        if self.NDIM==3:
            k_x=np.fft.fftfreq(self.DIM) #gives frequency spectrum corresponding to the Fourier decomposition for N1 cells
            a=np.power(k_x,2)[:, np.newaxis] + np.power(k_x,2) 
            b=a[:,:,np.newaxis] + np.power(k_x,2)
            k_norm=np.power(b,1/2)
        elif self.NDIM==2:
            k_x=np.fft.fftfreq(self.DIM) #gives frequency spectrum corresponding to the Fourier decomposition for N1 cells
            a=np.power(k_x,2)[:, np.newaxis] + np.power(k_x,2) 
            k_norm=np.power(a,1/2)

        PS = self.power_spectrum

        nbins=300
        delta_k=k_x[1]
        k_bin=np.logspace(np.log10(delta_k),np.log10(k_norm.max()),nbins)
        count=np.zeros(nbins)
        k_store=np.zeros(nbins)
        P_store=np.zeros(nbins)

        for m in range(0,nbins-1):
            mask= (k_bin[m]<=k_norm) & (k_norm<k_bin[m+1])
            if ~np.any(mask): #cannot compute mean if zero elements
                k_store[m]=0
                P_store[m]=0
                count[m]=0
            else:
                k_store[m]=np.mean(k_norm[mask])
                P_store[m]=np.mean(PS[mask])
                count[m]=k_norm[mask].size

        mask2 = (count!=0) & (k_store!=0)
        k_store=k_store[mask2]
        P_store=P_store[mask2]
        #count=count[mask2]

        #k3PS=np.multiply(np.power(k_store,3),P_store)
 
        plt.figure()
        plt.loglog(k_store/h,P_store)

        plt.xlabel(r'$k$ [$h$ Mpc$^{-1}$]')
        plt.ylabel(r'Matter power spectrum $\mathcal{P}(k)$ [($h$ Mpc$^{-1}$)$^3$]')

        plt.savefig('Power_spectrum_bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')

        return k_store, P_store

    # def phase_increment_spectrum(self):

    #     import numpy as np
    #     import matplotlib.pyplot as plt
    #     plt.ion()

    #     if self.NDIM==3:
    #         k_x=np.fft.fftfreq(self.DIM) #gives frequency spectrum corresponding to the Fourier decomposition for N1 cells
    #         a=np.power(k_x,2)[:, np.newaxis] + np.power(k_x,2) 
    #         b=a[:,:,np.newaxis] + np.power(k_x,2)
    #         k_norm=np.power(b,1/2)
    #     elif self.NDIM==2:
    #         k_x=np.fft.fftfreq(self.DIM) #gives frequency spectrum corresponding to the Fourier decomposition for N1 cells
    #         a=np.power(k_x,2)[:, np.newaxis] + np.power(k_x,2) 
    #         k_norm=np.power(a,1/2)

    #     k_gradx=np.array((k_norm[1:,:,:]+k_norm[:-1,:,:])/2)[:,:-1,:-1]
    #     k_grady=np.array((k_norm[:,1:,:]+k_norm[:,:-1,:])/2)[:-1,:,:-1]
    #     k_gradz=np.array((k_norm[:,:,1:]+k_norm[:,:,:-1])/2)[:-1,:-1,:]
    #     k_norm2=np.power(np.power(k_gradx,2)+np.power(k_grady,2)+np.power(k_gradz,2),1/2)/sqrt(3)

    #     _, phase_k, _ = self.phase_analysis(self.box)

    #     Dx=np.diff(phase_k,axis=0)[:,:-1,:-1] #removes last columns since they won't be used in the overall increment
    #     Dy=np.diff(phase_k,axis=1)[:-1,:,:-1]
    #     Dz=np.diff(phase_k,axis=2)[:-1,:-1,:]
    #     #overall phase increment
    #     Dk=np.power(np.power(Dx,2)+np.power(Dy,2)+np.power(Dz,2),1/2)/sqrt(3)
        
    #     k_bin=np.linspace(k_norm2.min(),k_norm2.max(),250)
    #     count=np.zeros(250)
    #     k_store=np.zeros(250)
    #     Dk_store=np.zeros(250)

    #     #phase_k2=np.where(phase_k<0, phase_k+2*pi, phase_k) #shift spectrum from -pi,pi to 0,2pi
    #     for m in range(0,250-1):
    #         mask= (k_bin[m]<=k_norm2) & (k_norm2<k_bin[m+1])
    #         k_store[m]=np.mean(k_norm2[mask])
    #         Dk_store[m]=np.mean(Dk[mask])
    #         count[m]=int(k_norm2[mask].size)
    #     #avoids zero bins for log scale
    #     mask = (count!=0) & (k_store!=0)
    #     k_store=k_store[mask]
    #     Dk_store=Dk_store[mask]
    #     count=count[mask]

    #     plt.figure(facecolor='white')
    #     plt.plot(k_store,Dk_store)

    #     #plt.xscale('log')
    #     #plt.yscale('log')
    #     plt.xlabel(r'$k$ ')
    #     plt.ylabel(r'$\phi_{k+1} - \phi_k$ ')
    #     plt.title('Phase increment spectrum of simulated box with QHII='+str(self.fillingfraction))

    #     plt.savefig('Phase_increment_spectrum_bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'.png')

    def phase_statistics(self):

        import scipy.stats
        import numpy as np
        import matplotlib.pyplot as plt

        if self.NDIM==2:
            raise Exception("Cannot perform this analysis in 2D")

        plt.ioff()
        _, phase_k, _ = self.phase_analysis(self.box)

        #KURTOSIS
        kurtosis_x=scipy.stats.kurtosis(phase_k,axis=0)
        kurtosis_y=scipy.stats.kurtosis(phase_k,axis=1)
        kurtosis_z=scipy.stats.kurtosis(phase_k,axis=2)

        #VARIANCE
        variance_x=np.var(phase_k,axis=0)
        variance_y=np.var(phase_k,axis=1)
        variance_z=np.var(phase_k,axis=2)

        #ENTROPY

        from scipy.special import xlogy
        plogp=np.multiply(abs(phase_k),xlogy(np.sign(abs(phase_k)),abs(phase_k))) #if sign(phase)=0 ie if phase=0 then output is 0, otherwise log(phi)
        entropy_x=-np.sum(plogp,axis=0)
        entropy_y=-np.sum(plogp,axis=1)
        entropy_z=-np.sum(plogp,axis=2)

        #PLOT

        plt.ion()
        fig, axes= plt.subplots(nrows=4, ncols=self.NDIM, figsize=(11,8))


        kurx = axes[0,0].imshow(kurtosis_x,cmap='Greys')
        fig.colorbar(kurx,ax=axes[0,0])
        axes[0,0].text(0.5, 1.1,'Along x', ha='center',va='bottom',transform=axes[0,0].transAxes,size=12)

        kury = axes[0,1].imshow(kurtosis_y,cmap='Greys')
        fig.colorbar(kury,ax=axes[0,1])
        axes[0,1].text(0.5, 1.1,'Along y', ha='center',va='bottom',transform=axes[0,1].transAxes,size=12)

        kurz = axes[0,2].imshow(kurtosis_z,cmap='Greys')
        fig.colorbar(kurz,ax=axes[0,2])
        axes[0,0].text(0.5, 1.1,'Along z', ha='center',va='bottom',transform=axes[0,2].transAxes,size=12)


        varx = axes[1,0].imshow(variance_x,cmap='Greys')
        fig.colorbar(varx,ax=axes[1,0])

        vary = axes[1,1].imshow(variance_y,cmap='Greys')
        fig.colorbar(vary,ax=axes[1,1])

        varz = axes[1,2].imshow(variance_z,cmap='Greys')
        fig.colorbar(varz,ax=axes[1,2])


        entx = axes[2,0].imshow(entropy_x,cmap='Greys')
        fig.colorbar(entx,ax=axes[2,0])

        enty = axes[2,1].imshow(entropy_y,cmap='Greys')
        fig.colorbar(enty,ax=axes[2,1])

        entz = axes[2,2].imshow(entropy_z,cmap='Greys')
        fig.colorbar(entz,ax=axes[2,2])


        skewx = axes[3,0].imshow(skew_x,cmap='Greys')
        fig.colorbar(entx,ax=axes[3,0])

        skewy = axes[3,1].imshow(skew_y,cmap='Greys')
        fig.colorbar(enty,ax=axes[3,1])

        skewz = axes[3,2].imshow(skew_z,cmap='Greys')
        fig.colorbar(entz,ax=axes[3,2])


        axes[0,0].text(-0.8, 0.5,'Kurtosis', ha='center',va='bottom',transform=axes[0,0].transAxes,size=12)
        axes[1,0].text(-0.8, 0.5,'Variance', ha='center',va='bottom',transform=axes[1,0].transAxes,size=12)
        axes[2,0].text(-0.8, 0.5,'Entropy', ha='center',va='bottom',transform=axes[2,0].transAxes,size=12)
        axes[3,0].text(-0.8, 0.5,'Skewness', ha='center',va='bottom',transform=axes[3,0].transAxes,size=12)

        plt.suptitle('Phase statistics in Fourier space for simulated box with QHII='+str(self.fillingfraction))
        #plt.tight_layout()
        plt.show()

        plt.savefig('Statistics_bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D.png')


    def write_overdensity(self):

        meand=self.box.mean()
        overdensity=(self.box-meand)/meand

        #Fourier transform
        rfilechain='Overdensity_bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D_realpart.txt'
        ifilechain='Overdensity_bubble_box_QHII='+str(self.fillingfraction)+'_nooverlap='+str(self.nooverlap)+'_dim'+str(self.DIM)+'_'+str(self.NDIM)+'D_imagpart.txt'
        if self.NDIM==3:
            overdensity_k=np.fft.fftn(overdensity,axes=(0,1,2))
            out_real=np.real(overdensity_k[:,:,0])
            for i in range(1,self.DIM):
                out_real=np.r_[out_real,np.real(overdensity_k[:,:,i])]
            np.savetxt(rfilechain,out_real,delimiter=' ',fmt='%-6.4f')
            out_imag=np.imag(overdensity_k[:,:,0])
            for i in range(1,self.DIM):
                out_imag=np.r_[out_imag,np.imag(overdensity_k[:,:,i])]
            np.savetxt(ifilechain,out_imag,delimiter=' ',fmt='%-6.4f')

        elif self.NDIM==2:
            overdensity_k=np.fft.fftn(overdensity,axes=(0,1))
            out_real=np.real(overdensity_k)
            out_imag=np.imag(overdensity_k)
            np.savetxt(rfilechain,out_real,delimiter=' ',fmt='%-6.4f')
            np.savetxt(ifilechain,out_imag,delimiter=' ',fmt='%-6.4f')
      

def von_mises(x,mu,K):
    """ Outputs von Mises distribution function for distribution x, mean mu and dispersion K"""
    from scipy import special
    from math import pi,sqrt
    out=np.zeros(len(x))
    for i in range(0,len(x)):
        out[i]= ((2*pi*special.iv(0,K))**(-1))* np.exp(K*np.cos(x[i]-mu))
    return out