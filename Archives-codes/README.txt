#########################################################################################
Here is a very small instruction to use the program via USER.cpp. 
For a more detailed explanation of the code, please refer to the comments in project.cpp.
#########################################################################################

System initialization:
*Modifies the total number of particles by changing the value of "Nx".
*Modifies the time step by changing the value of "h".
*Modifies the critical distance by changing the of "rc" and the additionnal distance by changing delta_r (for integration_neighbours onyl).
*Modifies the initial velocity norm by changing the value of "v0".
*Modifies the size of the ferrofluid droplet by changing the value of "size": the radius of the droplet will be sqrt(size*L)

Options:
-sys.mag: if true, the magnetic interaction will be taken into account, if false, they don't.
-sys.cutoff: It is the maximal value of the force in order to delay the velocity divergence.
             Take a large value of sys.cutoff (eg. 1E30) if you want no cutoff.

******

Creation of the files:
The files necessary for the display are created. Other files used to record various information can be added by using the same synthaxe.

******

Simulation and file filling:
*Modifies the number of main loops by changing the value of "N_it1".
/!\ N_it1 needs to be greater than 20 because of the terminal display of the pourcentage.
*Modifies the number of iteration in each main loop by changing the value of "N_it2".
/!\ The recording of trajectories is made in each main loop. Hence, if you perform long simulation, take care to take a small N_it1 and big N_it2.
    Otherwise you will have a lot of data and a very long and bulky display.

Rmk: the total number of iteration is N_it1 x N_it2

******

Graphic interface:
This part takes care of the display in the terminal. There is actually nothing to modify.

******

/!\ At the end of USER.cpp, a Python program is automatically called to create and save a .GIF animation. 
    Comment this line if you do not want an animation.
