Masanori Shimono, 20/05  wrote
Motoki Kajiwara,  20/11  translated

 ---------------  
# Example
  You can change data_index　to any value.
  data_index = 1; run_and_save0
 
 ---------------  
 # Program dependencies
  run_and_save0.m is the main code to be ran, which cites two other m-files, model_ver2020_simple_ver7.m and Post2Conmat_rev.m.

---------------  
# Steps of simulation

   At the beginning of the simulation, the movement is slow and the activity status will seem to be strange. However the figure 
   keeps appearing and gradually shifts to the activity closer to the behavior of real neuronal networks. Finally, a folder called 
   "data" is produced saving the results of recording the second half of the simulation. Also the data of the connection matrix 
   called conmat will be saved.
  
---------------  
# Figure monitoring activities

  In the upper left of the figure, the horizontal axis is the time and the vertical axis is the time series of the cells' index. 
  A blue dot is pointed where you are active In the upper right of the figure, cells that are active in a certain time cross 
  section are represented by yellow markers. The connection matrix (binary) is displayed at the bottom left of the figure. 
  The bottom right panel shows distributions of connectivity weights.


