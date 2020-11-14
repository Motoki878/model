Masanori Shimono, 20/05/29  wrote
Motoki Kajiwara,  20/11/06  translated

 ---------------  
# Example
  You can change data_index　to any value.
  data_index = 1; run_and_save0
 
 ---------------  
 # Program dependencies
  run_and_save0.m is the main code to be ran, which cites two other m-files, model_ver2020_simple_ver7.m and Post2Conmat_rev.m.

---------------  
# Steps of simulation

  At the beginning, the movement is slow and the activity status is strange.
  However the figure keeps appearing and  gradually shifts to the activity closer to reality.
  Finally, a folder called "data" is created to save the results of recording the second half of the simulation.
  Also the data of the connection matrix called conmat is saved.
  
---------------  
# Figure

  In the upper left of the figure, the horizontal axis is the time and the vertical axis is the time series of the cell index.
  A blue dot is pointed where you are active
  In the upper right of the figure, cells that are active in a certain time cross section are represented by yellow markers.
  The connection matrix (binary) is displayed at the bottom left of the figure.
  At the bottom right of the figure, you don't have to worry about it now.

