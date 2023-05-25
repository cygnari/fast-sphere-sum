import numpy as np
import netCDF4 as nc

output_data = nc.Dataset('vr45to5_parallel_output.nc',"r")
ssh = np.asarray(output_data["ssh"]) # Sea-surface height, eta
ssh_sal = np.asarray(output_data["ssh_sal"]) # SSH perturbation due to SAL, eta_SAL

print(np.shape(ssh))

print(np.sum(np.sum(ssh)))
print(np.sum(np.sum(ssh_sal)))
print(np.sum(np.sum(ssh==0)))
print(np.sum(np.sum(ssh_sal==0)))

# Data is output at the end of each time step
# This means that the ssh from one time step should
# give the ssh_sal from the subsequent time step
# Both Arrays have a shape of (4, 115484)
# So, e.g.,  doing the SAL calculation on ssh from (0,:) 
# should lead to the ssh_sal results from (1,:)
