
"""
mesh_geometry.jl  

Objective : Define  the  Mesh Geomety
    STEP 1 : first create a differential 2D Control Volume with 
    length along x = delta_x
    length along y = delta_y
    n = required number of such differential control volumes
    required to construct full CV

    STEP 2 : Compute the computational nodes for each
             differential control volume
"""

# parameters
H = 0.5

# STEP 1 ##### 

delta_x = 1.0 # change to get variable geometry 
delta_y = 1.0

# set n to be even number 
n = 30 

# initialize An undefined array called "points"
# with size  nxn 
# and its element is itself an array of type {Array{Float64}}
# with following characteristics
#      element[1] = x_cordinate
#      element[2] = y_cordinate
#      element[3] = Temperature (initially set to 0)
# 

points = Array{Array{Float64}}(undef,n,n)

# create the CV grid >>> 
# starting point of coordintaes
start_point_x = 0
start_point_y = 0

for i = 1:n 
    for j = 1:n 
        global start_point_x
        global start_point_y
        # set the temp initially to zero
        points[i,j] = [start_point_x start_point_y 0]
        global start_point_y += delta_y 
    end
    #  reset the start_point_y to 0 for new x cordinates 
    global start_point_y = 0 

    # increase the start_point_x by delta_x
    global start_point_x += delta_x 
end
# <<<

# now the points array consists of CV grids with initial temp
# assigned to be  0 

# STEP 2 ##### 

# Compute  computational nodes for each cv >>>
nodes = Array{Array{Float64}}(undef,n-1,n-1)
for i = 1:n-1 ,j = 1:n-1
    nodes[i,j] = [points[i,j][1] + delta_x/2 , points[i,j][2] +
        delta_y/2 , 0]
end
# <<<
