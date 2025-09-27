using Plots

# parameters
H = 0.5

"""
Objective : Define  the  Mesh Geomety
    STEP 1 : first create a differential 2D Control Volume with 
    length along x = delta_x
    length along y = delta_y
    n = required number of such differential control volumes
    required to construct full CV

    STEP 2 : Compute the computational nodes for each
             differential control volume
"""

# STEP 1 ##### 

delta_x = 0.01 # change to get variable geometry 
delta_y = 0.01

# set n to be even number 
n = 10 

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


"""
Objectives : COMPUTATIONS
        STEP 1 : identify the boundary nodes and apply the
                 boundary conditions  

                 Boundary 1 : T_1 = 15
                 Boundary 2 : T_2 = 10 
                 Boundary 3 : T_3 = 5(1-y/H) + 15 * sin(pi*y/H)

        STEP 2 : write the equation for boundary 4

        STEP 3 : write the equation for internal nodes

        STEP 4 : setup the conditons for tolerance
                 Approach :
                 1. pick n random nodes from grid
                 2. save temperature before each iteration
                 3. find the temperature after iteratrion 
                 4. diff = after_iteration_temperature - before_temperature
                 5. elementwise square each difference
                    diff.^2
                 6. max(diff.^2) < tolerance   

        STEP 5 : perform the computations 
"""


# STEP 1 #####  

# identify the boundary nodes
T_4 = nodes[1   ,   :]
T_2 = nodes[n-1 ,   :]
T_3 = nodes[:   , n-1]
T_1 = nodes[:   ,   1]

# apply the boundary conditions to the suitable nodes
# T_1 set to 15 <<<
for i = 1:length(T_4)
    T_1[i][3] = 15
end
# >>>

# T_2 is function of y >>>
function T2(y)
    global H
    5(1-y/H) + 15 * sin(pi*y/H)
end

for i = 1:length(T_2)
    T_2[i][3] = T2(T_2[i][2])
end
# <<<

# T_3 set to 10 >>> 
for i = 1:length(T_3)
    T_3[i][3] = 10
end
# <<< 


# STEP 2 ##### Boundary 4

# equations for Boundary 4 <<<
function B4Temp(node_element, T_E , T_N)
    global delta_x , delta_y
    y = node_element[2]
    k  = 16(y/H + 1)

    a_E  = 1 / delta_x
    a_N  = 1 / delta_y
    (a_E* T_E + a_N * T_N - 5000/k)/(a_E + a_N)
end
# <<<


# STEP 3 ##### equation for internal nodes 

# equations for internal nodes >>>
function internalTemp(node_element , T_W, T_E , T_S , T_N)
    # use global parameters 
    global H , delta_x, delta_y 

    # TODO stretch mesh 
    # calculate node geometry
    y = node_element[2]
    A = delta_x * delta_y
    k  = 16(y/H + 1)

    # TODO change for stretch mesh

    # TODO clarification needed on S_u and S_p
    S_p = 0
    S_u = -1.5 * A

    a_W = k * A/delta_x 
    a_E = k * A/delta_x 
    a_S = k * A/delta_y 
    a_N = k * A/delta_y 
    a_P = a_W + a_E + a_S + a_N - S_p

    # calculate at T_P
    T_P = (a_W* T_W + a_E * T_E + a_S * T_S +
        a_N * T_N + S_u) / a_P
end
# <<<

# STEP 4 ##### test for convergence 

# set the tolerance here
tolerance = 0.001



# function to return mid node for plotting
# no of loops required for achieving required tolerance
function mid_node(nodes)
    # nodes : grid of nodes
    # returns mid node in the grid
    nodes[Int(n/2), Int(n/2)][3]
end

# STEP 5 ##### Computations
function calculate_temp!(nodes)
    for i = 1:n-2 , j = 2:n-2
        current_node = nodes[i,j] 

        if i == 1 # check for Boundary 4 
            E_node = nodes[i+1 , j  ]
            N_node = nodes[i   , j+1]
            T_E = E_node[3]
            T_N = N_node[3]
            new_temp = B4Temp(current_node, T_E , T_N)

        else # internal nodes 
            current_node = nodes[i,j] 
            W_node = nodes[i-1 , j  ]
            S_node = nodes[i   , j-1]
            E_node = nodes[i+1 , j  ]
            N_node = nodes[i   , j+1]

            T_W = W_node[3] 
            T_S = S_node[3]
            T_E = E_node[3]
            T_N = N_node[3]

            new_temp = internalTemp(current_node , T_W,
                                        T_E , T_S , T_N)
        end

        # update the current node
        current_node[3] = new_temp 
    end
end


loop_counter = 0 
print(" >>> COMPUTING >>>" )

# returns n  random number  between 1 and n^2
# n random nodes are selected from [1, n^2] nodes
# for purpose of testing convergence 
test_nodes_id = rand(1:(n-1)^2 , n)
init_nodes_temp = Array{Float64}(undef,n)  
updated_nodes_temp = Array{Float64}(undef,n)  

while true  
    global init_nodes_temp, updated_nodes_temp
    ## before each iterations 
    for i = 1:length(test_nodes_id) 
        init_nodes_temp[i] = nodes[test_nodes_id[i]][3]
    end

    ###################
    calculate_temp!(nodes)
    ###################

    # after iteration temp is updated
    for i = 1:length(test_nodes_id) 
        updated_nodes_temp[i] = nodes[test_nodes_id[i]][3]
    end

    test = abs(maximum((updated_nodes_temp -
        init_nodes_temp).^2)) < tolerance

    @show test

    if  test
        break
    end

    init_nodes_temp = updated_nodes_temp[:]
    global loop_counter += 1
    print(" no of loop  " ,loop_counter, "\n " )
end

print(" <<< COMPLETED <<< ")

# <<< Plotting <<<<
using Plots

function professional_contour(data)
    # Get unique coordinates
    x_vals = sort(unique([p[1] for p in data]))
    y_vals = sort(unique([p[2] for p in data]))
    
    # Create temperature matrix
    temp_matrix = zeros(length(y_vals), length(x_vals))
    for point in data
        i = findfirst(y -> y ≈ point[2], y_vals)
        j = findfirst(x -> x ≈ point[1], x_vals)
        temp_matrix[i, j] = point[3]
    end
    
    # Professional contour plot
    c = contourf(x_vals, y_vals, temp_matrix,
            title="Temperature Distribution",
            xlabel="X Coordinate (m)", 
            ylabel="Y Coordinate (m)",
            colorbar_title="Temperature (°C)",
            aspect_ratio=:equal,
            levels=25,
            color=:hot,
            size=(700, 600),
            dpi=300,
            linewidth=0,
            fontfamily="Computer Modern")
    display(c)
end

# Usage:
professional_contour(nodes)

# write the data to csv file for analysing
function export_csv(data, filename="temp.csv")
    open(filename, "w") do f
        println(f, "X,Y,Temperature")
        for point in data
            println(f, "$(point[1]),$(point[2]),$(point[3])")
        end
    end
end

# Usage:
export_csv(nodes)
