
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
                 and

                 prepare the required helper function for computation

        STEP 5 : perform the computations 
"""

include("mesh_geometry.jl")

# STEP 1 ##### Boundary nodes

# identify the boundary nodes and apply Boundary condition
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


# STEP 2 ##### equation for Boundary 4

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

# STEP 4 ##### test for convergence  and required helper functions

function calculate_temp!(nodes)
    """
    Calculates Temperature at
    1. Boundary 4
    2. All internal nodes
    mutates the nodes array 
    """
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

# n random nodes are selected from [1, n-1^2] nodes
# for purpose of testing convergence 
test_nodes_id = rand(1:(n-1)^2 , n)
init_nodes_temp = Array{Float64}(undef,n)  
updated_nodes_temp = Array{Float64}(undef,n)  

function experiment!(tolerance)
    """
    # Approach
     Approach is to set a tolerance value and
    continue the iteration untill the temperature at nodes
    start converging ,

    - After each iteration test at random nodes if the difference
    of convergence is with in an acceptable range defined by
    parameter tolerance  

    - uses calculate_temp! which mutates the original nodes 

    # Arguments
    -`tolerance::Float`: set up the tolerance

    # Returns and Operations

    - call calculate_temp! function which mutates the nodes array
    with new temperatures and continues till the value start to converge
    at randomly selected n nodes

    -`loop_counter::Int` : return the number of loops required

    """
    # setup the loop counter for plotting
    loop_counter = 0 

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
        loop_counter += 1
    end
    @show loop_counter
    loop_counter
end

# STEP 5 ##### Computations

# # set the tolerance here 
# tolerance = 0.001
# # start the computations
# experiment(tolerance)
