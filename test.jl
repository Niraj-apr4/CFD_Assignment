using Plots

# mesh geometry 
H = 0.5

# set n
n = 10
delta_x = 0.001 
delta_y = 0.001

points = Array{Array{Float64}}(undef,n,n)
start_point_x = 0
start_point_y = 0

for i = 1:n 
    for j = 1:n 
        global start_point_x
        global start_point_y
        points[i,j] = [start_point_x start_point_y 0]
        global start_point_y += delta_y 
    end
    global start_point_y = 0 
    global start_point_x += delta_x 
end

# write down computational nodes for each cv

# # for plotting  change it 
# x_cords = [p[1] for p in points]
# y_cords = [p[2] for p in points]
# plot(x_cords, y_cords , marker = :line)
# # for plotting


nodes = Array{Array{Float64}}(undef,n-1,n-1)
for i = 1:n-1 ,j = 1:n-1
    nodes[i,j] = [points[i,j][1] + delta_x/2 , points[i,j][2] +
        delta_y/2 , 0]
end

# boundary conditions
# identify the boundary nodes
T_4 = nodes[1   ,   :]
T_2 = nodes[n-1 ,   :]
T_3 = nodes[:   , n-1]
T_1 = nodes[:   ,   1]

# apply the boundary conditions to nodes

# T_1 set to 15
for i = 1:length(T_4)
    T_1[i][3] = 15
end

# T_3 set to 10
for i = 1:length(T_3)
    T_3[i][3] = 10
end


# T_2 is function of y
function T2(y)
    global H
    5(1-y/H) + 15 * sin(pi*y/H)
end

for i = 1:length(T_2)
    T_2[i][3] = T2(T_2[i][2])
end

# TODO change to  appropriate  Boundary condition 
# T_4 set to 11 (for testing)
function B4Temp(node_element, T_E , T_N)
    global delta_x , delta_y
    y = node_element[2]
    k  = 16(y/H + 1)

    a_E  = 1 / delta_x
    a_N  = 1 / delta_y
    (a_E* T_E + a_N * T_N - 5000/k)/(a_E + a_N)
end

for i = 1:length(T_4)
    T_4[i][3] = 11 
end

# internal nodes
# calculation at internal nodes
function internalTemp(node_element , T_W, T_E , T_S , T_N)
    # use global parameters 
    global H , delta_x, delta_y , k

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

# test
# till the solution dont converge iterate
# TODO find out how to test  convergence 
# pick the mid nodes and see if it converges 

# function to return mid node
function mid_node(nodes)
    # nodes : grid of nodes
    # returns mid node in the grid
    nodes[Int(n/2), Int(n/2)][3]
end


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
while true  
    # capture the mid value and save it 
    init_mid_temp = mid_node(nodes)

    calculate_temp!(nodes)

    # after iteration temp is updated
    # recalculate the mid value 
    updated_mid_temp = mid_node(nodes)

    # compare 
    if abs(init_mid_temp - updated_mid_temp) < 0.00001 
        break
    end
    global loop_counter += 1
    print(" no of loop  " ,loop_counter, "\n " )
end

print(" <<< COMPLETED <<< ")

# <<< Plotting <<<<

using Plots
# # Direct plotting:
# x = [p[1] for p in nodes]
# y = [p[2] for p in nodes] 
# t = [p[3] for p in nodes]
# s = scatter(x, y, zcolor=t,aspect_ratio=:equal,
#             title="Temperature")
# display(s)
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

# # write the data to csv file for analysing
# function export_csv(data, filename="temp.csv")
#     open(filename, "w") do f
#         println(f, "X,Y,Temperature")
#         for point in data
#             println(f, "$(point[1]),$(point[2]),$(point[3])")
#         end
#     end
# end

# # Usage:
# export_csv(nodes)
