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
