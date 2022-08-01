using DelimitedFiles
total_set = 8
histo_data_vec = Array{Any,1}(undef,total_set)
path = "data/PP7-GAL10_CellProfiler"

histo_data_vec[1] = vec(readdlm("$path/Dataset#1/20170511_YTL047A1_GAL_PP7-cy3_1_cellular_count_cy3_G1cells.txt", ',', Float64))
histo_data_vec[2] = vec(readdlm("$path/Dataset#2/20170511_YTL047A2_GAL_PP7-cy3_1_cellular_count_cy3_G1cells.txt", ',', Float64))
histo_data_vec[3] = vec(readdlm("$path/Dataset#3/20170824_YTL047A1_gal_PP7-cy3_1_cellular_count_cy3_G1cells.txt", ',', Float64))
histo_data_vec[4] = vec(readdlm("$path/Dataset#4/20170824_YTL047A2_gal_PP7-cy3_1_cellular_count_cy3_G1cells.txt", ',', Float64))

histo_data_vec[5] = vec(readdlm("$path/Dataset#1/20170511_YTL047A1_GAL_PP7-cy3_1_cellular_count_cy3_G2cells.txt", ',', Float64))
histo_data_vec[6] = vec(readdlm("$path/Dataset#2/20170511_YTL047A2_GAL_PP7-cy3_1_cellular_count_cy3_G2cells.txt", ',', Float64))
histo_data_vec[7] = vec(readdlm("$path/Dataset#3/20170824_YTL047A1_gal_PP7-cy3_1_cellular_count_cy3_G2cells.txt", ',', Float64))
histo_data_vec[8] = vec(readdlm("$path/Dataset#4/20170824_YTL047A2_gal_PP7-cy3_1_cellular_count_cy3_G2cells.txt", ',', Float64))

histo_data_vec
histo_data_vec_merged = Array{Any,1}(undef,4)
histo_data_vec_merged[1] = vec(readdlm("$path/Dataset#1/20170511_YTL047A1_GAL_PP7-cy3_1_cellular_count_cy3_AllCells.txt", ',', Float64))
histo_data_vec_merged[2] = vec(readdlm("$path/Dataset#2/20170511_YTL047A2_GAL_PP7-cy3_1_cellular_count_cy3_AllCells.txt", ',', Float64))
histo_data_vec_merged[3] = vec(readdlm("$path/Dataset#3/20170824_YTL047A1_gal_PP7-cy3_1_cellular_count_cy3_AllCells.txt", ',', Float64))
histo_data_vec_merged[4] = vec(readdlm("$path/Dataset#4/20170824_YTL047A2_gal_PP7-cy3_1_cellular_count_cy3_AllCells.txt", ',', Float64))

[length.(histo_data_vec);length.(histo_data_vec_merged)]