using DelimitedFiles
total_set = 8
histo_data_vec = Array{Any,1}(undef,total_set)
histo_data_vec[1] = vec(readdlm("data/Signal_data/Dataset#1/20170511_YTL047A1_GAL_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G1cells.txt", ',', Float64))
histo_data_vec[2] = vec(readdlm("data/Signal_data/Dataset#2/20170511_YTL047A2_GAL_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G1cells.txt", ',', Float64))
histo_data_vec[3] = vec(readdlm("data/Signal_data/Dataset#3/20170824_YTL047A1_gal_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G1cells.txt", ',', Float64))
histo_data_vec[4] = vec(readdlm("data/Signal_data/Dataset#4/20170824_YTL047A2_gal_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G1cells.txt", ',', Float64))

histo_data_vec[5] = vec(readdlm("data/Signal_data/Dataset#1/20170511_YTL047A1_GAL_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G2cells.txt", ',', Float64))
histo_data_vec[6] = vec(readdlm("data/Signal_data/Dataset#2/20170511_YTL047A2_GAL_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G2cells.txt", ',', Float64))
histo_data_vec[7] = vec(readdlm("data/Signal_data/Dataset#3/20170824_YTL047A1_gal_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G2cells.txt", ',', Float64))
histo_data_vec[8] = vec(readdlm("data/Signal_data/Dataset#4/20170824_YTL047A2_gal_PP7-cy3_1_intensity_brightest_txpnsite_cy3_G2cells.txt", ',', Float64))

histo_data_vec_merged = Array{Any,1}(undef,4)
histo_data_vec_merged[1] = readdlm("data/DAPI_content/20170511_YTL047A1_GAL_PP7-cy3_1_dapiContent_nascentTSintensity_allcells.txt")[:,2]
histo_data_vec_merged[2] = readdlm("data/DAPI_content/20170511_YTL047A2_GAL_PP7-cy3_1_dapiContent_nascentTSintensity_allcells.txt")[:,2]
histo_data_vec_merged[3] = readdlm("data/DAPI_content/20170824_YTL047A1_gal_PP7-cy3_1_dapiContent_nascentTSintensity_allcells.txt")[:,2]
histo_data_vec_merged[4] = readdlm("data/DAPI_content/20170824_YTL047A2_gal_PP7-cy3_1_dapiContent_nascentTSintensity_allcells.txt")[:,2]


normalizing_valuev1 = vec(readdlm("data/Signal_data/cyto_normalization_values.txt",Float64))[1:4]
normalizing_value = repeat(normalizing_valuev1,2);
histo_data_vec_normalized = histo_data_vec./normalizing_value;


# length.(histo_data_vec_normalized)

for i in 1:4
    filter!(x->x.!="N/A",histo_data_vec_merged[i])
end

normalizing_valuev1 = vec(readdlm("data/Signal_data/cyto_normalization_values.txt",Float64))[1:4];
normalizing_value = repeat(normalizing_valuev1,1);
histo_data_vec_merged_normalized = histo_data_vec_merged./normalizing_value;