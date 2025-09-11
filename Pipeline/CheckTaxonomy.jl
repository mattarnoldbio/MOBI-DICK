#%%
using DataFrames
using CSV
using NCBITaxonomy
using AbstractTrees
using StatsBase

#%%

viral_hits = CSV.read("/Users/arno01m/OneDrive - University of Glasgow/MarineMammalVirome/SRA_mining/Results/second_run/vertebrate_viral_hits.csv", DataFrame)


# is_caudoviricetes = []

# for (i, tx) in enumerate(viral_hits.species)
#     println("Processing hit $i")
#     if viral_hits.genus[i] == "none" || viral_hits.family[i] == "none"
#         print("a")
#         lin  = taxon(tx, strict=false)
#         print("a_")
#         lin = lineage(lin)
#         print("b")
#         if  ncbi"Caudoviricetes" in filter(t -> rank(t) == :class , lin)
#             print("c")
#             push!(is_caudoviricetes, true)
#         else
#             print("d")
#             push!(is_caudoviricetes, false)
#         end
#     else
#         push!(is_caudoviricetes, false)
#     end
# end

# viral_hits.is_caudoviricetes = is_caudoviricetes


#%%

viral_hits = CSV.read("/Users/arno01m/OneDrive - University of Glasgow/MarineMammalVirome/SRA_mining/Results/second_run/vertebrate_viral_hits_updated.csv", DataFrame)

#%%

# Drop rows where is_caudoviricetes is true
viral_hits = viral_hits[.!viral_hits.is_caudoviricetes, :]

#%%

CSV.write("/Users/arno01m/OneDrive - University of Glasgow/MarineMammalVirome/SRA_mining/Results/second_run/vertebrate_viral_hits_updated.csv", viral_hits)
