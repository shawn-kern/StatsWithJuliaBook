using DataFrames, CSV, Query, CategoricalArrays, DataStructures, FreqTables, LinearAlgebra

raw_data = CSV.read("hfacs/HFACS_FY19-22_5YR.csv", DataFrame, 
types=Dict(:Type=>String, :NanoCode=>Symbol))

dropmissing!(raw_data, [:Type, :NanoCode], disallowmissing=true)

raw_data.Type = categorical(raw_data.Type, ordered=true)
levels!(hfacs_data.Type, ["Non-8621 Reportable", "Close Call",    
    "Mishap - Type D", "Mishap - Type C", "Mishap - Type B", "Mishap - Type A"])

hfacs_data = @from row in raw_data begin
  @select { row.ID, Center=split(row.Center)[1], row.Type, 
#    Injury=row.Injury=="Yes" ? 1 : 0, Damage=row.Damage=="Yes" ? 1 : 0, 
    Injury=row.Injury=="Yes", Damage=row.Damage=="Yes", 
    row.NanoCode, Days=row.OSHA_Days_Away,
    Cost=row.Final_Cost }
  @collect DataFrame
end

nanocodes = @from row in hfacs_data begin
    @group row by row.ID into grp
    @select { ID=key(grp), codes=Set(grp.NanoCode) }
    @collect DataFrame
end

costs = @from row in hfacs_data begin
    @group row by row.ID into grp
    @select { ID=key(grp), Days=maximum(grp.Days), Cost=maximum(grp.Cost) }
    @collect DataFrame
end

hfacs_data = hfacs_data[:, Not([:NanoCode, :Days, :Cost])]
unique!(hfacs_data)
hfacs_data = innerjoin(hfacs_data, costs, on=:ID)

#nanocode_index = SortedSet(union(nanocodes.codes...))
nanocode_index = sort(collect(union(nanocodes.codes...)))

codes = Array{Bool, 2}(undef, nrow(nanocodes), length(nanocode_index))
for (i, code) in enumerate(nanocode_index)
    codes[:, i] = code .∈ nanocodes.codes
end
hfacs_data = hcat(hfacs_data, DataFrame(codes, collect(nanocode_index)))

# simple analyses
for level in levels(hfacs_data.Type)
    println("Mishap totals for $level: all centers $(sum(hfacs_data.Type .== level)), Armstrong $(sum(hfacs_data.Type[hfacs_data.Center .== "Armstrong"] .== level))")
end

for code in nanocode_index
    println(code, " ", sum(hfacs_data[:, code]), " ", sum(hfacs_data[hfacs_data.Center .== "Armstrong", code]))
end

full_matrix = [ dot(hfacs_data[:, nanocode_index[i]], hfacs_data[:, nanocode_index[j]]) 
    for i in 1:length(nanocode_index), j in 1:length(nanocode_index) ]

Armstrong_matrix = [ dot(hfacs_data[hfacs_data.Center .== "Armstrong", 
    nanocode_index[i]], hfacs_data[hfacs_data.Center .== "Armstrong", nanocode_index[j]]) 
    for i in 1:length(nanocode_index), j in 1:length(nanocode_index) ]

χ2_pvalues = [ pvalue(ChisqTest(freqtable(hfacs_data[:, i], hfacs_data[:,j])), tail=:right)
     for i in nanocode_index, j in nanocode_index ]
V = χ2

for i in 1:length(nanocode_index)
    for j in 1:length(nanocode_index)
        if i==j continue end
        if χ2_pvalues[i, j] < 0.001
            println("$(nanocode_index[i]) & $(nanocode_index[j]) have $(χ2_pvalues[i, j]).")
        end
    end
end

function categorical_correlation(index1, index2)
    table = freqtable(hfacs_data[:, index1], hfacs_data[:, index2])
    X2 = 0
    row_margins = hcat(table, sum(table, dims=2))
    margins = vcat(row_margins, sum(row_margins, dims=1))
    for i=1:2
        for j=1:2
            expected = margins[i,3]*margins[3,j]/margins[3,3]
            X2 += (margins[i,j]-expected)^2/expected
        end
    end
    V = sqrt(X2/margins[3,3])
    pvalue = ccdf(Chisq(1),X2)
    return V, X2, pvalue
end