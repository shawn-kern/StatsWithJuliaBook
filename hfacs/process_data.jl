using LinearAlgebra, DataStructures, DataFrames, CSV, Query
using Distributions, CategoricalArrays, FreqTables, HypothesisTests

raw_data = CSV.read("hfacs/HFACS_FY19-22_5YR.csv", DataFrame, 
types=Dict(:Type=>String, :NanoCode=>Symbol))

dropmissing!(raw_data, [:Type, :NanoCode], disallowmissing=true)

raw_data.Type = categorical(raw_data.Type, ordered=true)
levels!(raw_data.Type, ["Non-8621 Reportable", "Close Call",    
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

nanocode_index = sort(collect(union(nanocodes.codes...)))

nanocode_names = Dict()
for row in eachrow(unique(raw_data[:, [:NanoCode, :NanoCode_Name]]))
    nanocode_names[row.NanoCode] = row.NanoCode_Name
end

frequency=Dict()
for code in nanocode_index
    frequency[code, "All"] = sum(hfacs_data[:, code])
    frequency[code, "Armstrong"] = sum(hfacs_data[hfacs_data.Center .== "Armstrong", code])
end

codes = Array{Bool, 2}(undef, nrow(nanocodes), length(nanocode_index))
for (i, code) in enumerate(nanocode_index)
    codes[:, i] = code .∈ nanocodes.codes
end
hfacs_data = hcat(hfacs_data, DataFrame(codes, collect(nanocode_index)))

# analyses
for level in levels(hfacs_data.Type)
    println("Mishap totals for $level: all centers $(sum(hfacs_data.Type .== level)), Armstrong $(sum(hfacs_data.Type[hfacs_data.Center .== "Armstrong"] .== level))")
end

X2_pvalues = [ pvalue(ChisqTest(freqtable(hfacs_data[:, i], hfacs_data[:,j])), tail=:right)
     for i in nanocode_index, j in nanocode_index ]

function categorical_correlation(index1, index2)
    table = freqtable(hfacs_data[:, index1], hfacs_data[:, index2])
    nr, nc = size(table)
    X2=0
    row_margins = hcat(table, sum(table, dims=2))
    margins = vcat(row_margins, sum(row_margins, dims=1))
    for i=1:nr
        for j=1:nc
            expected = margins[i,nc+1]*margins[nr+1,j]/margins[nr+1,nc+1]
            X2 += (margins[i,j]-expected)^2/expected
        end
    end
    V = sqrt(X2/margins[nr+1,nc+1])
    pvalue = ccdf(Chisq((nr-1)*(nc-1)),X2)
    return V, X2, pvalue
end

# find correlated pairs of nanocodes that occur more than once
function find_index_pairs(; α=0.05, center="All", min_n=1)
    results = DataFrame(Index1=Symbol[], Index2=Symbol[], V=Float64[], pvalue=Float64[], N1=Int[], N2=Int[])
    for i in nanocode_index
        for j in nanocode_index
            if i ≥ j continue end
            V, X2, pvalue = categorical_correlation(i, j)
            if pvalue < α && (frequency[i, center] > min_n && frequency[j, center] > min_n)
                push!(results, (i, j, V, pvalue, frequency[i, center], frequency[j, center]))
            end
        end
    end
    return sort(results, [:pvalue])
end

# find nanocodes correlated with specific index
function find_index_match(index=:Type; α=0.05, center="All", min_n=1)
    results = DataFrame(Index1=Symbol[], Index2=Symbol[], V=Float64[], pvalue=Float64[], N=Int[])
    for j in nanocode_index
        if j == index continue end
        V, X2, pvalue = categorical_correlation(index, j)
        if pvalue < α && frequency[j, center] > min_n
            push!(results, (index, j, V, pvalue, frequency[j, center]))
        end
    end
    return sort(results, [:pvalue])
end
