using LinearAlgebra, DataStructures, DataFrames, CSV, Query
using Distributions, CategoricalArrays, FreqTables, HypothesisTests

raw_data = CSV.read("hfacs/HFACS_ALL_5YR.csv", DataFrame, 
types=Dict(:Type=>String, :NanoCode=>Symbol))

dropmissing!(raw_data, [:Type, :NanoCode], disallowmissing=true)

raw_data.Type = categorical(raw_data.Type, ordered=true)
levels!(raw_data.Type, ["Non-8621 Reportable", "Close Call",    
    "Mishap - Type D", "Mishap - Type C", "Mishap - Type B", "Mishap - Type A"])

base_data = @from row in raw_data begin
    @where row.Type != "Non-8621 Reportable"
    @select { row.ID, Center=split(row.Center)[1], row.Type, 
    Injury=row.Injury=="Yes", Damage=row.Damage=="Yes", 
    row.NanoCode, NanoCode_Reason=row.HFACS_Rationale, Days=row.OSHA_Days_Away, 
    Cost=row.Final_Cost, Desc=row.Description }
  @collect DataFrame
end
droplevels!(base_data.Type)

nanocodes = @from row in base_data begin
    @group row by row.ID into grp
    @select { ID=key(grp), codes=Set(grp.NanoCode) }
    @collect DataFrame
end

costs = @from row in base_data begin
    @group row by row.ID into grp
    @select { ID=key(grp), Days=maximum(grp.Days), Cost=maximum(grp.Cost) }
    @collect DataFrame
end

hfacs_data = base_data[:, Not([:NanoCode, :NanoCode_Reason, :Days, :Cost])]
unique!(hfacs_data)
hfacs_data = innerjoin(hfacs_data, costs, on=:ID)

nanocode_index = sort(collect(union(nanocodes.codes...)))

nanocode_names = Dict()
for row in eachrow(unique(raw_data[:, [:NanoCode, :NanoCode_Name]]))
    nanocode_names[row.NanoCode] = row.NanoCode_Name
end

codes = Array{Bool, 2}(undef, nrow(nanocodes), length(nanocode_index))
for (i, code) in enumerate(nanocode_index)
    codes[:, i] = code .∈ nanocodes.codes
end
hfacs_data = hcat(hfacs_data, DataFrame(codes, collect(nanocode_index)))

frequency=Dict()
for code in nanocode_index
    frequency[code, "All"] = sum(hfacs_data[:, code])
    frequency[code, "Armstrong"] = sum(hfacs_data[hfacs_data.Center .== "Armstrong", code])
end

function get_desc(ID::String)
    desc = hfacs_data[hfacs_data.ID .== ID, :Desc]
    println(desc)
    return desc
end

function get_hfacs_desc(ID::String, code::Symbol)
    desc = base_data[(base_data.ID .== ID) .& (base_data.NanoCode .== code), :NanoCode_Reason]
    println(desc)
    return desc
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
function find_index_pairs(; α=0.05, center="All", min_n=2)
    results = DataFrame(Index1=Symbol[], Index2=Symbol[], V=Float64[], pvalue=Float64[], N1=Int[], N2=Int[])
    for i in nanocode_index
        for j in nanocode_index
            if i ≥ j continue end
            V, X2, pvalue = categorical_correlation(i, j)
            if pvalue < α && (frequency[i, center] ≥ min_n && frequency[j, center] ≥ min_n)
                push!(results, (i, j, V, pvalue, frequency[i, center], frequency[j, center]))
            end
        end
    end
    return sort(results, [:pvalue])
end

# find nanocodes correlated with specific index
function find_index_match(index=:Type; α=0.05, center="All", min_n=2)
    results = DataFrame(Index1=Symbol[], Index2=Symbol[], V=Float64[], pvalue=Float64[], N=Int[])
    for j in nanocode_index
        if j == index continue end
        V, X2, pvalue = categorical_correlation(index, j)
        if pvalue < α && frequency[j, center] ≥ min_n
            push!(results, (index, j, V, pvalue, frequency[j, center]))
        end
    end
    return sort(results, [:pvalue])
end

# analyses

# mishap totals
println("Mishap totals")
for level in levels(hfacs_data.Type)
    println("$level:\n \t\tAll Centers \t$(sum(hfacs_data.Type .== level)), \tArmstrong \t$(sum(hfacs_data.Type[hfacs_data.Center .== "Armstrong"] .== level))")
end

#find dirty dozen nanocodes
results = DataFrame(Code=Symbol[], Center=String[], Count=Int[])
for index in nanocode_index
    push!(results, (index, "All", frequency[index, "All"]))
    push!(results, (index, "Armstrong", frequency[index, "Armstrong"]))
end
sort!(results, [:Center, :Count], rev=true)
nasa_dirty_d = Set(first(results.Code[results.Center .== "All", :],12))
armstrong_dirty_d = Set(first(results.Code[results.Center .== "Armstrong", :],12))

println("\nDirty dozen")
println("All Centers")
for (i, row) in enumerate(eachrow(first(results[results.Center .== "All", :],12)))
    println("$i. $(nanocode_names[row.Code]) ($(row.Code)) (n=$(row.Count))")
end

println("Armstrong")
for (i, row) in enumerate(eachrow(first(results[results.Center .== "Armstrong", :],12)))
    println("$i. $(nanocode_names[row.Code]) ($(row.Code)) (n=$(row.Count))")
end

#find damage associations
println("\nDamage association")
println("All Centers")
results = find_index_match(:Damage; α=0.05, center="All", min_n=4)
for (i, row) in enumerate(eachrow(results))
    println("$i. $(nanocode_names[row.Index2]) (p=$(round(row.pvalue, sigdigits=3)), n=$(row.N))")
end

println("\nArmstrong")
results = find_index_match(:Damage; α=0.05, center="Armstrong", min_n=4)
for (i, row) in enumerate(eachrow(results))
    println("$i. $(nanocode_names[row.Index2]) (p=$(round(row.pvalue, sigdigits=3)), n=$(row.N))")
end

#find injury associations
println("\nInjury association")
println("All Centers")
results = find_index_match(:Injury; α=0.05, center="All", min_n=4)
for (i, row) in enumerate(eachrow(results))
    println("$i. $(nanocode_names[row.Index2]) (p=$(round(row.pvalue, sigdigits=3)), n=$(row.N))")
end

println("\nArmstrong")
results = find_index_match(:Injury; α=0.05, center="Armstrong", min_n=4)
for (i, row) in enumerate(eachrow(results))
    println("$i. $(nanocode_names[row.Index2]) (p=$(round(row.pvalue, sigdigits=3)), n=$(row.N))")
end

#find pairs
println("\nAssociated pairs")
println("All Centers")
pairs = first(find_index_pairs(; α=0.001, center="All", min_n=4),20)
for (i, row) in enumerate(eachrow(pairs))
    println("$i. $(nanocode_names[row.Index1]) ($(row.Index1)) & $(lowercase(nanocode_names[row.Index2])) ($(row.Index2)) (p=$(round(row.pvalue, sigdigits=3)), n1=$(row.N1), n2=$(row.N2))")
end

println("\nArmstrong")
pairs = first(find_index_pairs(; α=0.001, center="Armstrong", min_n=4),20)
for (i, row) in enumerate(eachrow(pairs))
    println("$i. $(nanocode_names[row.Index1]) ($(row.Index1)) & $(lowercase(nanocode_names[row.Index2])) ($(row.Index2)) (p=$(round(row.pvalue, sigdigits=3)), n1=$(row.N1), n2=$(row.N2))")
end

pc007_av002_pp104_si001 = @from row in hfacs_data begin
    @where row.PC007 == true && row.AV002 == true && row.PP104 == true && row.SI001 == true
    @select { row.ID, row.Days, row.Cost }
    @collect DataFrame
end

trng_pubs_system = @from row in hfacs_data begin
    @where row.OP005 == true && row.OP006 == true && row.PT009
    @select { row.ID, row.Center, row.Days, row.Cost }
    @collect DataFrame
end

most_cost_all = @from row in hfacs_data begin
    @where !isna(row.Cost)
    @join row2 in nanocodes on row.ID equals row2.ID
    @select { row.ID, row.Center, row.Days, row.Cost, Count=length(intersect(row2.codes,nasa_dirty_d)), Match=intersect(row2.codes,nasa_dirty_d) }
    @collect DataFrame
end
first(sort!(most_cost_all, [:Cost], rev=true),6)

most_cost = @from row in hfacs_data begin
    @where row.Center == "Armstrong" && !isna(row.Cost)
    @join row2 in nanocodes on row.ID equals row2.ID
    @select { row.ID, row.Center, row.Days, row.Cost, Count=length(intersect(row2.codes,nasa_dirty_d)), Match=intersect(row2.codes,nasa_dirty_d) }
    @collect DataFrame
end
first(sort!(most_cost, [:Cost], rev=true),6)

most_days_all = @from row in hfacs_data begin
    @where !isna(row.Days)
    @join row2 in nanocodes on row.ID equals row2.ID
    @select { row.ID, row.Center, row.Days, row.Cost, Count=length(intersect(row2.codes,nasa_dirty_d)), Match=intersect(row2.codes,nasa_dirty_d) }
    @collect DataFrame
end
first(sort!(most_days_all, [:Days], rev=true),6)

most_days = @from row in hfacs_data begin
    @where row.Center == "Armstrong" && !isna(row.Days)
    @join row2 in nanocodes on row.ID equals row2.ID
    @select { row.ID, row.Center, row.Days, row.Cost, Count=length(intersect(row2.codes,nasa_dirty_d)), Match=intersect(row2.codes,nasa_dirty_d) }
    @collect DataFrame
end
first(sort!(most_days, [:Days], rev=true),6)
