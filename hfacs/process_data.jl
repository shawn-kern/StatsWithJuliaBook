using DataFrames, CSV, Query, CategoricalArrays, DataStructures

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

nanocode_index = SortedSet(union(nanocodes.codes...))
codes = Array{Bool, 2}(undef, nrow(nanocodes), length(nanocode_index))
for (i, code) in enumerate(nanocode_index)
    codes[:, i] = code .∈ nanocodes.codes
end
hfacs_data = hcat(hfacs_data, DataFrame(codes, collect(nanocode_index)))

# simple analyses
