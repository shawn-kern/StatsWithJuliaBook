using DataFrames, CSV, Query

raw_data = CSV.read("hfacs/HFACS_FY19-22_5YR.csv", DataFrame, types=Dict(:Type=>String, :NanoCode=>Symbol), copycols=true)

function convert_mishap_type(value)
    if value == "Close Call" return "L"
    elseif value == "Non-8621 Reportable" return "O"
    else return value[end-1:end]
    end
end

base_data = @from row in raw_data begin
  @where !isna(row.NanoCode)
  @select { row.ID, Center=split(row.Center)[1], Type=convert_mishap_type(row.Type), 
    Injury=row.Injury=="Yes", Damage=row.Damage=="Yes", row.NanoCode, Days=row.OSHA_Days_Away,
   Cost=row.Final_Cost }
  @collect DataFrame
end

nanocodes = @from row in base_data begin
    @group row by row.ID into grp
    @select { ID=key(grp), codes=collect(grp.NanoCode) }
    @collect DataFrame
end

newcol = Set{Symbol}[]
for row in eachrow(nanocodes)
    S = Set{Symbol}()
    for item in row.codes
        push!(S, get(item))
    end
    push!(newcol, S)
end
nanocodes.codes = newcol

