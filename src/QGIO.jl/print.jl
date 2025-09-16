function _print_ancestries(x) 
    tmp = sort(combine(
        groupby(x, "population"),
        [
            "population" => (x -> sum(x .!= "unknown")) => "col1",
            "invcf" => (x -> sum(x .> 0)) => "col2",
            "omit" => (x -> sum(x)) => "col3",
            ["population", "invcf", "omit"] => ((x, y, z) -> sum((y .> 0) .& (.~z))) => "used"
        ]
    ))
  
    println("")
    println("Number of haplotypes per population")
    append!(tmp, DataFrame(population = "Any", col1 = sum(tmp.col1), col2 = sum(tmp.col2), col3 = sum(tmp.col3), used = sum(tmp.used)))
    pretty_table(tmp, header=["Reference population", "Ancestries", "Haplotypes", "Omitted", "Included"], hlines =[0,1,nrow(tmp),nrow(tmp)+1])
end
