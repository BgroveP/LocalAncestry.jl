
function header(path)
    file = open_vcf(path)
    o = Vector{String}()
    s = readline(file)
    while ~eof(file) & (s[1:2] == "##")
        push!(o, s)
        s = readline(file)
    end
    close(file)
    return o
end