function printinput(chromosome, referencepath, targetpath, ancestrypath, threshold, printlevel, maf, totcm, maxblockfraction)

    println("")
    println("LocalAncestry.jl v$(pkgversion(LocalAncestry))")
    println("Started at $(Dates.format(now(), "HH:MM:SS on u d, Y")  )")
    println("")
    println("Files")
    println("   Reference        $(referencepath)")
    println("   Target           $(targetpath)")
    println("   Ancestry         $(ancestrypath)")
    println("")
    println("Parameters")
    println("   Chromosome       $(chromosome)")
    println("   MAF threshold    $(maf)")
    println("   IA threshold     $(threshold)")
    println("   Genetic length   $(totcm) cM")
    println("   Max block frac.  $(maxblockfraction)")
    println("   N threads        $(nthreads())")

    if printlevel == "debug"
        println("")
        println("Debug")
        println("   Hostname         $(gethostname())")
        println("   CPU              $(_print_cpu())")
   end

end

function _print_cpu()
    try
        if Sys.islinux()
            return readchomp(pipeline(`cat /proc/cpuinfo`, `grep "model name"`))
        elseif Sys.isapple()
            return readchomp(`sysctl -n machdep.cpu.brand_string`)
        elseif Sys.iswindows()
            return replace(split(readchomp(pipeline(`wmic cpu get name`)), '\n')[2], '\r' => "")
        else
            return "Unsupported OS for CPU name detection."
        end
    catch
        return "CPU name not available"
    end
end
