using Documenter
using LocalAncestry

makedocs(
    sitename="LocalAncestry.jl",
    format=Documenter.HTML(),
    modules=[LocalAncestry],
    clean=true,
    doctest=true,
    highlightsig=true,
    pages=[
        "Overview" => "index.md",
        "Citation" => "citation.md",
        "Release notes" => "releasenotes.md",
        "Troubleshooting" => "troubleshooting.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "github.com/BgroveP/LocalAncestry.jl"
)
