using RNAInferenceTool
using Documenter

DocMeta.setdocmeta!(RNAInferenceTool, :DocTestSetup, :(using RNAInferenceTool); recursive=true)

makedocs(;
    modules=[RNAInferenceTool],
    authors="Xiaoming Fu",
    repo="https://github.com/palmtree2013/RNAInferenceTool.jl/blob/{commit}{path}#{line}",
    sitename="RNAInferenceTool.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://palmtree2013.github.io/RNAInferenceTool.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/palmtree2013/RNAInferenceTool.jl",
    devbranch="main",
)
