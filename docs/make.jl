using FourierSplines
using Documenter

DocMeta.setdocmeta!(FourierSplines, :DocTestSetup, :(using FourierSplines); recursive=true)

makedocs(;
    modules=[FourierSplines],
    authors="Adam Zaidan <azaidan2020@fau.edu>",
    repo="https://github.com/NuclearPy/FourierSplines.jl/blob/{commit}{path}#{line}",
    sitename="FourierSplines.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://NuclearPy.github.io/FourierSplines.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/NuclearPy/FourierSplines.jl",
    devbranch="master",
)
