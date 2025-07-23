using PSDTFE
using Documenter

DocMeta.setdocmeta!(PSDTFE, :DocTestSetup, :(using PSDTFE); recursive=true)

makedocs(;
    modules=[PSDTFE],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="PSDTFE.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/PSDTFE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/PSDTFE.jl",
    push_preview = true,  # important for GitHub Actions!
    devbranch="main",
)
