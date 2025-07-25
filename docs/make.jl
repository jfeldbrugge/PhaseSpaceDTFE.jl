ENV["JULIA_PROGRESS"] = "false"

using PhaseSpaceDTFE
using Documenter

DocMeta.setdocmeta!(PhaseSpaceDTFE, :DocTestSetup, :(using PhaseSpaceDTFE); recursive=true)

makedocs(;
    modules=[PhaseSpaceDTFE],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="PhaseSpaceDTFE.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/PhaseSpaceDTFE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/PhaseSpaceDTFE.jl",
    push_preview = true,  # important for GitHub Actions!
    devbranch="main",
)
