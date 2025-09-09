"""
    $TYPEDEF

Acts as an interface between a given `InputFields` container and one or more `AbstractInputSource`s
that are responsible for initializing and updating the input `Field`s at each time step.
"""
struct InputProvider{
    NF,
    Grid<:AbstractLandGrid{NF},
    Inputs<:InputFields{NF, Grid},
    SourceType<:AbstractInputSource
}
    "Input variable `Field` states"
    fields::Inputs

    "A list of input sources to be invoked by `update_inputs!`"
    sources::Vector{SourceType}
end

InputProvider(grid::AbstractLandGrid, sources::AbstractInputSource...) = InputProvider(InputFields(grid), collect(sources))

function update_inputs!(provider::InputProvider, clock::Clock)
    for source in provider.sources
        update_inputs!(provider.fields, source, clock)
    end
end
