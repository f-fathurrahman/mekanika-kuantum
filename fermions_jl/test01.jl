using Infiltrator
using LinearAlgebra: eigen, Hermitian

function BasisStates(
        numLevels::Int64, 
        totOccReq::Vector{Int64},
        magzReq::Vector{Int64},
        localCriteria::Function
    )
    @assert !isempty(totOccReq) && !isempty(magzReq)
    basis = Dict{BitVector,Float64}[]
    for decimalNum in 0:2^numLevels-1
        config = digits(decimalNum, base=2, pad=numLevels) |> reverse
        if !isempty(totOccReq)
            totOcc = sum(config)
            if totOcc ∉ totOccReq
                continue
            end
        end
        if !isempty(magzReq)
            magz = sum(config[1:2:end]) - sum(config[2:2:end])
            if magz ∉ magzReq
                continue
            end
        end
        if localCriteria(config)
            push!(basis, Dict(BitVector(config) => 1.0))
        end
    end
    return basis
end


function BasisStates(
        numLevels::Int64;
        totOccReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        magzReq::Union{Vector{Int64}, Int64, Nothing}=nothing,
        localCriteria::Function=x -> true
    )
    if isnothing(totOccReq)
        totOccReq = collect(0:numLevels)
    elseif typeof(totOccReq) == Int64
        totOccReq = [totOccReq]
    end
    if isnothing(magzReq)
        magzReq = collect(-div(numLevels, 2):numLevels - div(numLevels, 2))
    elseif typeof(magzReq) == Int64
        magzReq = [magzReq]
    end
    return BasisStates(numLevels, totOccReq, magzReq, localCriteria)
end

# designing the tight-binding Hamiltonian
function TightBindHamiltonian(numSites)

    # define the array of tuples to store all terms in the Hamiltonian 
    hamiltonianTerms = Tuple{String, Vector{Int64}, Float64}[]

    # loop over all lattice sites (skip the end, because we have c^\dagger_i+1 
    for i in 1:numSites-1

        # the term c^\dagger_i c_i+1
        term1 = ("+-", [i, i+1], 1.0)

        # the term c^\dagger_i+1 c_i
        term2 = ("+-", [i+1, i], 1.0)

        # add both terms to the Hamiltonian
        push!(hamiltonianTerms, term1)
        push!(hamiltonianTerms, term2)
    end
    return hamiltonianTerms
end

function TransformBit(qubit::Bool, operator::Char)
    @assert operator in ('n', 'h', '+', '-')
    if operator == 'n'
        return 0 + qubit, 0 + qubit
    elseif operator == 'h'
        return 0 + qubit, 1 - qubit
    elseif (operator == '+' && qubit == 0) || (operator == '-' && qubit == 1)
        return 1 - qubit, 1
    else
        return 0 + qubit, 0
    end
end

function StateOverlap(
        state1::Dict{BitVector,Float64}, 
        state2::Dict{BitVector,Float64}
    )
    overlap = 0.
    keys2 = keys(state2)
    for (key, val) in state1
        if key ∈ keys2
            overlap += val * state2[key]
        end
    end
    return overlap
end


function Spectrum(
    operator::Vector{Tuple{String,Vector{Int64},Float64}},
    basisStates::Vector{Dict{BitVector,Float64}};
    diagElements::Vector{Float64}=Float64[],
    tolerance::Float64=1e-14,
    assumeHerm::Bool=true,
    )
    matrix = OperatorMatrix(basisStates, operator; tolerance=tolerance)
    if assumeHerm
        @assert maximum(abs.(matrix .- matrix')) < tolerance
    end
    if !isempty(diagElements)
        matrix += diagm(diagElements)
    end

    eigenVals, eigenVecs = assumeHerm ? eigen(Hermitian(matrix)) : eigen(matrix)
    eigenStates = [TransformState(collect(vector), basisStates; tolerance=tolerance)
                   for vector in eachcol(eigenVecs)]
    return eigenVals, eigenStates
end


function OperatorMatrix(
        basisStates::Vector{Dict{BitVector,Float64}},
        operator::Vector{Tuple{String,Vector{Int64},Float64}};
        tolerance::Float64=1e-16,
    )
    operatorMatrix = zeros(length(basisStates), length(basisStates))
    newStates = [ApplyOperator(operator, incomingState; tolerance=tolerance)
                        for incomingState in basisStates]
    for incomingIndex in findall(!isempty, newStates)
        for outgoingIndex in eachindex(basisStates)
            operatorMatrix[outgoingIndex, incomingIndex] = StateOverlap(basisStates[outgoingIndex], newStates[incomingIndex])
        end
    end
    return operatorMatrix
end


function TransformState(
        vector::Union{Vector{ComplexF64}, Vector{Float64}},
        basisStates::Vector{Dict{BitVector,Float64}};
        tolerance::Float64=1e-16,
    )
    transformedState = Dict{BitVector, vector |> eltype}()
    keysArr = keys.(basisStates)
    valuesArr = values.(basisStates)
    for (i, c_i) in enumerate(vector)
        if abs(c_i) < tolerance
            continue
        end
        mergewith!(+, transformedState, Dict(keysArr[i] .=> c_i .* valuesArr[i]))
    end
    return transformedState
end


function ApplyOperator(
        operator::Vector{Tuple{String,Vector{Int64},Float64}},
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    @assert !isempty(operator)
    @assert maximum([maximum(positions) for (_, positions, _) in operator]) ≤ length.(keys(incomingState))[1]

    return mergewith(+, [ApplyOperatorChunk(opType, opMembers, opStrength, copy(incomingState); tolerance=tolerance) 
                                for (opType, opMembers, opStrength) in operator]...)

    return outgoingState
end


function ApplyOperatorChunk(
        opType::String,
        opMembers::Vector{Int64},
        opStrength::Float64,
        incomingState::Dict{BitVector,Float64};
        tolerance::Float64=1e-16
    )
    outgoingState = Dict{BitVector,Float64}()
    for i in eachindex(opMembers)[end:-1:2]
        if opMembers[i] ∈ opMembers[i+1:end]
            continue
        end
        if opType[i] == '+' || opType[i] == 'h'
            filter!(p -> p[1][opMembers[i]] == 0, incomingState)
        else
            filter!(p -> p[1][opMembers[i]] == 1, incomingState)
        end
    end
    for (incomingBasisState, coefficient) in incomingState

        newCoefficient = coefficient
        outgoingBasisState = copy(incomingBasisState)

        # for each basis state, obtain a modified state after applying the operator tuple
        for (siteIndex, operator) in zip(reverse(opMembers), reverse(opType))
            newQubit, factor = TransformBit(outgoingBasisState[siteIndex], operator)
            if factor == 0
                newCoefficient = 0
                break
            end
            # calculate the fermionic exchange sign by counting the number of
            # occupied states the operator has to "hop" over
            exchangeSign = ifelse(operator in ['+', '-'], (-1)^sum(outgoingBasisState[1:siteIndex-1]), 1)

            outgoingBasisState[siteIndex] = newQubit
            newCoefficient *= exchangeSign * factor
        end

        if abs(newCoefficient) > tolerance
            if haskey(outgoingState, outgoingBasisState)
                outgoingState[outgoingBasisState] += opStrength * newCoefficient
            else
                outgoingState[outgoingBasisState] = opStrength * newCoefficient
            end
        end
    end
    return outgoingState
end

function test_main()
    #=
    numSites = 2 # number of lattice sites
    occupancy = 1
    basis = BasisStates(numSites, totOccReq=occupancy)
    =#

    # quite large number
    numSites = 10 #50
    basis = BasisStates(numSites, totOccReq=1)
    

    hamiltonian = TightBindHamiltonian(numSites)
    eigvals, eigvecs = Spectrum(hamiltonian, basis)

    @infiltrate

end