"""
Write GameTracer `.gam` text for `x` to `io`.

Step IR boundary:
- `p = GamPayoffVector(x)` must yield:
  - `p.nums_actions::NTuple{N,Int}`
  - `p.payoffs::Vector{<:Real}` already in `.gam` flat order.

This writer does not reorder; it only serializes.
No trailing newline is written.
"""
function write_gam(io::IO, x)
    p = GamPayoffVector(x)  # supports NormalFormGame and GamPayoffVector

    N = length(p.nums_actions)

    # Header
    print(io, N, '\n')
    join(io, p.nums_actions, ' ')
    print(io, "\n\n")

    # Payoffs (already `.gam`-ordered)
    join(io, p.payoffs, ' ')

    return nothing
end


"""
    read_gam(io::IO) -> NormalFormGame

Read a GameTracer `.gam` payload from an IO stream, funneling through the IR:

    tokens -> (N, nums_actions, payoffs) -> GAMPayoffVector -> NormalFormGame

This uses the same numeric token rule as Python `_str2num`:

- if token contains `.` => parse as Float64
- else => parse as Int

and mimics NumPy promotion: if *any* payoff token contains `.`, then the
entire payoff vector is Float64 (with integer tokens parsed as Int then
converted).
"""
function read_gam(io::IO)
    tokens = split(read(io, String))
    isempty(tokens) && throw(ArgumentError("empty .gam input"))

    # Parse N at the dynamic boundary (data-dependent)
    tokN = tokens[1]
    N = try
        parse(Int, tokN)
    catch
        throw(ArgumentError("invalid N token: $(repr(tokN))"))
    end
    N > 0 || throw(ArgumentError("N must be a positive integer"))

    # Minimal split (function barrier): specialize downstream on N via Val(N)
    return _read_gam(Val(N), tokens)
end

# Minimal helper: once we're here, N is a type parameter, so nums_actions can be NTuple{N,Int}
function _read_gam(::Val{N}, tokens::AbstractVector{<:AbstractString}) where {N}
    pos = 2  # tokens[1] was N

    # Header: nums_actions
    if length(tokens) < pos + N - 1
        got = max(0, length(tokens) - pos + 1)
        throw(ArgumentError("incomplete header: expected $N action counts, got $got"))
    end

    nums_actions = ntuple(Val(N)) do i
        tok = tokens[pos + i - 1]
        try
            parse(Int, tok)
        catch
            throw(ArgumentError("invalid action count token in header: $(repr(tok))"))
        end
    end
    pos += N

    # Payoffs: player-major blocks, profile order as in .gam (we preserve token order)
    payoff_tokens = @view tokens[pos:end]

    # Promotion rule consistent with Python np.array([...]) over _str2num:
    # if any token contains '.', allocate Float64; else allocate Int.
    has_float = any(tok -> occursin(".", tok), payoff_tokens)

    if has_float
        payoffs = Vector{Float64}(undef, length(payoff_tokens))
        @inbounds for i in eachindex(payoff_tokens)
            tok = payoff_tokens[i]
            if occursin(".", tok)
                payoffs[i] = parse(Float64, tok)
            else
                payoffs[i] = Float64(parse(Int, tok))  # preserves Python's rejection of "1e3"
            end
        end
        pv = GAMPayoffVector(nums_actions, payoffs)   # => GAMPayoffVector{N,Float64}
        return NormalFormGame(pv)                     # => NormalFormGame{N,Float64}
    else
        payoffs = Vector{Int}(undef, length(payoff_tokens))
        @inbounds for i in eachindex(payoff_tokens)
            payoffs[i] = parse(Int, payoff_tokens[i])
        end
        pv = GAMPayoffVector(nums_actions, payoffs)   # => GAMPayoffVector{N,Int}
        return NormalFormGame(pv)                     # => NormalFormGame{N,Int}
    end
end

# --- Thin wrappers (Step API) ---

read_gam(path::AbstractString) = open(read_gam, path)

parse_gam(text::AbstractString) = read_gam(IOBuffer(text))

function read_gam_url(url::AbstractString)
    import Downloads
    buf = IOBuffer()
    Downloads.request(url; output=buf)
    seekstart(buf)
    return read_gam(buf)
end
