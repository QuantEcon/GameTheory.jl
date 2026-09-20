# GAMPayoffVector #

"""
    GAMPayoffVector{N,T}

Intermediate representation that stores payoffs in a single flat vector.

Payoff values are ordered as in the GameTracer .gam format:
1. Player-major blocks: player 1, ..., player N.
2. Within each block, action profiles are ordered with player 1 varying fastest,
   then player 2, ..., player N (i.e., Fortran/column-major order).

# Fields

- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `payoffs::Vector{T}` : Vector storing payoffs in .gam order.
"""
struct GAMPayoffVector{N,T<:Real}
    nums_actions::NTuple{N,Int}
    payoffs::Vector{T}

    function GAMPayoffVector{N,T}(
        nums_actions::NTuple{N,Int}, payoffs::Vector{T}
    ) where {N,T<:Real}
        any(n -> n <= 0, nums_actions) &&
            throw(ArgumentError("all nums_actions must be positive"))
        expected = prod(nums_actions) * N
        length(payoffs) == expected || throw(ArgumentError(
            "payoffs length mismatch: expected $expected, got $(length(payoffs))"
        ))
        return new(nums_actions, payoffs)
    end
end

num_players(::GAMPayoffVector{N}) where {N} = N

GAMPayoffVector(
    nums_actions::NTuple{N,Int}, payoffs::Vector{T}
) where {N,T<:Real} = GAMPayoffVector{N,T}(nums_actions, payoffs)

GAMPayoffVector(
    ::Type{T}, nums_actions::NTuple{N,Int}, payoffs::AbstractVector
) where {N,T<:Real} =
    GAMPayoffVector{N,T}(nums_actions, convert(Vector{T}, payoffs))
GAMPayoffVector(
    nums_actions::NTuple{N,Int}, payoffs::AbstractVector{T}
) where {N,T<:Real} = GAMPayoffVector(T, nums_actions, payoffs)


# Forward: (i, i+1, ..., N, 1, ..., i-1)
@inline _perm_fwd(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(i + k - 1, N), Val(N))

# Backward: inverse rotation
@inline _perm_back(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(k - i + 1, N), Val(N))


"""
    GAMPayoffVector([T], g)

Construct a GAMPayoffVector (of eltype `T` if specified) from a NormalFormGame
`g`.

# Examples

```julia
julia> player1 = Player([1 4; 2 5; 3 6]);

julia> player2 = Player([7 8 9; 10 11 12]);

julia> g = NormalFormGame(player1, player2)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)

julia> p = GameTheory.GAMPayoffVector(g);

julia> @show p.payoffs;
p.payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
```
"""
function GAMPayoffVector(::Type{T}, g::NormalFormGame{N}) where {N,T<:Real}
    nums_actions = g.nums_actions
    na = prod(nums_actions)
    payoffs = Vector{T}(undef, na*N)

    ntuple(Val(N)) do i
        copyto!(
            reshape(view(payoffs, na*(i-1)+1:na*i), nums_actions),
            PermutedDimsArray(g.players[i].payoff_array, _perm_back(Val(N), Val(i)))
        )
        nothing
    end

    return GAMPayoffVector{N,T}(nums_actions, payoffs)
end

GAMPayoffVector(g::NormalFormGame{N,T}) where {N,T<:Real} = GAMPayoffVector(T, g)


"""
    NormalFormGame([T], p)

Construct a NormalFormGame (of eltype `T` if specified) from a GAMPayoffVector
`p`.

# Examples

```julia
julia> nums_actions = (3, 2);

julia> payoffs = collect(1:12);

julia> @show payoffs;
payoffs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

julia> p = GameTheory.GAMPayoffVector(nums_actions, payoffs);

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)
```
"""
function NormalFormGame(::Type{T}, p::GAMPayoffVector{N}) where {N,T<:Real}
    nums_actions = p.nums_actions
    na = prod(nums_actions)

    players = ntuple(Val(N)) do i
        Player(
            T,
            PermutedDimsArray(
                reshape(view(p.payoffs, na*(i-1)+1:na*i), nums_actions),
                _perm_fwd(Val(N), Val(i))
            )
        )
    end

    return NormalFormGame{N,T}(players, nums_actions)
end

NormalFormGame(p::GAMPayoffVector{N,T}) where {N,T<:Real} = NormalFormGame(T, p)


# .gam reader and writer #

# The GameTracer .gam format is a whitespace-separated text format: the number
# of players N, the N numbers of actions, and then the prod(nums_actions) * N
# payoffs in the order described in the docstring of `GAMPayoffVector`.
# Reference: B. Blum, D. Koller, and C. Shelton, "Game Theory: GameTracer",
# http://dags.stanford.edu/Games/gametracer.html

"""
    read_gam([T], io)
    read_gam([T], path)

Read a normal form game in the GameTracer .gam format from the stream `io` or
the file at `path`, and return it as a `NormalFormGame`. See
[`GAMPayoffVector`](@ref) for the ordering of the payoffs in the format, and
[`parse_gam`](@ref) for reading from a string.

# Arguments

- `T::Type` : Element type of the payoffs, where `T<:Real`. If omitted, `Int`
  when every payoff in the input is an integer and `Float64` otherwise.
- `io::IO` : Input stream.
- `path::AbstractString` : Path to the file to read.

# Returns

- `::NormalFormGame{N,T}` : The game described by the input.

# Examples

```julia
julia> g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]));

julia> path = tempname();

julia> write_gam(path, g)

julia> read_gam(path)
3×2 NormalFormGame{2, Int64}:
 (3, 2)  (1, 3)
 (0, 6)  (4, 0)
 (2, 1)  (5, 4)

julia> read_gam(Float64, path)
3×2 NormalFormGame{2, Float64}:
 (3.0, 2.0)  (1.0, 3.0)
 (0.0, 6.0)  (4.0, 0.0)
 (2.0, 1.0)  (5.0, 4.0)
```

A file at a URL can be read with `read_gam(Downloads.download(url))`.
"""
read_gam(io::IO) = _read_gam(_parse_payoffs, io)
read_gam(::Type{T}, io::IO) where {T<:Real} =
    _read_gam(tokens -> _parse_payoffs(T, tokens), io)

read_gam(path::AbstractString) = open(read_gam, path)
read_gam(::Type{T}, path::AbstractString) where {T<:Real} =
    open(io -> read_gam(T, io), path)

"""
    parse_gam([T], text)

Parse the string `text` in the GameTracer .gam format and return the game as a
`NormalFormGame`. See [`read_gam`](@ref) for the meaning of `T` and for
reading from a stream or a file.

# Arguments

- `T::Type` : Element type of the payoffs, where `T<:Real`. If omitted, `Int`
  when every payoff in `text` is an integer and `Float64` otherwise.
- `text::AbstractString` : String in the .gam format.

# Returns

- `::NormalFormGame{N,T}` : The game described by `text`.

# Examples

```julia
julia> s = \"\"\"
       2
       3 2

       3 2 0 3 5 6 3 2 3 2 6 1
       \"\"\";

julia> parse_gam(s)
3×2 NormalFormGame{2, Int64}:
 (3, 3)  (3, 2)
 (2, 2)  (5, 6)
 (0, 3)  (6, 1)

julia> parse_gam(Float64, s)
3×2 NormalFormGame{2, Float64}:
 (3.0, 3.0)  (3.0, 2.0)
 (2.0, 2.0)  (5.0, 6.0)
 (0.0, 3.0)  (6.0, 1.0)
```
"""
parse_gam(text::AbstractString) = read_gam(IOBuffer(text))
parse_gam(::Type{T}, text::AbstractString) where {T<:Real} =
    read_gam(T, IOBuffer(text))

function _read_gam(parse_payoffs, io::IO)
    tokens = split(read(io, String))
    isempty(tokens) && throw(ArgumentError("empty .gam input"))

    # Header: N, then the N numbers of actions
    N = parse(Int, tokens[1])
    N > 0 || throw(ArgumentError("number of players must be positive"))
    length(tokens) >= N + 1 || throw(ArgumentError(
        "incomplete header: expected $N numbers of actions, got $(length(tokens)-1)"
    ))
    nums_actions = ntuple(i -> parse(Int, tokens[i+1]), N)

    # Payoffs, in .gam order; the length is checked by GAMPayoffVector
    payoffs = parse_payoffs(@view tokens[N+2:end])

    return NormalFormGame(GAMPayoffVector(nums_actions, payoffs))
end

# Number parsing, shared by the readers #

# A number token is an integer, a decimal with an optional exponent, or a
# rational `n/d`, each with an optional sign. When no element type is given,
# it is Int if every token is an integer, Rational{BigInt} if the other tokens
# are all rationals, and Float64 otherwise.
function _parse_payoffs(tokens)
    payoffs = Vector{Int}(undef, length(tokens))
    for (i, tok) in enumerate(tokens)
        x = tryparse(Int, tok)
        if x === nothing
            isexact(t) = occursin('/', t) || tryparse(Int, t) !== nothing
            T = all(isexact, @view tokens[i:end]) ? Rational{BigInt} : Float64
            return _parse_payoffs(T, tokens)
        end
        payoffs[i] = x
    end
    return payoffs
end

_parse_payoffs(::Type{T}, tokens) where {T<:Real} =
    T[_parse_payoff(T, tok) for tok in tokens]

_parse_payoff(::Type{T}, tok) where {T<:Integer} = parse(T, tok)
_parse_payoff(::Type{T}, tok) where {T<:AbstractFloat} =
    occursin('/', tok) ? T(_parse_exact(tok)) : parse(T, tok)
# Other types, such as `Rational`, get the exact value of the token
_parse_payoff(::Type{T}, tok) where {T<:Real} = T(_parse_exact(tok))

# Exact value of a number token: "3", "-1/3", "0.1", "-12.5e-3"
function _parse_exact(tok::AbstractString)
    if occursin('/', tok)
        n, d = split(tok, '/'; limit=2)
        d = parse(BigInt, d; base=10)
        iszero(d) && throw(ArgumentError("zero denominator in $(repr(tok))"))
        return parse(BigInt, n; base=10) // d
    end
    mantissa, ex = occursin(r"[eE]", tok) ? split(tok, r"[eE]"; limit=2) : (tok, "0")
    int, frac = occursin('.', mantissa) ? split(mantissa, '.'; limit=2) : (mantissa, "")
    num = parse(BigInt, int * frac; base=10)
    e = parse(Int, ex) - length(frac)
    return e >= 0 ? num * big(10)^e // 1 : num // big(10)^(-e)
end


# Number printing, shared by the writers #

# Element types that the writers accept
const _PayoffNumber = Union{Integer,AbstractFloat,Rational}

_print_payoff(io::IO, x::Real) = print(io, x)
# A rational is written as `n/d`, or as `n` if the denominator is 1
function _print_payoff(io::IO, x::Rational)
    print(io, numerator(x))
    isone(denominator(x)) || print(io, '/', denominator(x))
    return nothing
end


"""
    write_gam(io, g)
    write_gam(path, g)

Write the game `g` to the stream `io` or the file at `path` in the GameTracer
.gam format. Each payoff is written with `print`, so the element type of `g`
must be an `Integer` or an `AbstractFloat` type; convert first otherwise, e.g.
with `NormalFormGame(Float64, g)`. See [`gam_string`](@ref) for writing to a
string.

# Arguments

- `io::IO` : Output stream.
- `path::AbstractString` : Path to the file to write; an existing file is
  overwritten.
- `g::Union{NormalFormGame,GAMPayoffVector}` : Game to write.

# Examples

```julia
julia> g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]));

julia> write_gam(stdout, g)
2
3 2

3 0 2 1 4 5 2 6 1 3 0 4

julia> write_gam("game.gam", g)
```
"""
function write_gam(io::IO, p::GAMPayoffVector{N,T}) where {N,T<:_PayoffNumber}
    print(io, N, '\n')
    join(io, p.nums_actions, ' ')
    print(io, "\n\n")  # blank line between the header and the payoffs
    for (k, x) in enumerate(p.payoffs)
        k > 1 && print(io, ' ')
        _print_payoff(io, x)
    end
    print(io, '\n')
    return nothing
end

write_gam(io::IO, g::NormalFormGame{N,T}) where {N,T<:_PayoffNumber} =
    write_gam(io, GAMPayoffVector(g))

# Same bound on the element type as the methods for `io`, so that an
# unsupported game is rejected before the file is opened
write_gam(
    path::AbstractString, g::Union{NormalFormGame{N,T},GAMPayoffVector{N,T}}
) where {N,T<:_PayoffNumber} = open(io -> write_gam(io, g), path, "w")

"""
    gam_string(g)

Return the GameTracer .gam representation of the game `g` as a string. See
[`write_gam`](@ref) for the requirement on the element type of `g`.

# Arguments

- `g::Union{NormalFormGame,GAMPayoffVector}` : Game to write.

# Returns

- `::String` : The .gam representation of `g`.

# Examples

```julia
julia> g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]));

julia> gam_string(g)
"2\\n3 2\\n\\n3 0 2 1 4 5 2 6 1 3 0 4\\n"

julia> print(gam_string(g))
2
3 2

3 0 2 1 4 5 2 6 1 3 0 4
```
"""
gam_string(g::Union{NormalFormGame,GAMPayoffVector}) = sprint(write_gam, g)
