# PayoffVector #

"""
    PayoffLayout

Abstract supertype of the singleton types that specify the ordering of the
payoff values in a [`PayoffVector`](@ref).
"""
abstract type PayoffLayout end

"""
    PlayerMajor <: PayoffLayout

Player-major ordering, as in the GameTracer .gam format: all the payoffs of
player 1 come first, then those of player 2, and so on. Within each block,
action profiles are ordered with player 1 varying fastest, then player 2, ...,
player N (i.e., column-major order).
"""
struct PlayerMajor <: PayoffLayout end

"""
    ProfileMajor <: PayoffLayout

Profile-major ordering, as in the Gambit .nfg format: the payoffs of players
1, ..., N at the first action profile come first, then those at the second
action profile, and so on. Action profiles are ordered with player 1 varying
fastest, then player 2, ..., player N (i.e., column-major order).
"""
struct ProfileMajor <: PayoffLayout end

"""
    PayoffVector{L,N,T}

Intermediate representation that stores the payoffs of an `N`-player
game in a single flat vector of eltype `T`, ordered according to the layout
`L<:PayoffLayout`. See [`GAMPayoffVector`](@ref) and
[`NFGPayoffVector`](@ref) for the two layouts available.

Viewed as a `prod(nums_actions) × N` matrix whose `[a, i]` entry is the payoff
to player `i` at the `a`-th action profile in column-major order, `payoffs` is
that matrix vectorized in column-major order for `PlayerMajor`, and in
row-major order for `ProfileMajor`.

# Fields

- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `payoffs::Vector{T}` : Vector storing payoffs in the order specified by `L`.
"""
struct PayoffVector{L<:PayoffLayout,N,T<:Real}
    nums_actions::NTuple{N,Int}
    payoffs::Vector{T}

    function PayoffVector{L,N,T}(
        nums_actions::NTuple{N,Int}, payoffs::Vector{T}
    ) where {L<:PayoffLayout,N,T<:Real}
        N > 0 || throw(ArgumentError("nums_actions must be non-empty"))
        any(n -> n <= 0, nums_actions) &&
            throw(ArgumentError("all nums_actions must be positive"))
        expected = prod(nums_actions) * N
        length(payoffs) == expected || throw(ArgumentError(
            "payoffs length mismatch: expected $expected, got $(length(payoffs))"
        ))
        return new(nums_actions, payoffs)
    end
end

# NOTE: The parameter bounds of an alias must be identical to those of
# `PayoffVector` (`T<:Real`), so that e.g. `GAMPayoffVector === PayoffVector{PlayerMajor}`
# holds and the `PayoffVector{L}(...)` constructors below apply to the alias.

"""
    GAMPayoffVector{N,T}

Alias for `PayoffVector{PlayerMajor,N,T}`: payoff values are ordered as in the
GameTracer .gam format:
1. Player-major blocks: player 1, ..., player N.
2. Within each block, action profiles are ordered with player 1 varying fastest,
   then player 2, ..., player N (i.e., column-major order).
"""
const GAMPayoffVector{N,T<:Real} = PayoffVector{PlayerMajor,N,T}

"""
    NFGPayoffVector{N,T}

Alias for `PayoffVector{ProfileMajor,N,T}`: payoff values are ordered as in the
Gambit .nfg format:
1. Profile-major blocks: action profiles are ordered with player 1 varying
   fastest, then player 2, ..., player N (i.e., column-major order).
2. Within each block, the payoffs to player 1, ..., player N.
"""
const NFGPayoffVector{N,T<:Real} = PayoffVector{ProfileMajor,N,T}

num_players(::PayoffVector{L,N}) where {L,N} = N

PayoffVector{L}(
    nums_actions::NTuple{N,Int}, payoffs::Vector{T}
) where {L<:PayoffLayout,N,T<:Real} = PayoffVector{L,N,T}(nums_actions, payoffs)

PayoffVector{L}(
    ::Type{T}, nums_actions::NTuple{N,Int}, payoffs::AbstractVector
) where {L<:PayoffLayout,N,T<:Real} =
    PayoffVector{L,N,T}(nums_actions, convert(Vector{T}, payoffs))
PayoffVector{L}(
    nums_actions::NTuple{N,Int}, payoffs::AbstractVector{T}
) where {L<:PayoffLayout,N,T<:Real} = PayoffVector{L}(T, nums_actions, payoffs)


# Forward: (i, i+1, ..., N, 1, ..., i-1)
@inline _perm_fwd(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(i + k - 1, N), Val(N))

# Backward: inverse rotation
@inline _perm_back(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(k - i + 1, N), Val(N))

@inline _colons(::Val{N}) where {N} = ntuple(_ -> Colon(), Val(N))

"""
    _player_block(p, i)

Return a view of `p.payoffs` holding the payoffs to player `i`, as an `N`-dim
array indexed by the action profile `(a_1, ..., a_N)`. This is the only place
where the layout `L` matters; no copy is made.
"""
@inline _player_block(p::PayoffVector{PlayerMajor,N}, i::Int) where {N} =
    view(reshape(p.payoffs, (p.nums_actions..., N)), _colons(Val(N))..., i)
@inline _player_block(p::PayoffVector{ProfileMajor,N}, i::Int) where {N} =
    view(reshape(p.payoffs, (N, p.nums_actions...)), i, _colons(Val(N))...)


"""
    PayoffVector{L}([T], g)

Construct a `PayoffVector` of layout `L` (and of eltype `T` if specified) from
a NormalFormGame `g`. `GAMPayoffVector([T], g)` and `NFGPayoffVector([T], g)`
are the versions for the two layouts.

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

julia> p = GameTheory.NFGPayoffVector(g);

julia> @show p.payoffs;
p.payoffs = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]
```
"""
function PayoffVector{L}(
    ::Type{T}, g::NormalFormGame{N}
) where {L<:PayoffLayout,N,T<:Real}
    nums_actions = g.nums_actions
    payoffs = Vector{T}(undef, prod(nums_actions)*N)
    p = PayoffVector{L,N,T}(nums_actions, payoffs)

    ntuple(Val(N)) do i
        copyto!(
            _player_block(p, i),
            PermutedDimsArray(g.players[i].payoff_array, _perm_back(Val(N), Val(i)))
        )
        nothing
    end

    return p
end

PayoffVector{L}(g::NormalFormGame{N,T}) where {L<:PayoffLayout,N,T<:Real} =
    PayoffVector{L}(T, g)


"""
    PayoffVector{L}([T], p)

Construct a `PayoffVector` of layout `L` (and of eltype `T` if specified) from
a PayoffVector `p` of possibly another layout. The payoffs are copied.

# Examples

```julia
julia> p = GameTheory.GAMPayoffVector((3, 2), collect(1:12));

julia> p_nfg = GameTheory.NFGPayoffVector(p);

julia> @show p_nfg.payoffs;
p_nfg.payoffs = [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]
```
"""
function PayoffVector{L}(
    ::Type{T}, p::PayoffVector{L1,N}
) where {L<:PayoffLayout,L1<:PayoffLayout,N,T<:Real}
    payoffs = Vector{T}(undef, length(p.payoffs))
    p_new = PayoffVector{L,N,T}(p.nums_actions, payoffs)

    ntuple(Val(N)) do i
        copyto!(_player_block(p_new, i), _player_block(p, i))
        nothing
    end

    return p_new
end

PayoffVector{L}(p::PayoffVector{L1,N,T}) where {L<:PayoffLayout,L1<:PayoffLayout,N,T<:Real} =
    PayoffVector{L}(T, p)
PayoffVector{L,N,T}(
    p::PayoffVector{L1,N}
) where {L<:PayoffLayout,N,T<:Real,L1<:PayoffLayout} = PayoffVector{L}(T, p)

# As for `Player` and `NormalFormGame`: `p` itself if it already has the
# layout and eltype, a copy otherwise
Base.convert(::Type{T}, p::PayoffVector) where {T<:PayoffVector} =
    p isa T ? p : T(p)


"""
    NormalFormGame([T], p)

Construct a NormalFormGame (of eltype `T` if specified) from a PayoffVector
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

julia> p = GameTheory.NFGPayoffVector(nums_actions, payoffs);

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 2)  (7, 8)
 (3, 4)  (9, 10)
 (5, 6)  (11, 12)
```
"""
function NormalFormGame(::Type{T}, p::PayoffVector{L,N}) where {L<:PayoffLayout,N,T<:Real}
    players = ntuple(Val(N)) do i
        Player(
            T,
            PermutedDimsArray(_player_block(p, i), _perm_fwd(Val(N), Val(i)))
        )
    end

    return NormalFormGame{N,T}(players, p.nums_actions)
end

NormalFormGame(p::PayoffVector{L,N,T}) where {L<:PayoffLayout,N,T<:Real} =
    NormalFormGame(T, p)


# .gam reader and writer #

# The GameTracer .gam format is a whitespace-separated text format: the number
# of players N, the N numbers of actions, and then the prod(nums_actions) * N
# payoffs in the order described in the docstring of `PlayerMajor`.
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
  when every payoff in the input is written as an integer (`BigInt` if one
  does not fit in `Int`) and `Float64` otherwise.
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
  when every payoff in `text` is written as an integer (`BigInt` if one does
  not fit in `Int`) and `Float64` otherwise.
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
    N < length(tokens) || throw(ArgumentError(
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
# it is Int if every token is an integer, or BigInt if one of them does not
# fit in Int; Rational{BigInt} if the other tokens are all rationals; and
# Float64 otherwise.
function _parse_payoffs(tokens)
    payoffs = Vector{Int}(undef, length(tokens))
    for (i, tok) in enumerate(tokens)
        x = tryparse(Int, tok)
        if x === nothing
            isint(t) = tryparse(BigInt, t) !== nothing
            isexact(t) = occursin('/', t) || isint(t)
            rest = @view tokens[i:end]
            T = all(isint, rest) ? BigInt :
                all(isexact, rest) ? Rational{BigInt} : Float64
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
    # A sign in `frac` would be moved to a valid position by the concatenation
    all(isdigit, frac) ||
        throw(ArgumentError("cannot parse $(repr(tok)) as a number"))
    num = parse(BigInt, int * frac; base=10)
    e = parse(Int, ex) - length(frac)
    return e >= 0 ? num * big(10)^e // 1 : num // big(10)^(-e)
end


# Number printing, shared by the writers #

# Element types that the writers accept
const _PayoffNumber = Union{Integer,AbstractFloat,Rational}

_print_payoff(io::IO, x::Real) = print(io, x)
# `print` would write `true` or `false`
_print_payoff(io::IO, x::Bool) = print(io, Int(x))
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
- `g::Union{NormalFormGame,PayoffVector}` : Game to write. A `PayoffVector` of
  any layout is accepted; one that is not player-major is converted first.

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
    # `print` would round floats if the caller's context has `:compact => true`
    io = IOContext(io, :compact => false)
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

# Reached only for layouts other than PlayerMajor, which the method above
# handles; the conversion copies the payoffs into player-major order
write_gam(
    io::IO, p::PayoffVector{L,N,T}
) where {L<:PayoffLayout,N,T<:_PayoffNumber} = write_gam(io, GAMPayoffVector(p))

write_gam(io::IO, g::NormalFormGame{N,T}) where {N,T<:_PayoffNumber} =
    write_gam(io, GAMPayoffVector(g))

# Same bound on the element type as the methods for `io`, so that an
# unsupported game is rejected before the file is opened
write_gam(
    path::AbstractString,
    g::Union{NormalFormGame{N,T},PayoffVector{<:PayoffLayout,N,T}}
) where {N,T<:_PayoffNumber} = open(io -> write_gam(io, g), path, "w")

"""
    gam_string(g)

Return the GameTracer .gam representation of the game `g` as a string. See
[`write_gam`](@ref) for the requirement on the element type of `g`.

# Arguments

- `g::Union{NormalFormGame,PayoffVector}` : Game to write.

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
gam_string(g::Union{NormalFormGame,PayoffVector}) = sprint(write_gam, g)


# .nfg reader and writer #

# The Gambit .nfg format is a text format with a prologue (the keyword NFG,
# the version, R or D, the title, the names of the players, the numbers of
# actions or the names of the actions, and an optional comment) and a body in
# one of two versions: the payoffs in the order described in the docstring of
# `ProfileMajor`, or a list of outcomes (a name and N payoffs each) followed by
# the index of the outcome at each action profile, 0 for zero payoffs.
# Reference: The Gambit Project, "Game representation formats",
# https://gambitproject.readthedocs.io/en/latest/formats.html

# A token is a quoted string (with `\"` for a quote inside), a brace, or a run
# of other characters; commas are separators
const _NFG_TOKEN = r"\"(?:[^\"\\]|\\.)*\"|[{}]|[^\s{}\",]+"

# Return the item starting at `tokens[pos]` and the position after it: a
# nested vector for a braced group, the token itself otherwise (the Lisp
# reader)
function _read_tree(tokens, pos)
    if tokens[pos] == "{"
        items = Any[]
        pos += 1
        while tokens[pos] != "}"
            item, pos = _read_tree(tokens, pos)
            push!(items, item)
        end
        return items, pos + 1
    end
    return tokens[pos], pos + 1
end

"""
    read_nfg([T], io)
    read_nfg([T], path)

Read a normal form game in the Gambit .nfg format from the stream `io` or the
file at `path`, and return it as a `NormalFormGame`. Both the payoff version
and the outcome version of the format are read; the title, the names of the
players, of the actions, and of the outcomes, and the comment are ignored. See
[`ProfileMajor`](@ref) for the ordering of the payoffs in the format, and
[`parse_nfg`](@ref) for reading from a string.

# Arguments

- `T::Type` : Element type of the payoffs, where `T<:Real`. If omitted, `Int`
  when every payoff in the input is written as an integer (`BigInt` if one
  does not fit in `Int`), `Rational{BigInt}` when the others are written as
  rationals `n/d`, and `Float64` otherwise.
- `io::IO` : Input stream.
- `path::AbstractString` : Path to the file to read.

# Returns

- `::NormalFormGame{N,T}` : The game described by the input.

# Examples

```julia
julia> g = NormalFormGame(Player([3 1; 0 4; 2 5]), Player([2 6 1; 3 0 4]));

julia> path = tempname();

julia> write_nfg(path, g)

julia> read_nfg(path)
3×2 NormalFormGame{2, Int64}:
 (3, 2)  (1, 3)
 (0, 6)  (4, 0)
 (2, 1)  (5, 4)

julia> read_nfg(Float64, path)
3×2 NormalFormGame{2, Float64}:
 (3.0, 2.0)  (1.0, 3.0)
 (0.0, 6.0)  (4.0, 0.0)
 (2.0, 1.0)  (5.0, 4.0)
```

A file at a URL can be read with `read_nfg(Downloads.download(url))`.
"""
read_nfg(io::IO) = _read_nfg(_parse_payoffs, io)
read_nfg(::Type{T}, io::IO) where {T<:Real} =
    _read_nfg(tokens -> _parse_payoffs(T, tokens), io)

read_nfg(path::AbstractString) = open(read_nfg, path)
read_nfg(::Type{T}, path::AbstractString) where {T<:Real} =
    open(io -> read_nfg(T, io), path)

"""
    parse_nfg([T], text)

Parse the string `text` in the Gambit .nfg format and return the game as a
`NormalFormGame`. See [`read_nfg`](@ref) for the meaning of `T` and for
reading from a stream or a file.

# Arguments

- `T::Type` : Element type of the payoffs, where `T<:Real`. If omitted, it is
  determined as in [`read_nfg`](@ref).
- `text::AbstractString` : String in the .nfg format.

# Returns

- `::NormalFormGame{N,T}` : The game described by `text`.

# Examples

```julia
julia> s = \"\"\"
       NFG 1 R "" { "1" "2" } { 3 2 }

       3 2 0 6 2 1 1 3 4 0 5 4
       \"\"\";

julia> parse_nfg(s)
3×2 NormalFormGame{2, Int64}:
 (3, 2)  (1, 3)
 (0, 6)  (4, 0)
 (2, 1)  (5, 4)
```
"""
parse_nfg(text::AbstractString) = read_nfg(IOBuffer(text))
parse_nfg(::Type{T}, text::AbstractString) where {T<:Real} =
    read_nfg(T, IOBuffer(text))

function _read_nfg(parse_payoffs, io::IO)
    tokens = SubString{String}[
        m.match for m in eachmatch(_NFG_TOKEN, read(io, String))
    ]
    (!isempty(tokens) && tokens[1] == "NFG") ||
        throw(ArgumentError("not in the .nfg format"))

    # Prologue: NFG, version, R or D, title, players, actions (the numbers of
    # actions, or the lists of their names), and an optional comment
    pos = 4
    _, pos = _read_tree(tokens, pos)  # title
    _, pos = _read_tree(tokens, pos)  # players
    actions, pos = _read_tree(tokens, pos)
    startswith(tokens[pos], '"') && (pos += 1)  # comment
    nums_actions = ntuple(length(actions)) do i
        a = actions[i]
        a isa Vector ? length(a) : parse(Int, a)
    end

    if tokens[pos] == "{"
        # Outcome version: a list of outcomes, each a name and N payoffs, then
        # the index of the outcome at each action profile, 0 meaning zero
        # payoffs
        outcomes, pos = _read_tree(tokens, pos)
        N = length(nums_actions)
        values = parse_payoffs(
            SubString{String}[x for o in outcomes for x in o[2:end]]
        )
        table = hcat(zeros(eltype(values), N), reshape(values, N, :))
        indices = [parse(Int, tok) for tok in @view tokens[pos:end]]
        payoffs = vec(table[:, indices .+ 1])
    else
        # Payoff version: the payoffs at each action profile
        payoffs = parse_payoffs(@view tokens[pos:end])
    end

    return NormalFormGame(NFGPayoffVector(nums_actions, payoffs))
end
