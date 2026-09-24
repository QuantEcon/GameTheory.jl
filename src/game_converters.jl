# PayoffProfileMatrix #

"""
    PayoffLayout

Abstract supertype of the singleton types that specify the order in which the
payoffs of a game are listed in a flat vector; see [`PlayerMajor`](@ref) and
[`ProfileMajor`](@ref).
"""
abstract type PayoffLayout end

"""
    PlayerMajor()

Player-major order, as in the GameTracer .gam format: all the payoffs to
player 1, then all the payoffs to player 2, ..., then all the payoffs to player
N. Within each block, action profiles are ordered with player 1 varying
fastest, then player 2, ..., player N (i.e., column-major order).
"""
struct PlayerMajor <: PayoffLayout end

"""
    ProfileMajor()

Profile-major order, as in the Gambit .nfg format: the payoffs to players 1,
..., N at the first action profile, then those at the second action profile,
and so on. Action profiles are ordered with player 1 varying fastest, then
player 2, ..., player N (i.e., column-major order).
"""
struct ProfileMajor <: PayoffLayout end

"""
    PayoffProfileMatrix{N,T,TM}

Intermediate representation of the payoffs of an `N`-player game: the matrix
of size `prod(nums_actions) × N` whose `[a, i]` entry is the payoff to player
`i` at the `a`-th action profile, where action profiles are ordered with player
1 varying fastest, then player 2, ..., player N (i.e., column-major order).
This is the `payoff_profile_array` of the game with the action-profile axes
flattened into one.

The .gam and .nfg formats list the payoffs in a flat vector, in the
[`PlayerMajor`](@ref) and [`ProfileMajor`](@ref) orders respectively, which
are the column-major flattenings of this matrix and of its transpose. A
`PayoffProfileMatrix` constructed from such a vector stores a matrix backed by
the vector without copying: a `Matrix` for `PlayerMajor()`, and the
`Transpose` of a `Matrix` for `ProfileMajor()`. See [`as_vector`](@ref) for
the converse.

# Fields

- `nums_actions::NTuple{N,Int}` : Tuple of the numbers of actions, one for each
  player.
- `payoffs::TM` : Matrix of size `prod(nums_actions) × N` storing the
  payoffs, where `TM<:AbstractMatrix{T}`.
"""
struct PayoffProfileMatrix{N,T<:Real,TM<:AbstractMatrix{T}}
    nums_actions::NTuple{N,Int}
    payoffs::TM

    function PayoffProfileMatrix{N,T,TM}(
        nums_actions::NTuple{N,Int}, payoffs::TM
    ) where {N,T<:Real,TM<:AbstractMatrix{T}}
        _check_nums_actions(nums_actions)
        expected = (prod(nums_actions), N)
        size(payoffs) == expected || throw(ArgumentError(
            "payoffs size mismatch: expected $expected, got $(size(payoffs))"
        ))
        return new(nums_actions, payoffs)
    end
end

function _check_nums_actions(nums_actions)
    any(n -> n <= 0, nums_actions) &&
        throw(ArgumentError("all nums_actions must be positive"))
    return nothing
end

num_players(::PayoffProfileMatrix{N}) where {N} = N

PayoffProfileMatrix(
    nums_actions::NTuple{N,Int}, payoffs::AbstractMatrix{T}
) where {N,T<:Real} =
    PayoffProfileMatrix{N,T,typeof(payoffs)}(nums_actions, payoffs)

# The matrix backed by the vector `v` listing the payoffs in the given order
_as_matrix(v::Vector, nums_actions::NTuple{N,Int}, ::PlayerMajor) where {N} =
    reshape(v, prod(nums_actions), N)
_as_matrix(v::Vector, nums_actions::NTuple{N,Int}, ::ProfileMajor) where {N} =
    transpose(reshape(v, N, prod(nums_actions)))

"""
    PayoffProfileMatrix([T], nums_actions, payoffs, layout)

Construct a PayoffProfileMatrix (of eltype `T` if specified) from the vector
`payoffs` listing the payoffs of a game in the order `layout`, a
[`PlayerMajor`](@ref) or a [`ProfileMajor`](@ref). `payoffs` is converted to a
`Vector{T}`, which makes no copy if it already is one, and the matrix stored is
backed by that vector: a `Matrix` if `layout` is `PlayerMajor()`, and the
`Transpose` of a `Matrix` if `ProfileMajor()`.

# Examples

```julia
julia> payoffs = collect(1:12);

julia> p = GameTheory.PayoffProfileMatrix((3, 2), payoffs, GameTheory.PlayerMajor());

julia> p.payoffs
6×2 Matrix{Int64}:
 1   7
 2   8
 3   9
 4  10
 5  11
 6  12

julia> p = GameTheory.PayoffProfileMatrix((3, 2), payoffs, GameTheory.ProfileMajor());

julia> p.payoffs
6×2 transpose(::Matrix{Int64}) with eltype Int64:
  1   2
  3   4
  5   6
  7   8
  9  10
 11  12
```
"""
function PayoffProfileMatrix(
    ::Type{T}, nums_actions::NTuple{N,Int}, payoffs::AbstractVector,
    layout::PayoffLayout
) where {N,T<:Real}
    _check_nums_actions(nums_actions)
    expected = prod(nums_actions) * N
    length(payoffs) == expected || throw(ArgumentError(
        "payoffs length mismatch: expected $expected, got $(length(payoffs))"
    ))
    v = convert(Vector{T}, payoffs)
    return PayoffProfileMatrix(nums_actions, _as_matrix(v, nums_actions, layout))
end

PayoffProfileMatrix(
    nums_actions::NTuple{N,Int}, payoffs::AbstractVector{T}, layout::PayoffLayout
) where {N,T<:Real} = PayoffProfileMatrix(T, nums_actions, payoffs, layout)

# `p.payoffs` or its transpose, whichever lists the payoffs in the order
# `layout` when iterated (in column-major order); no copy is made
_in_layout(p::PayoffProfileMatrix, ::PlayerMajor) = p.payoffs
_in_layout(p::PayoffProfileMatrix, ::ProfileMajor) = transpose(p.payoffs)

"""
    as_vector(p, layout)

Return the payoffs of the PayoffProfileMatrix `p` as a `Vector` listing them in
the order `layout`, a [`PlayerMajor`](@ref) or a [`ProfileMajor`](@ref). The
vector shares memory with `p.payoffs` if that is possible without copying,
i.e., if `p.payoffs` is a `Matrix` and `layout` is `PlayerMajor()`, or if
`p.payoffs` is the `Transpose` of a `Matrix` and `layout` is `ProfileMajor()`;
otherwise it is a copy.

# Examples

```julia
julia> p = GameTheory.PayoffProfileMatrix((3, 2), collect(1:12), GameTheory.PlayerMajor());

julia> GameTheory.as_vector(p, GameTheory.PlayerMajor())'
1×12 adjoint(::Vector{Int64}) with eltype Int64:
 1  2  3  4  5  6  7  8  9  10  11  12

julia> GameTheory.as_vector(p, GameTheory.ProfileMajor())'
1×12 adjoint(::Vector{Int64}) with eltype Int64:
 1  7  2  8  3  9  4  10  5  11  6  12
```
"""
as_vector(p::PayoffProfileMatrix{N,T}, layout::PayoffLayout) where {N,T} =
    convert(Vector{T}, vec(_in_layout(p, layout)))


# Forward: (i, i+1, ..., N, 1, ..., i-1)
@inline _perm_fwd(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(i + k - 1, N), Val(N))

# Backward: inverse rotation
@inline _perm_back(::Val{N}, ::Val{i}) where {N,i} =
    ntuple(k -> mod1(k - i + 1, N), Val(N))

# The payoffs to player i as an N-dim array indexed by the action profile
# (a_1, ..., a_N); shares memory with `p.payoffs`
_player_block(p::PayoffProfileMatrix{N,T,Matrix{T}}, i::Int) where {N,T} =
    selectdim(reshape(p.payoffs, (p.nums_actions..., N)), N+1, i)
_player_block(p::PayoffProfileMatrix{N}, i::Int) where {N} =
    reshape(view(p.payoffs, :, i), p.nums_actions)


"""
    PayoffProfileMatrix([T], g)

Construct a PayoffProfileMatrix (of eltype `T` if specified) from a
NormalFormGame `g`. The matrix stored is a `Matrix{T}`.

# Examples

```julia
julia> player1 = Player([1 4; 2 5; 3 6]);

julia> player2 = Player([7 8 9; 10 11 12]);

julia> g = NormalFormGame(player1, player2)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)

julia> p = GameTheory.PayoffProfileMatrix(g);

julia> p.payoffs
6×2 Matrix{Int64}:
 1   7
 2   8
 3   9
 4  10
 5  11
 6  12
```
"""
function PayoffProfileMatrix(::Type{T}, g::NormalFormGame{N}) where {N,T<:Real}
    nums_actions = g.nums_actions
    payoffs = Matrix{T}(undef, prod(nums_actions), N)
    p = PayoffProfileMatrix(nums_actions, payoffs)

    ntuple(Val(N)) do i
        copyto!(
            _player_block(p, i),
            PermutedDimsArray(g.players[i].payoff_array, _perm_back(Val(N), Val(i)))
        )
        nothing
    end

    return p
end

PayoffProfileMatrix(g::NormalFormGame{N,T}) where {N,T<:Real} =
    PayoffProfileMatrix(T, g)


"""
    NormalFormGame([T], p)

Construct a NormalFormGame (of eltype `T` if specified) from a
PayoffProfileMatrix `p`. The payoffs are copied.

# Examples

```julia
julia> payoffs = collect(1:12);

julia> p = GameTheory.PayoffProfileMatrix((3, 2), payoffs, GameTheory.PlayerMajor());

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 7)  (4, 10)
 (2, 8)  (5, 11)
 (3, 9)  (6, 12)

julia> p = GameTheory.PayoffProfileMatrix((3, 2), payoffs, GameTheory.ProfileMajor());

julia> NormalFormGame(p)
3×2 NormalFormGame{2, Int64}:
 (1, 2)  (7, 8)
 (3, 4)  (9, 10)
 (5, 6)  (11, 12)
```
"""
function NormalFormGame(::Type{T}, p::PayoffProfileMatrix{N}) where {N,T<:Real}
    players = ntuple(Val(N)) do i
        Player(
            T,
            PermutedDimsArray(_player_block(p, i), _perm_fwd(Val(N), Val(i)))
        )
    end

    return NormalFormGame{N,T}(players, p.nums_actions)
end

NormalFormGame(p::PayoffProfileMatrix{N,T}) where {N,T<:Real} =
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
[`PlayerMajor`](@ref) for the ordering of the payoffs in the format, and
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

    # Payoffs, in .gam order; the length is checked by PayoffProfileMatrix
    payoffs = parse_payoffs(@view tokens[N+2:end])

    return NormalFormGame(
        PayoffProfileMatrix(nums_actions, payoffs, PlayerMajor())
    )
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
- `g::Union{NormalFormGame,PayoffProfileMatrix}` : Game to write.

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
function write_gam(io::IO, p::PayoffProfileMatrix{N,T}) where {N,T<:_PayoffNumber}
    # `print` would round floats if the caller's context has `:compact => true`
    io = IOContext(io, :compact => false)
    print(io, N, '\n')
    join(io, p.nums_actions, ' ')
    print(io, "\n\n")  # blank line between the header and the payoffs
    # Iterated in column-major order, i.e., in player-major order; no copy
    for (k, x) in enumerate(_in_layout(p, PlayerMajor()))
        k > 1 && print(io, ' ')
        _print_payoff(io, x)
    end
    print(io, '\n')
    return nothing
end

write_gam(io::IO, g::NormalFormGame{N,T}) where {N,T<:_PayoffNumber} =
    write_gam(io, PayoffProfileMatrix(g))

# Same bound on the element type as the methods for `io`, so that an
# unsupported game is rejected before the file is opened
write_gam(
    path::AbstractString, g::Union{NormalFormGame{N,T},PayoffProfileMatrix{N,T}}
) where {N,T<:_PayoffNumber} = open(io -> write_gam(io, g), path, "w")

"""
    gam_string(g)

Return the GameTracer .gam representation of the game `g` as a string. See
[`write_gam`](@ref) for the requirement on the element type of `g`.

# Arguments

- `g::Union{NormalFormGame,PayoffProfileMatrix}` : Game to write.

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
gam_string(g::Union{NormalFormGame,PayoffProfileMatrix}) = sprint(write_gam, g)
