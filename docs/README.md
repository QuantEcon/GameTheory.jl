## Generate docs automatically

Run
```
$julia make.jl
```

The generated docs can be accessed in `build/index.html`.

You may modify `Structure` to arrange the structure of the website to be generated. For example:

```
Order: normal_form_game, Computing Nash equilibria, Repeated Game
Computing Nash equilibria: pure_nash, support_enumeration
Repeated Game: repeated_game_util, repeated_game
```

The first line sets the order of pages. You can add a file name if you want it to generate a single page, or the name of section you want to create. In the following you need to declare which files are included in each section. Note that you shouldn't put a comma or dot at the end of each line.

It is not necessary to declare structure for every file included in the Module. Files not mentioned in `Structure` will simply generate a single page automatically.
