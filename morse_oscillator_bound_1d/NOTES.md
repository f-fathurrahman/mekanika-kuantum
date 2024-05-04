## TODO

How to initialize Hamiltonian matrix?

How to initialize the potential?


## Normal run

This is the usual sequence call:
```matlab
qm_setup()
qm_init() % initialize system
qm_bound()
```

## Get global variables

```matlab
global global_var_name;
```

## Working with data types

Get field names of a struct variable, using function `fieldnames`:

```matlab
fieldnames(atomic)
```

How to get type of a variable, using function `class`:

```
>> class(atomic.m)

ans =

    'struct'
```

Examples:
```
>> class(hamilt)
ans = struct
>> class(hamilt.pot)
ans = cell
>> class(hamilt.pot{1})
ans = pot.morse
```