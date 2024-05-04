## Normal run

qm_setup()

qm_init() % initialize system
qm_bound()



## Get global variables

```matlab
global global_var_name;
```

## Get field names of a struct variable

Using function `fieldnames`:

```matlab
fieldnames(atomic)
```

## How to get type of a variable

Using function `class`:

```
>> class(atomic.m)

ans =

    'struct'
```