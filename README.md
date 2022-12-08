# Regular Expression Partitioner

Written by Mikel Mcdaniel

## What is it?

regex_partitioner is a Python3 library and script designed to solve a useless academic problem:
take all the strings that a regular expression accepts and split them up into equally sized groups.

### Simple Example

The regular expression `[abc]*` accepts any string made of the letters a, b, and c.
There are infinitely many strings like that: "", "a", "b", "c", "aa", "ab", ...

Even with infinitely many strings, 1/3 of them would come before the string "b" in sorted order, 1/3 would be between "b" and "c", and the last 1/3 would come after "c".
The strings "b" and "c" partition the strings that the regular expression `[abc]*` accepts into 3 equally sized partitions.

That's exactly the kind of problem this code solves... mostly:

```commandline
$ python3 regex_partitioner.py --partitions 3 "[abc]*"
The regular expression '[abc]*' accepts 88573 strings with length <= 10.

These strings can be evenly split, within a 0.100000% tolerance by the strings:
  '[abc]*' accepts 59049 (66.66704300407574%) of strings after 'accccccccc' (or 'accccccccc').
  '[abc]*' accepts 29525 (33.334086008151466%) of strings after 'bccccccccc' (or 'bccccccccc').
```

### Limitations

By default, regex_partitioner will not give an exact answer and will instead give an answer that is within some `--tolerance_ratio`.
If you want an exact answer, you can use a tolerance ratio of 0.

This code does not partition the set of *all* strings that the regular expression can accept;
It partitions the set of all strings less than `--max_str_len`. However, you can set the maximum string length to be as large as you want.

Example:

```commandline
$ python3 regex_partitioner.py --partitions 3 --tolerance_ratio 0 --max_str_len 50 "[abc]*"
The regular expression '[abc]*' accepts 1076846981537778883155373 strings with length <= 50.

These strings can be evenly split, within a 0.000000% tolerance by the strings:
  '[abc]*' accepts 717897987691852588770249 (66.66666666666667%) of strings after 'accccccccccccccccccccccccccccccccccccccccccccccccc' (or 'accccccccccccccccccccccccccccccccccccccccccccccccc').
  '[abc]*' accepts 358948993845926294385125 (33.333333333333336%) of strings after 'bccccccccccccccccccccccccccccccccccccccccccccccccc' (or 'bccccccccccccccccccccccccccccccccccccccccccccccccc').
```

### Other Examples

#### Lower and Upper Bounds

```commandline
$ python3 regex_partitioner.py --lower_bound ab --upper_bound bc "[abc]*"
The regular expression '[abc]*' accepts 88573 strings with length <= 10.
Between 'ab' and 'bc', '[abc]*' accepts 39365 with len <= 10.
...
```

#### Uhhh

How many floating point representations are there in 64 characters or less? Alot.

```commandline
$ python3 regex_partitioner.py --max_str_len 64 "[+-]?[1-9][0-9]*\.?[0-9]*([eE][+-]?[0-9]+)?"
The regular expression '[+-]?[1-9][0-9]*\\.?[0-9]*([eE][+-]?[0-9]+)?' accepts 130565925925925925925925925925925925925925925925925925925925925936 strings with length <= 64.

These strings can be evenly split, within a 0.100000% tolerance by the strings:
  '[+-]?[1-9][0-9]*\\.?[0-9]*([eE][+-]?[0-9]+)?' accepts 97923881810699588477366255144032921810699588477366255144032921807 (74.99956908072235%) of strings after '+549.' (or '+549-.E26E6e6--7E05261E4e5e.9374.222-EE477ee66.35+8E9160844.1335').
  '[+-]?[1-9][0-9]*\\.?[0-9]*([eE][+-]?[0-9]+)?' accepts 65282962962962962962962962962962962962962962962962962962962962968 (50.0%) of strings after '-1' (or '+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+e+Ee').
  '[+-]?[1-9][0-9]*\\.?[0-9]*([eE][+-]?[0-9]+)?' accepts 32640918847736625514403292181069958847736625514403292181069958839 (24.999569080722353%) of strings after '-549.' (or '-549-.e539+E6712E-E729+3595299EE761.e4E5eE544e.33243E0E.46e0-317').
```

## How do I use it?

As of 2022-11, this code is Python 3.10 compatible.
It's probably compatible with any Python 3.x version, but I haven't test it.

If you want to use this inside some of your own code, then have a look at regex_partitioner.main and the unit tests in regex_partitioner_test.

If you want to use it on the command line, then look at the examples in this README.

## Why does it exist?

It's a fun side project. I like regular expressions, finite automata, and puzzles.

## How does it work?

Ask me!