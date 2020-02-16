

# GUI version

# Command line version

## Arguments

```bash
swEOS -D1 -VT -R0/2/1000 -P316 -X0.032
```

### Required arguments

* -D[0|1|2|3]: calculation model (**D**imension) selection, 0~3 represent scatter point, 1D, 2D, 3D calculation, respectively

* -V[PTX|PHX|T|P|X|H|PT|PX|TX|PH|HX]: **V**ariables selection. **Note** that different `-D` parameter correspond to different options for `-V`.

    1) -D0: [PTX|PHX]
    2) -D1: [T|P|X|H]
    3) -D2: [PT|PX|TX|PH|HX]
    4) -D3: [PTX|PHX]

* -P: fixed pressure value

* -T: fixed temperature value

* -X: fixed salinity value

* -H: fixed enthalpy value

* -R: only for 1~3 dimension. Corresponding to `-V` parameter, e.g. `-VPT` and `-R10/1/600/0/1/1000` means pressure in range of [10, 600] bar and interval is 1 bar, temperature in range of [0, 1000] deg. C and interval is 1 deg. C.

* -G: input file for multiple points calculation

* -O: output file for 1D, 2D, and 3D calculation

### Optional arguments

* -t: number of threads
