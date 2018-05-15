# FESR 
QCD analysis of e+e- annihilation into Hadrons.

## Install
```
mkdir ./build
cd ./build
cmake ..
make
```

## Usage for single fit
1. Edit the `./configuration.json`
2. run `./FESR`

## Usage for multiple fits
1. Add arbitray number of `configuration.json` in `./configuration` folder
2. run `./run.sh`

## Configuration.json example
```
{
    "parameters": {
        "nc": 3,
        "nf": 3,
        "order": 5,
        "alphaLoops": 5,
        "s0Set": [3.1572314596, 3.0, 2.8, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0],
        "weight": 1,
        "RVANormalization": 0.99743669
    },
    "constants": {
        "sTau": 3.1572314596,
        "be": 17.815,
        "dBe": 0.023
    },
    "variables": {
        "astau": {
            "fixed": false,
            "value": 0.3,
            "stepSize": 0.002
        },
        "aGGInv": {
            "fixed": true,
            "value": 0.021,
            "stepSize": 0.01
        },
        "rhoVpA": {
            "fixed": false,
            "value": -0.30949,
            "stepSize": 0.1
        },
        "c8VpA": {
            "fixed": false,
            "value": -0.030869,
            "stepSize": 0.3
        }
    },
    "adler": {
        "D0": true,
        "D2": false,
        "D4": true,
        "D68": true,
        "PionPole": true
    }
}

```


## Testing
`cmake Dtest=ON` tells cmake to include testing sourcefiles and generate `runUnitTests` executable. 
```
cd ./build
cmake -Dtest=ON ..
make
./runUnitTests
```

## Developer
### C++ autocomplete 
Generate **.clang_complete** file in the build folder.
```
cd ./build
CXX="cc_args.py clang++" cmake ..
make
cp .clang_complete ../.clang_complete
```

## Dependencies 
* [json](https://github.com/nlohmann/json)
* [gtest](https://github.com/google/googletest)
* [ROOT](https://root.cern.ch/)
* [BOOST](https://www.boost.org/)
* [GSL](https://www.gnu.org/software/gsl/doc/html/index.html)
