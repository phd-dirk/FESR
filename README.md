# FESR 
QCD analysis of e+e- annihilation into Hadrons.

## Usage
```
mkdir ./build
cd ./build
cmake ..
make
./FESR
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
