# Introduction

This is an implementation of the KV-CBO method proposed by Massimo Fornasier, Hui Huang, Lorenzo Pareschi and Philippe Sünnen.

This implementation is open source and is distributed under the MIT license.

## License

Copyright (c) 2020, Massimo Fornasier, Hui Huang, Loranzo Pareschi, Philippe Sünnen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Third party code

The folder ext contains code from thrid party developers. The precise citations are:

@online{SphericalDistributionsRand,
    author = {Yu-Hui Chen},
    title = {SphericalDistributionsRand,
    url = {https://github.com/yuhuichen1015/SphericalDistributionsRand},
} 

@online{RotMatrix,
    author = {Jan Simon},
    title = {Rotation Matrix},
    url = {https://de.mathworks.com/matlabcentral/fileexchange/66446-rotation-matrix},
}

@online{HyperSphere,
    author = {Gianluca Dorini},
    title = {HyperSphere},
    url = {https://de.mathworks.com/matlabcentral/fileexchange/5397-hypersphere},
}

The file linesearch.m relies extensively on third party developers. The file KVCBO.m contains third party contributions (see section "Statistics" in the code). The precise citation is:

@Article{manopt,
    author = {Boumal, N. and Mishra, B. and Absil, P.-A. and Sepulchre, R.},
    journal = {Journal of Machine Learning Research},
    title = {{M}anopt, a {M}atlab Toolbox for Optimization on Manifolds},
    year = {2014},
    number = {42},
    pages = {1455--1459},
    volume = {15},
    url = {https://www.manopt.org},
} 
