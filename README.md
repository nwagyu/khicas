# KhiCAS for NumWorks

[![Build](https://github.com/nwagyu/khicas/actions/workflows/build.yml/badge.svg)](https://github.com/nwagyu/khicas/actions/workflows/build.yml)

This apps lets you run [Xcas](https://www-fourier.ujf-grenoble.fr/~parisse/giac.html) on your [NumWorks calculator](https://www.numworks.com).

## Install the app

Installing is rather easy:
1. Download the latest `khicas.nwa` file from the [Releases](https://github.com/nwagyu/khicas/releases) page
2. Get a [Nwagra](https://www.nwagyu.com/pages/extended-memory/) pill
3. Head to [my.numworks.com/apps](https://my.numworks.com/apps) to send the `nwa` file on your calculator

## Dependencies

This programs uses four libraries:

|Library|Version|
|-|-|
|[gmp](https://gmplib.org/)|6.2.1|
|[mpfr](https://www.mpfr.org/)|4.1.0|
|[mpfi](http://perso.ens-lyon.fr/nathalie.revol/software.html)|1.5.4|
|[giac](https://www-fourier.ujf-grenoble.fr/~parisse/giac.html)|1.9.0-21|

## Build the app

To build this sample app, you will need to install the [embedded ARM toolchain](https://developer.arm.com/Tools%20and%20Software/GNU%20Toolchain) and [nwlink](https://www.npmjs.com/package/nwlink).

```shell
brew install numworks/tap/arm-none-eabi-gcc node # Or equivalent on your OS
npm install -g nwlink
make clean && make build
```
