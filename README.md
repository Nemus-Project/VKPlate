# VKPlate

Implementation of a Föppl–von Kármán (VK) plate.

| ![16-mode modal plate doing with initial excitation](./img/modeplate.gif) |
| :-----------------------------------------------------------------------: |
|             16-mode modal plate doing with initial excitation             |

- [VKPlate](#vkplate)
  - [About](#about)
  - [Getting started](#getting-started)
    - [Cloning this Repository](#cloning-this-repository)
    - [Directory setup](#directory-setup)
    - [Using main.m](#using-mainm)
    - [Parameter File Format](#parameter-file-format)
  - [Further Goals](#further-goals)
  - [References](#references)
  - [Troubleshooting](#troubleshooting)



## About

Simulation of a von-Karman plate written in modal formalism, using scalar auxiliary variable (SAV) method.

## Getting started

### Cloning this Repository

This project uses the [`magpie-matlab`](https://github.com/Nemus-Project/magpie-matlab) project as a submodule. To correctly clone this repository yo will need to use the command

```sh
git clone --recurse-submodules https://github.com/Nemus-Project/VKPlate
```

Alternatively, clone via the [GitHub Desktop](https://github.com/apps/desktop) application or [CLI Client](https://cli.github.com).

**NOTE:** Downloading the zip will skip the download of the submodules. 

### Directory setup

There are some assets that need to be generated or downloaded before the project is useable.

The `main.m` script needs loads a `.mat` file containing a specific set of variables. 

You can either 

1. generate it with `genparams.m` 
2. download the [`param.zip`](https://github.com/Nemus-Project/VKPlate/releases/download/0.2.0/param.zip)

The contents of `param.zip` should be put in a directory named `param/` at the top level of the repository. The data is quite weighty (~2GB).

The `genparams.m` will create the `param/` directory by default.

At the end the repository directory should look like:

```tree
VKPlate
├── README.md
├── genparams.m
├── magpie
├── main.m
├── img
│   └── ...
├── param
│   └── ...
├── private
│   ├── DxBuild.m
│   ├── DxxBuild.m
│   ├── DxyBuild.m
│   ├── DyBuild.m
│   ├── DyyBuild.m
│   ├── eigensign.m
│   ├── eigenMAC.m
│   ├── trapzIntcalc.m
│   └── vkOperator.m
└── test
    ├── Hcalc.m
    ├── ...
    └── SplitDuffing.m
```


### Using main.m

Once the `param/` is populates, run the `main.m`.

### Parameter File Format

To avoid needless and time consuming recalculation, the central `main.m` script loads 
a `.mat` with all the require variables.

| name         | description |
| ------------ | ----------- |
| `rho`        |             |
| `E`          |             |
| `nu`         |             |
| `Lz`         |             |
| `Lx`         |             |
| `Ly`         |             |
| `Nmodes`     |             |
| `Phi`        |             |
| `Om`         |             |
| `Psi`        |             |
| `Om2`        |             |
| `Nx`         |             |
| `Ny`         |             |
| `h`          |             |
| `X`          |             |
| `Y`          |             |
| `zetafourth` |             |
| `BCsPhi`     |             |
| `BCsPsi`     |             |
| `Hv`         |             |

Use the `genparams.m` script to create the `.mat` files or download from the latest release (See [Directory setup](#directory-setup))

## Further Goals

Final checks on SAV, different boundary condition, plate under tension

## References

TBA

## Troubleshooting

TBA