# VKPlate

Implementation of a Föppl–von Kármán (VK) plate.

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
├── img
│   └── modeplate.gif
├── param
│   ├── Test_NL_Fullclamp_1.mat
│   ├── Test_NL_Fullclamp_2.mat
│   ├── Test_NL_Fullclamp_3.mat
│   ├── Test_NL_Fullclamp_4.mat
│   ├── Test_NL_Fullclamp_5.mat
│   └── Test_NL_Fullclamp_6.mat
├── src
│   ├── DxBuild.m
│   ├── DxxBuild.m
│   ├── DxyBuild.m
│   ├── DyBuild.m
│   ├── DyyBuild.m
│   ├── eigenMAC.m
│   ├── eigensign.m
│   ├── genparams.m
│   ├── magpie(submodule)
│   ├── main.m
│   ├── trapzIntcalc.m
│   └── vkOperator.m
└── test    
    └── *.m
```


### Using main.m

Once the `param/` is populates, run the `src/main.m`.

### Example Output

| ![16-mode modal plate doing with initial excitation](./img/modeplate.gif) |
| :-----------------------------------------------------------------------: |
|             16-mode modal plate doing with initial excitation             |

## Further Goals

Final checks on SAV, different boundary condition, plate under tension

## References

TBA

## Troubleshooting

TBA