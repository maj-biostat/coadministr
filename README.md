# coadministr

Trial simulation package.
Depends on [coadministr.stanc](https://github.com/maj-biostat/coadministr.stanc).

## Build

```
git clone https://github.com/maj-biostat/coadministr
R CMD build coadministr
```

## Install

```
git clone https://github.com/maj-biostat/coadministr
R CMD INSTALL --preclean coadministr
```

or after build:

```
R CMD INSTALL --preclean coadministr_1.0.tar.gz
```

## Environment

Brief environment information:

+ Ubuntu 20.04.2 LTS

Tree:

```
coadministr
├── DESCRIPTION
├── man
│   ├── coadministr-package.Rd
│   ├── get_data.Rd
│   ├── pr_adverse_evt.Rd
│   └── trial.Rd
├── NAMESPACE
├── R
│   ├── data.R
│   ├── trial.R
│   ├── utils.R
│   └── zzz.R
└── README.md
```



