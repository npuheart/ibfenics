

![](logo.jpg)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](https://swanlab.cn/@SimCardiac/ideal_valve_2D/overview)
# AFSI: Automated Fluid-Structure Interaction Solver
(WIP)


## Installation

> Assume that FEniCS has be installed.

```
source .bashrc
python3 -m build && pip3 install dist/*.whl --force-reinstall
pip3 install -r requirements.txt
```

## Eulerian Solver
(WIP)

## Lagrangian Solver
(WIP)

## Eulerian-Lagrangian Interaction
(WIP)

## Documentation

```bash
cd docs && jupyter-book build . && python3 -m http.server 8000 -d _build/html
```

The address for the documentaion online is available on https://npuheart.github.io/ibfenics.



## Contents

```{tableofcontents}
```
