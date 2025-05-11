
<!-- 代码一旦写好了，就不应该再打开。所有变量，要么用配置文件，要么用命令行。如果你需要打开代码进行修改，那你的代码就不合格。 -->
<!-- 如果想走学术这条路，就得启动这个大循环：发论文->申项目->升职称->带学生，而且循环越快越好，最好实现自发运转。  -->



# IB-FEniCS

(WIP)

## 安装

```
source .bashrc
python3 -m build && pip3 install dist/*.whl --force-reinstall
pip3 install -r requirements.txt
```

## Test Eulerian-Lagrangian Interaction



## 可能存在问题
1. 函数 iterate_grid_3D 中的 index_type 可能出现负值。


# Welcome to your Jupyter Book

This is a small sample book to give you a feel for how book content is
structured.
It shows off a few of the major file types, as well as some sample content.
It does not go in-depth into any particular topic - check out [the Jupyter Book documentation](https://jupyterbook.org) for more information.

Check out the content pages bundled with this sample book to see more.

```{tableofcontents}
```

This is an inline math expression: $E = mc^2$

And a block math expression:

$$
\int_a^b f(x) \, dx
$$
