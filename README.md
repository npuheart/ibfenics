
<!-- 代码一旦写好了，就不应该再打开。所有变量，要么用配置文件，要么用命令行。如果你需要打开代码进行修改，那你的代码就不合格。 -->
<!-- 如果想走学术这条路，就得启动这个大循环：发论文->申项目->升职称->带学生，而且循环越快越好，最好实现自发运转。  -->
## 依赖
- FEniCS 2019.1
- pybind11, cmake: 暴露 C++ 程序 Python 接口
- github submodule 网格文件(WIP)



## 在用户目录下安装此库
```bash
python3 setup.py install --user
```

## 测试
```bash
python3 tests/test_duality.py
```

## 通过以下环境测试
- gcc 9.4.0
- cmake 3.27.8
- Python 3.8.10 
- pybind11 2.4.3




```bash
pip3 show serverchan-sdk
```

## 虚拟环境
python3 -m venv myenv
source myenv/bin/activate
pip3 install -r requirements.txt



## 可能存在问题
1. 函数 iterate_grid_3D 中的 index_type 可能出现负值。