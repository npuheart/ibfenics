

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
