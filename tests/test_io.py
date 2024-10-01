# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : 检查文件创建文件夹，生成唯一文件名


import os
import ibfe.io 

print(__file__)
print(os.path.basename(__file__))
print(ibfe.io.unique_filename(os.path.basename(__file__), "tag", ".txt"))
