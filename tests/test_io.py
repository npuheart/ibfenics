# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics1 (https://github.com/npuheart/ibfenics1)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : mapengfei@mail.nwpu.edu.cn
#
# brief : 检查文件创建文件夹，生成唯一文件名. 创建日志文件.


import os
import ibfenics1.io
from loguru import logger
from ibfenics1.io import unique_filename


print(__file__)
print(os.path.basename(__file__))
print(ibfenics1.io.unique_filename(os.path.basename(__file__), "tag", ".txt"))


note = "none"
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)
logger.info(f"file_log_name = {file_log_name}")
