# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email : ibfenics@pengfeima.cn
#
# brief : 检查文件创建文件夹，生成唯一文件名. 创建日志文件. 


import os
import ibfenics.io 
from loguru import logger
from ibfenics.io import unique_filename


print(__file__)
print(os.path.basename(__file__))
print(ibfenics.io.unique_filename(os.path.basename(__file__), "tag", ".txt"))


note = "none"
file_log_name = unique_filename(os.path.basename(__file__), note, "/info.log")
logger.add(file_log_name)
logger.info(f"file_log_name = {file_log_name}")
