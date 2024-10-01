# Copyright (C) 2024 Pengfei Ma
#
# This file is part of ibfenics (https://github.com/npuheart/ibfenics)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# email: ibfenics@pengfeima.cn
# brief: 如果想要将数据存储到大文件夹, 可以创建一个软链接:
#           mkdir /mnt/large0/data0/ibfenics
#           ln -s /mnt/large0/data0/ibfenics data/large

import os
import pandas as pd
from datetime import datetime
from dolfin import XDMFFile

output_path = "data/large/"

def check_path(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"{path} 文件夹已创建。")
    else:
        print(f"{path} 文件夹已存在。")

def unique_filename(current_file_name, tag, extension):
    check_path(output_path)
    note =  os.path.splitext(current_file_name)[0]  # 去掉扩展名
    file_id = f"{output_path}{note}/{tag}/" + datetime.now().strftime('%Y%m%d-%H%M') + extension
    return file_id

def create_xdmf_file(mpi_comm, filename):
    file = XDMFFile(mpi_comm, filename)
    file.parameters['rewrite_function_mesh'] = False
    file.parameters["functions_share_mesh"] = True
    file.parameters["flush_output"] = True
    return file


# 指定要保存的文件名和表单名称
def write_excel(volume_list,excel_file1,sheet_name = 'v'):
    df = pd.DataFrame(volume_list)
    df.to_excel(excel_file1, sheet_name=sheet_name, index=False)


__all__ = ['unique_filename', 'create_xdmf_file', 'write_excel']

