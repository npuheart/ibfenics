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
from dolfin import *
from mshr import *
import json

output_path = "data/large/"


def check_path(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"{path} 文件夹已创建。")
    else:
        print(f"{path} 文件夹已存在。")


def unique_filename(current_file_name, tag, extension):
    check_path(output_path)
    note = os.path.splitext(current_file_name)[0]
    file_id = (
        f"{output_path}{note}/{tag}/"
        + datetime.now().strftime("%Y%m%d-%H%M%S")
        + extension
    )
    return file_id


def create_xdmf_file(mpi_comm, filename):
    file = XDMFFile(mpi_comm, filename)
    file.parameters["rewrite_function_mesh"] = False
    file.parameters["functions_share_mesh"] = True
    file.parameters["flush_output"] = True
    return file


def send_to_xiaomi(desp):
    # pip install serverchan-sdk
    sdk_key = "sctp1508t8kbmtrb5kat2nsazwhmqy2"
    from serverchan_sdk import sc_send

    response = sc_send(sdk_key, "ibfenics log", desp, {"tags": "ibfenics"})

    return response


# 指定要保存的文件名和表单名称
def write_excel(volume_list, excel_file1, sheet_name="v"):
    df = pd.DataFrame(volume_list)
    df.to_excel(excel_file1, sheet_name=sheet_name, index=False)


def write_paramters(filename, **params):
    params = params or {}
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(params, f, ensure_ascii=False, indent=4)


# 1. 时间步长为 total_time / total_steps
# 2. total_steps 太大了怎么办？ 用 fps (frame per second) 控制输出的步数，fps*total_time
# 3. total_time 特别小怎么办？ 忽略fps，直接按照 total_steps 输出
class TimeManager:
    def __init__(self, total_time, total_steps, fps=100):
        self.total_time = total_time
        self.total_steps = total_steps
        self.fps = fps
        self.time_per_step = total_time / total_steps
        fps_m_time = max(fps * total_time, 1)
        self.step_interval = max(int(total_steps / fps_m_time), 1)

    def should_output(self, current_step):
        # 判断当前步是否是输出步
        if current_step > self.total_steps:
            raise ValueError("当前步数超过总步数。")
        if current_step % int(self.step_interval) == 0:
            return True
        # 额外输出的步数：0,1,total_steps-1,total_steps, 保证一开始和最后都有输出。
        if current_step == 1:
            return True
        if current_step == self.total_steps:
            return True
        if current_step == self.total_steps - 1:
            return True
        return False


def write_mesh(_mesh, _domains, boundary, mesh_size, mesh_name, mesh_path):
    # output mesh
    mesh_file = XDMFFile(mesh_path + mesh_name + "_" + str(mesh_size) + ".xdmf")
    mesh_file.write(_mesh)
    mesh_file.close()

    # output domains
    domains_file = XDMFFile(
        mesh_path + mesh_name + "_" + str(mesh_size) + "_domains.xdmf"
    )
    domains_file.write(_domains)
    domains_file.close()

    # output domains xml
    File(mesh_path + mesh_name + "_" + str(mesh_size) + "_domains.xml") << _domains

    # output boundaries xml
    domains_file = XDMFFile(
        mesh_path + mesh_name + "_" + str(mesh_size) + "_boundaries.xdmf"
    )
    domains_file.write(boundary)
    domains_file.close()
    File(mesh_path + mesh_name + "_" + str(mesh_size) + "_boundaries.xml") << boundary


def write_bg_mesh(L, H, mesh_size, mesh_name, mesh_path):
    _domain = Rectangle(Point(0, 0), Point(L, H))
    _mesh = generate_mesh(_domain, mesh_size)

    # output boundaries
    class Boundary_1(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] < DOLFIN_EPS_LARGE and on_boundary

    class Boundary_4(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] > L * (1.0 - DOLFIN_EPS_LARGE) and on_boundary

    class Boundary_2(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] < DOLFIN_EPS_LARGE and on_boundary

    class Boundary_3(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] > H * (1.0 - DOLFIN_EPS_LARGE) and on_boundary

    boundary = MeshFunction("size_t", _mesh, 1)
    Boundary_1().mark(boundary, 1)
    Boundary_2().mark(boundary, 2)
    Boundary_3().mark(boundary, 3)
    Boundary_4().mark(boundary, 4)

    # output mesh
    mesh_file = XDMFFile(mesh_path + mesh_name + "_bg_" + str(mesh_size) + ".xdmf")
    mesh_file.write(_mesh)
    mesh_file.close()

    domains_file = XDMFFile(
        mesh_path + mesh_name + "_bg_" + str(mesh_size) + "_boundaries.xdmf"
    )
    domains_file.write(boundary)
    domains_file.close()
    File(
        mesh_path + mesh_name + "_bg_" + str(mesh_size) + "_boundaries.xml"
    ) << boundary


__all__ = [
    "unique_filename",
    "create_xdmf_file",
    "write_excel",
    "write_bg_mesh",
    "write_mesh",
]
